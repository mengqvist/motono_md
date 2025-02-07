import argparse
import math
import numpy as np
import os
import re
import sys
import traceback

from argparse import RawTextHelpFormatter
from contact_profile_calculator import _parse_region, _natural_keys
from log import get_logger

SIMILARITY_TANIMOTO = 'TANIMOTO'
SIMILARITY_EUCLIDEAN = 'EUCLIDEAN'


def drawer_main(contact_profile_diff_file: str,
         similarity_table_file: str,
         trajectory: str,
         exclude_region: str = '',
         similarity_mode: str = SIMILARITY_TANIMOTO,
         output_dir: str = None):
    """Append similarity information to PDB

    Args:
        contact_profile_diff_file (str): Contact retention rate delta file
        similarity_table_file (str): similarity table
        trajectory (str): multi-PDB / PDB
        no_end (bool): Exclusion flag for terminal residues. Defaults to False
        similarity_mode (str): Similarity referenced in the similarity table. 'TANIMOTO' or 'EUCLIDEAN'.
        output_dir (str, optional): Output destination. Defaults to None.
    """
    log_file_path = os.path.join(output_dir, "similarity_drawer.log")
    lg = get_logger(name=__name__, log_file_path=log_file_path, is_stdout=True)
    lg.info('START similarity_drawer')

    lg.info(
        'args: \n'
        'contact_profile_diff_file: {diff}\n'
        'similarity_table_file: {table}\n'
        'trajectory: {pdb}\n'
        'exclude_region: {exclude}\n'
        'similarity_mode: {s}'.format(
            diff=contact_profile_diff_file, table=similarity_table_file,
            pdb=trajectory, exclude=exclude_region, s=similarity_mode)
    )

    try:
        _validate(exclude_region)
        pdb_prefix = os.path.splitext(os.path.basename(trajectory))[0]

        # Residue regions excluded from processing
        exclude_region = _parse_region(exclude_region)
        # list of excluded residue IDs
        exclude_residue_list = list()

        # Get a dictionary that associates the residue number with the similarity of the contact profile
        similarity_table_dict = dict()
        with open(similarity_table_file, 'r') as table_rf:
            for row_line in table_rf.readlines():
                if(row_line.startswith("residue_id")):
                    continue
                cols = row_line.strip().split(',')
                # Residue ID format：[Chain][residue number][icode]
                residue_id = cols[0].strip()

                if(_is_include_region(residue_id, exclude_region)):
                    # Add excluded residues
                    exclude_residue_list.append(residue_id)

                elif(similarity_mode == SIMILARITY_TANIMOTO):
                    similarity_tanimoto = float(cols[1])
                    similarity_table_dict[residue_id] = similarity_tanimoto

                elif(similarity_mode == SIMILARITY_EUCLIDEAN):
                    similarity_euclidean = float(cols[2])
                    similarity_table_dict[residue_id] = similarity_euclidean

        # Sort by lowest similarity
        if(similarity_mode == SIMILARITY_TANIMOTO):
            # TANIMOTO coefficient: The less the two contact profiles match, the smaller the value
            # -> sort in ascending order
            similarity_table_dict = sorted(
                similarity_table_dict.items(), key=lambda x: x[1])

        elif(similarity_mode == SIMILARITY_EUCLIDEAN):
            # Euclidean distance: the greater the mismatch between the two contact profiles, the greater the value
            # -> sort in ascending order
            similarity_table_dict = sorted(
                similarity_table_dict.items(), key=lambda x: x[1], reverse=True)

        lg.info('exclude_residue_list = {}'.format(exclude_residue_list))

        dict_length = len(similarity_table_dict)

        def _sub_main(target_resi_list: list,
                      contact_profile_diff_file: str,
                      trajectory: str,
                      output_pdb_prefix: str,
                      output_dir: str,
                      exclude_list: list,
                      logger=None):
            """contact profile diff Extract interacting residues for each threshold value and output to PDB.

            Args:
                target_resi_list (list): Interesting residue list
                contact_profile_diff_file (str): Contact retention rate delta file
                trajectory (str): multi-PDB / PDB
                output_pdb_prefix (str): Output file prefix
                exclude_list (list): exclusion list.
                logger (logging.Logger, optional): Logger. Defaults to None.
            """
            for frequency_cutoff in [(-0.1, None), (None, 0.1), (-0.1, 0.1)]:
                cutoff_upper = None
                cutoff_lower = None
                file_prefix_atom = None
                pseudo_name = ' CA'

                if(frequency_cutoff[1] is None):
                    cutoff_lower = frequency_cutoff[0]
                    file_prefix_atom = 'neg'

                elif(frequency_cutoff[0] is None):
                    cutoff_upper = frequency_cutoff[1]
                    file_prefix_atom = 'pos'

                else:
                    cutoff_lower = frequency_cutoff[0]
                    cutoff_upper = frequency_cutoff[1]
                    file_prefix_atom = 'st'

                target_residue_pairs_dict = get_partner_residue(
                    contact_profile_diff_file=contact_profile_diff_file,
                    target_residue=target_resi_list,
                    frequency_cutoff_upper=cutoff_upper,
                    frequency_cutoff_lower=cutoff_lower,
                    exclude_list=exclude_list,
                    logger=logger)

                output_pdb = os.path.join(output_dir,
                                          "{}_{}.pdb".format(output_pdb_prefix, file_prefix_atom))

                save_pdb(input_pdb=trajectory,
                         output_pdb=output_pdb,
                         similarity_table_dict=dict(similarity_table_dict),
                         target_as_psude_residue_list=target_resi_list,
                         target_pair_dict=target_residue_pairs_dict,
                         pseudo_name=pseudo_name,
                         logger=logger)

        # similarity Top X% residues when sorting the values
        for x_per in [0.05, 0.1, 0.2, 0.5, 1]:
            ref_length = math.floor(dict_length * x_per)
            target_resi_list = [
                data[0] for data in similarity_table_dict[:ref_length]]

            lg.info("START get similarity data")
            lg.info("similarity condition: Top {} percent:({} cases/{} cases)".format(
                x_per * 100, ref_length, dict_length))
            lg.debug("target_resi_list: \n{}".format(target_resi_list))

            output_pdb_prefix = "{fix}_top{x}pct".format(fix=pdb_prefix, x=int(x_per * 100))

            _sub_main(target_resi_list=target_resi_list,
                      contact_profile_diff_file=contact_profile_diff_file,
                      trajectory=trajectory,
                      output_pdb_prefix=output_pdb_prefix,
                      output_dir=output_dir,
                      exclude_list=exclude_residue_list,
                      logger=lg)

        # Similarity Top X residues when sorting the values
        for x in [10, 20, 30, 40]:
            target_resi_list = [data[0] for data in similarity_table_dict[:x]]

            ref_length = len(target_resi_list)
            lg.info("similarity condition: ranktop {}:({} cases/{} cases)".format(
                x, ref_length, dict_length))
            lg.debug("target_resi_list: \n{}".format(target_resi_list))

            output_pdb_prefix = "{fix}_ranktop{x}".format(fix=pdb_prefix, x=x)

            _sub_main(target_resi_list=target_resi_list,
                      contact_profile_diff_file=contact_profile_diff_file,
                      trajectory=trajectory,
                      output_pdb_prefix=output_pdb_prefix,
                      output_dir=output_dir,
                      exclude_list=exclude_residue_list,
                      logger=lg)

        if(similarity_mode == SIMILARITY_TANIMOTO):
            # Residues with a TANIMOTO coefficient threshold smaller than the specified value
            for cutoff in [0.5, 0.7, 0.9]:
                target_resi_list = [
                    data[0] for data in similarity_table_dict if data[1] <= cutoff]

                ref_length = len(target_resi_list)
                lg.info("similarity condition: similarity <= {}:({} cases/{} cases)".format(
                    cutoff, ref_length, dict_length))
                lg.debug("target_resi_list: \n{}".format(target_resi_list))

                output_pdb_prefix = "{fix}_cutoff{x}".format(fix=pdb_prefix, x=cutoff)

                _sub_main(target_resi_list=target_resi_list,
                          contact_profile_diff_file=contact_profile_diff_file,
                          trajectory=trajectory,
                          output_pdb_prefix=output_pdb_prefix,
                          output_dir=output_dir,
                          exclude_list=exclude_residue_list,
                          logger=lg)

    except Exception as e:
        lg.error('{}\n{}'.format(e, traceback.format_exc()))
        sys.exit(1)

    finally:
        lg.info('END similarity_drawer')
        for handle in lg.handlers[::-1]:
            handle.close()
            lg.removeHandler(handle)


def _is_include_region(residue_id: str, region: dict) -> bool:
    """ある残基IDが残基領域に含まれているかを判定する関数

    Args:
        residue_id (str): 残基ID
        region (dict): Chainごとの残基領域がリスト化された辞書
                       ex) dict['A'] = ['100', '300_400C']
                       ex) dict['B'] = ['_200', '300_']

    Returns:
        bool: 残基領域に含まれる場合は True
    """
    chain = residue_id[0]
    resnum = residue_id[1:]

    if(chain not in region.keys()):
        return False

    for chain_region in region[chain]:
        if('_' in chain_region):
            # 残基領域が範囲として指定されている場合
            # 想定される残基領域文字列
            # 100_200A: 残基番号 100 から 残基番号 200 icode A
            # _100: 残基番号が 100以下
            # 100_: 残基番号が 100以上

            lower = chain_region.split('_')[0]
            upper = chain_region.split('_')[1]

            # <含まれているかの判定方法>
            # 残基領域の下端、上端と対象残基を自然順に並べ替え、その順番より判定
            if(len(lower) == 0):
                # 残基番号 upper 以下に resnum が含まれているかを判定
                sorted_list = sorted([resnum, upper], key=_natural_keys)

                if(resnum == sorted_list[0]):
                    return True

            elif(len(upper) == 0):
                # 残基番号 lower 以上に resnum が含まれているかを判定
                sorted_list = sorted([lower, resnum], key=_natural_keys)

                if(resnum == sorted_list[1]):
                    return True

            else:
                # 残基番号 lower 以上 upper 以下に resnum が含まれているかを判定
                sorted_list = sorted([lower, resnum, upper], key=_natural_keys)

                if(resnum == sorted_list[1]):
                    return True

        else:
            if(resnum == chain_region):
                return True
    return False


def get_partner_residue(contact_profile_diff_file: str, target_residue: list,
                        frequency_cutoff_upper: float = None,
                        frequency_cutoff_lower: float = None,
                        exclude_list: list = [], logger=None) -> dict:
    """対象残基に対して contact frezuency が閾値条件を満たす残基の残基IDを取得する関数

    残基ID: [chain][残基番号][icode]

    Args:
        contact_profile_diff_file (str): - contact profile mapのdifferenceタブのテーブルをCSV変換したファイル
                                         - 1行目と1列目に残基IDが記載
        target_residue (list): 対象残基の残基IDリスト
        frequency_cutoff_upper (float, optional): 閾値条件-上限値. Defaults to None.
        frequency_cutoff_lower (float, optional): 閾値条件-下限値. Defaults to None.
        exclude_list (list, optional): 探索除外リスト. Defaults to []
        logger (logger, optional): ロガー. Defaults to None.

    Returns:
        dict: 対象残基に対するパートナー残基のリスト
              key = 残基ID
              value = list([残基ID,...])
    """
    assert frequency_cutoff_upper is not None or frequency_cutoff_lower is not None,\
           "frequency_cutoff_upper/lower のいずれかに値が必要"

    if(logger is None):
        logger = get_logger(name=__name__, log_file_path="", is_stdout=False)

    logger.info("START get partner residue")
    logger.debug(
        "args: \n"
        "contact_profile_diff_file: {diff}\n"
        "target_residue: {target}\n"
        "contact_frequency_cutoff_upper: {upper}\n"
        "contact_frequency_cutoff_lower: {lower}\n"
        "exclude_list: {exclude}".format(
            diff=contact_profile_diff_file, target=target_residue,
            upper=frequency_cutoff_upper, lower=frequency_cutoff_lower,
            exclude=exclude_list
        )
    )

    if(all([frequency_cutoff_upper is not None,
            frequency_cutoff_lower is not None])):
        assert frequency_cutoff_lower < frequency_cutoff_upper
        logger.info("閾値条件: {} < contact_frequency < {}".format(
            frequency_cutoff_lower, frequency_cutoff_upper))

    elif(all([frequency_cutoff_upper is not None,
              frequency_cutoff_lower is None])):
        logger.info("閾値条件: {} <= contact_frequency".format(frequency_cutoff_upper))

    elif(all([frequency_cutoff_lower is not None,
              frequency_cutoff_upper is None])):
        logger.info("閾値条件: contact_frequency <= {}".format(frequency_cutoff_lower))

    target_residue_pairs = dict()
    with open(contact_profile_diff_file, 'r') as map_rf:
        file_header = True
        col_residue_idx = None
        for row_line in map_rf.readlines():

            # 列番号：残基IDをヘッダーから取得する。
            if(file_header):
                col_residue_idx = row_line.strip().split(',')[1:]
                file_header = False

            cols = row_line.strip().split(',')
            residue_id = cols[0].strip()

            if(residue_id not in target_residue):
                continue

            target_residue_pairs[residue_id] = list()

            for idx, val in enumerate(cols[1:]):

                if(len(val) == 0):
                    # contact profile diff map の値が空欄の場合はパートナーから除外する
                    # contact profile diff map が空白
                    # # → cp1, cp2 の contact profile の値が両方 0 であるという意味
                    continue

                # パートナー候補残基の残基ID
                candidate_id = col_residue_idx[idx].strip()

                if(candidate_id in residue_id):
                    # 対象残基自身をパートナー残基として取得しない
                    continue
                elif(candidate_id in exclude_list):
                    # パートナー候補残基が除外リストに含まれる場合パートナーから除外する
                    continue

                # contact profile diff 値
                val = float(val)

                if(all([frequency_cutoff_upper is not None,
                        frequency_cutoff_lower is not None])):

                    if(frequency_cutoff_lower < val < frequency_cutoff_upper):
                        target_residue_pairs[residue_id].append(candidate_id)

                        logger.debug("[type I] partner_residue_id: '{}'\t"
                                     "contact frequency: {}".format(candidate_id, val))

                elif(all([frequency_cutoff_upper is not None,
                          frequency_cutoff_lower is None]) and val >= frequency_cutoff_upper):
                    target_residue_pairs[residue_id].append(candidate_id)

                    logger.debug("[type II] partner_residue_id: '{}'\t"
                                 "contact frequency: {}".format(candidate_id, val))

                elif(all([frequency_cutoff_lower is not None,
                          frequency_cutoff_upper is None]) and val <= frequency_cutoff_lower):
                    target_residue_pairs[residue_id].append(candidate_id)

                    logger.debug("[type III] partner_residue_id: '{}'\t"
                                 "contact frequency: {}".format(candidate_id, val))

            logger.info("residue_pairs: target residue_id {} - total {} partners".format(
                residue_id, len(target_residue_pairs[residue_id])))
            logger.debug("residue_pairs: target_residue_pairs[{}] = {}".format(
                residue_id, target_residue_pairs[residue_id]))

    logger.info("END get partner residue")
    return target_residue_pairs


def save_pdb(input_pdb: str, output_pdb: str, similarity_table_dict: dict,
             target_as_psude_residue_list: list, target_pair_dict: dict,
             pseudo_name: str, logger=None):
    """相互作用情報（疑似ATOM行、CONECT行）を追加したPDBファイルを出力する関数

    Args:
        input_pdb (str): 入力PDB名
        output_pdb (str): 出力PDB名
        similarity_table_dict (dict): 残基IDとコンタクトプロファイルの類似度が対応した辞書
        target_as_psude_residue_list (list): 類似度が閾値条件を満たした残基の残基ID
        target_pair_dict (dict): 対象残基（target_as_psude_residue_list）に対して
                                 contact profile diff値が閾値条件を満たした残基の残基IDのリストが対応した辞書
        pseudo_name (str): 疑似ATOM行の原子名
        logger (logging.Logger, optional): ロガー. Defaults to None.
    """
    if(logger is None):
        logger = get_logger(name=__name__, log_file_path="", is_stdout=False)

    logger.info("START PDB ({}) への出力".format(output_pdb))

    logger.debug(
        "args: \n"
        "input_pdb: {}\n"
        "output_pdb: {}\n"
        "similarity_table_dict: {}\n"
        "target_as_psude_residue_list: {}\n"
        "target_pair_dict: {}\n"
        "pseudo_name: {}".format(
            input_pdb, output_pdb, similarity_table_dict, target_as_psude_residue_list,
            target_pair_dict, pseudo_name
        )
    )

    # 疑似ATOM行とCONECT行が追記されたPDBを作成
    target_atom_line_to_build_pseudo_atom = dict()
    last_atom_id_of_chains = dict()
    last_residue_id_of_chains = dict()

    # 疑似ATOM行作成対象となるパートナー遺伝子の一覧を取得
    partners_list = set()
    for partner_residues in target_pair_dict.values():
        if(len(partners_list) == 0):
            partners_list = set(partner_residues)
        else:
            partners_list = partners_list.union(set(partner_residues))

    logger.debug("partner_residue_list: \n{}".format(partners_list))

    with open(input_pdb, 'r') as rpdb,\
         open(output_pdb, 'w') as wpdb:

        # 各残基IDごとに中心座標を保持する変数
        coords = dict()

        for line in rpdb.readlines():
            if(line.startswith('ATOM')):
                name = line[12:16].strip()
                chain = line[21]
                resi = line[22:26].strip()
                icode = line[26].strip()

                x_coord = float(line[30:38].strip())
                y_coord = float(line[38:46].strip())
                z_coord = float(line[46:54].strip())

                element_symbol = line[76:78].strip()

                residue_id = '{}{}{}'.format(chain, resi, icode)

                if(residue_id in target_as_psude_residue_list
                   or residue_id in partners_list):
                    if(name == 'CA'):
                        # 疑似 ATOM行を作成する対象残基のCA行を取得
                        target_atom_line_to_build_pseudo_atom[residue_id] = line
                        logger.debug("PSUDE target line: {}".format(line.rstrip('\n')))

                    elif(element_symbol != 'H'):
                        # ACE, リガンド, イオンの場合は、CAが存在しない
                        # その場合、疑似ATOM行の原子座標を重原子の中心座標として出力するため
                        # 水素原子以外の原子座標を保持する

                        if(residue_id not in target_atom_line_to_build_pseudo_atom.keys()):
                            target_atom_line_to_build_pseudo_atom[residue_id] = line
                            logger.debug("PSUDE target line: {}".format(line.rstrip('\n')))

                        if(residue_id not in coords.keys()):
                            coords[residue_id] = dict()
                            coords[residue_id]['x_coord'] = list()
                            coords[residue_id]['y_coord'] = list()
                            coords[residue_id]['z_coord'] = list()

                        coords[residue_id]['x_coord'].append(x_coord)
                        coords[residue_id]['y_coord'].append(y_coord)
                        coords[residue_id]['z_coord'].append(z_coord)

                # 疑似 ATOM 行作成の為に末尾にあたる原子番号と残基番号をChainごとに取得
                # ATOM行の読み込み後に値が原子番号と残基番号の末尾となる
                atom_id = int(line[6:11].strip())

                last_atom_id_of_chains[chain] = atom_id
                if(len(icode) == 0):
                    last_residue_id_of_chains[chain] = int(resi)
                else:
                    last_residue_id_of_chains[chain] = resi + icode

                # b-factor の値を"0"に修正し出力
                wpdb.write("{front}{bfactor:>6.2f}{end}".format(
                    front=line[:60],
                    bfactor=0,
                    end=line[66:]
                ))

            elif(line.startswith('END')):
                # END行の前に疑似ATOM行とCONECT行を追加

                # 残基番号と疑似ATOM行の原子番号との対応表となる辞書
                residue2psude_atom_id_dict = dict()
                residue2partner_psude_atom_id_dict = dict()

                for target_atom_line in target_atom_line_to_build_pseudo_atom.values():
                    name = target_atom_line[12:16].strip()
                    chain = target_atom_line[21]
                    resi = target_atom_line[22:26].strip()
                    icode = target_atom_line[26].strip()

                    x_coord = float(target_atom_line[30:38].strip())
                    y_coord = float(target_atom_line[38:46].strip())
                    z_coord = float(target_atom_line[46:54].strip())

                    residue_id = '{}{}{}'.format(chain, resi, icode)

                    if(name != 'CA'):
                        # ACE, リガンド, イオンの場合は、CAが存在しない
                        # その場合、疑似ATOM行の原子座標を重原子の中心座標として出力
                        x_coord, y_coord, z_coord = _calc_center_coord(
                            x_coords=coords[residue_id]['x_coord'],
                            y_coords=coords[residue_id]['y_coord'],
                            z_coords=coords[residue_id]['z_coord'])

                    is_similarity_res = residue_id in target_as_psude_residue_list
                    is_partner_res = residue_id in partners_list

                    for key, val in {"similarity": is_similarity_res,
                                     "partner": is_partner_res}.items():
                        if(val is False):
                            continue

                        # 疑似ATOM行の原子番号とATOM番号を割り振る
                        psude_atom_id = last_atom_id_of_chains[chain] + 1
                        psude_resi = None
                        if(type(last_residue_id_of_chains[chain]) is int):
                            psude_resi = int(last_residue_id_of_chains[chain]) + 1
                        else:
                            # 挿入コードを含んだ残基番号の場合
                            # 挿入コードを除いた残基番号に+1した値を疑似ATOM行の残基番号とする
                            psude_resi = int(last_residue_id_of_chains[chain][:-1]) + 1

                        psude_res_name = "PSD" if key == "similarity" else "XXX"
                        pseudo_element = ""
                        if("CA" in pseudo_name):
                            pseudo_element = "C"
                        elif("O" in pseudo_name):
                            pseudo_element = "O"
                        elif("N" in pseudo_name):
                            pseudo_element = "N"

                        # 原子番号と残基番号の末尾となる値を更新
                        last_atom_id_of_chains[chain] = psude_atom_id
                        last_residue_id_of_chains[chain] = psude_resi

                        similarity_value = 0
                        if(residue_id in target_as_psude_residue_list):
                            # 類似度の値が閾値条件を満たした残基は
                            # b-factor の値を類似度で出力
                            similarity_value = similarity_table_dict[residue_id]

                        # 疑似ATOM行の出力
                        wpdb.write(_build_psd_atom_line(
                            atom_line=target_atom_line,
                            pseudo_atom_id=psude_atom_id,
                            pseudo_residue_id=psude_resi,
                            similarity=similarity_value,
                            pseudo_name=pseudo_name,
                            pseudo_res_name=psude_res_name,
                            pseudo_x_coord=x_coord,
                            pseudo_y_coord=y_coord,
                            pseudo_z_coord=z_coord,
                            pseudo_element=pseudo_element))

                        # 残基番号と疑似ATOM行の原子番号を対応付ける
                        if(key == "similarity"):
                            residue2psude_atom_id_dict[residue_id] = psude_atom_id
                        elif(key == "partner"):
                            residue2partner_psude_atom_id_dict[residue_id] = psude_atom_id

                logger.info("similarity残基と疑似ATOM行の原子番号の対応表:\n{}".format(residue2psude_atom_id_dict))
                logger.info("partner残基と疑似ATOM行の原子番号の対応表:\n{}".format(
                    residue2partner_psude_atom_id_dict))

                # 疑似原子の CONECT 行の出力
                for key in target_pair_dict.keys():
                    # CONECT 行の書式に合わせて、
                    # パートナー残基を４つごとに記載する。
                    for idx in range(0, len(target_pair_dict[key]), 4):
                        partners_resi = target_pair_dict[key][idx:idx + 4]
                        partner_atom_id = [residue2partner_psude_atom_id_dict[residue_id]
                                           for residue_id in partners_resi]

                        if(len(partner_atom_id) < 4):
                            partner_atom_id.extend([""] * (4 - len(partner_atom_id)))

                        wpdb.write(_build_conect_line(
                            atom_id=residue2psude_atom_id_dict[key],
                            bonded_id1=partner_atom_id[0],
                            bonded_id2=partner_atom_id[1],
                            bonded_id3=partner_atom_id[2],
                            bonded_id4=partner_atom_id[3]))

                        # パートナー残基のCONECT行を記載
                        for partner_id in partner_atom_id:
                            if(partner_id == ""):
                                break
                            wpdb.write(_build_conect_line(
                                atom_id=partner_id,
                                bonded_id1=residue2psude_atom_id_dict[key]
                            ))

                wpdb.write(line)
                break
            else:
                wpdb.write(line)
        logger.info("END PDB ({}) への出力".format(output_pdb))


def _calc_center_coord(x_coords: list, y_coords: list, z_coords: list) -> tuple:
    """3次元の中心座標を算出する関数

    Args:
        x_coords (list): x座標リスト
        y_coords (list): y座標リスト
        z_coords (list): z座標リスト

    Returns:
        tuple: 中心座標（x, y, z）
    """
    num_coords = len(x_coords)
    assert num_coords > 0
    assert len(y_coords) == num_coords
    assert len(z_coords) == num_coords

    coord = np.array([0.0, 0.0, 0.0])
    for x, y, z in zip(x_coords, y_coords, z_coords):
        coord += np.array([x, y, z])
    center_coord = coord / num_coords
    return center_coord[0], center_coord[1], center_coord[2]


def _build_psd_atom_line(atom_line: str, pseudo_atom_id: str, pseudo_residue_id: str,
                         similarity: float, pseudo_name: str = " CA",
                         pseudo_res_name: str = "PSD", pseudo_x_coord: float = 0.0,
                         pseudo_y_coord: float = 0.0, pseudo_z_coord: float = 0.0,
                         pseudo_element: str = "C") -> str:
    """疑似ATOM行を返す関数

    Args:
        atom_line (str): 疑似ATOM行の基になるATOM行
        pseudo_atom_id (str): 疑似ATOM行の原子番号
        pseudo_residue_id (str): 疑似ATOM行の残基番号
        similarity (float): 疑似ATOM行の原子名 tempFactor に代替する値
        pseudo_name (str, optional): 疑似ATOM行の原子名. Defaults to " X".
        pseudo_res_name (str, optional): 疑似ATOM行の残基名. Defaults to "PSD".

    Returns:
        str: 書式に則った疑似ATOM行文字列
    """
    psd_atom_line = "ATOM  {serial:>5} {name:<4}{altLoc:1}{res_name:>3} {chain:1}"\
                    "{res_seq:>4}{icode:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6}"\
                    "{tempfactor:>6.2f}{element:>12}{charge:>2}\n".format(
                        serial=pseudo_atom_id,
                        name=pseudo_name,
                        altLoc=atom_line[16],
                        res_name=pseudo_res_name,
                        chain=atom_line[21],
                        res_seq=pseudo_residue_id,
                        icode=atom_line[26],
                        x=pseudo_x_coord,
                        y=pseudo_y_coord,
                        z=pseudo_z_coord,
                        occupancy=atom_line[54:60],
                        tempfactor=similarity,
                        element=pseudo_element,
                        charge=atom_line[78:80]
                    )
    return psd_atom_line


def _build_conect_line(atom_id: str, bonded_id1: str, bonded_id2: str = "",
                       bonded_id3: str = "", bonded_id4: str = "") -> str:
    """指定された原子番号からCONECT行を作成し返す関数

    Args:
        atom_id (str): CONECT行 7-11 行に対応
        bonded_id1 (str): CONECT行 12-16 に対応
        bonded_id2 (str, optional): CONECT行 17-21 に対応. Defaults to "".
        bonded_id3 (str, optional): CONECT行 22-26 に対応. Defaults to "".
        bonded_id4 (str, optional): CONECT行 27-31 に対応. Defaults to "".

    Returns:
        str: 書式に則ったCONECT行文字列
    """
    conect_line = "CONECT{atom_id:>5}{bonded1:>5}{bonded2:>5}"\
                  "{bonded3:>5}{bonded4:>5}\n".format(
                      atom_id=atom_id,
                      bonded1=bonded_id1,
                      bonded2=bonded_id2,
                      bonded3=bonded_id3,
                      bonded4=bonded_id4)
    return conect_line


def _validate(region: str):
    """main 引数のバリデーションチェックをする関数

    Args:
        region (str): 計算対象除外残基領域
    """
    # リージョンの下端、上端
    RESI_MIN = -99
    RESI_MAX = 999

    if(type(region) is not str):
        raise TypeError('[TypeError] region "{}" is not str!!'.format(region))
    elif(region.endswith('\n')):
        raise ValueError('[ValueError] region "{}" contains line breaks!!'.format(region))

    region_dict = _parse_region(region)

    for chain in region_dict.keys():
        if(len(chain) == 0):
            continue

        if(re.match(r'^([A-Z]|[a-z])$', chain) is None):
            raise ValueError('[ValueError] "{}" is not chain!! in "{}"'.format(chain, region))

        for resi_region in region_dict[chain]:
            for resi in resi_region.split('_'):
                if(len(resi) == 0):
                    continue

                icode = None
                if(not resi[-1].isdecimal()):
                    icode = resi[-1]
                    resi = resi[:-1]

                if(resi is not None and re.match(r'(-|)[0-9]+', resi) is None):
                    raise ValueError('[ValueError] "{}" is not Interger!! in "{}"'.format(
                                                                    resi, region))

                elif(resi is not None and int(resi) < RESI_MIN):
                    raise ValueError('[ValueError] "{}" less than {}!! in "{}"'.format(
                                                                    resi, RESI_MIN, region))

                elif(resi is not None and int(resi) > RESI_MAX):
                    raise ValueError('[ValueError] "{}" more than {}!! in "{}"'.format(
                                                                    resi, RESI_MAX, region))

                if(icode is not None and re.match(r'^([A-Z]|[a-z])$', icode) is None):
                    raise ValueError('[ValueError] "{}" is not icode!! in "{}"'.format(
                                                                    icode, region))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        'contact_profile',
        help='contact profile diff file.'
    )
    parser.add_argument(
        'similarity_table',
        help='Similarity table file.'
    )
    parser.add_argument(
        'trajectory',
        help='multi-PDB / PDB'
    )
    parser.add_argument(
        '--noregion', '-n',
        help='Exclude residues.\n'
             'Format: (<chain>):(<resnum>) \n'
             'Multiple specifications can be specified by separating them with ","',
        default=''
    )
    parser.add_argument(
        '-s', '--similarity',
        help='Similarity coefficient to be visualized',
        type=str,
        choices=[SIMILARITY_TANIMOTO, SIMILARITY_EUCLIDEAN],
        default=SIMILARITY_TANIMOTO
    )

    args = parser.parse_args()
    contact_profile_diff = args.contact_profile
    similarity_table = args.similarity_table
    trajectory = args.trajectory
    exclude_region = args.noregion
    similarity_mode = args.similarity
    main(contact_profile_diff_file=contact_profile_diff,
         similarity_table_file=similarity_table,
         trajectory=trajectory,
         exclude_region=exclude_region,
         similarity_mode=similarity_mode,
         output_dir=os.getcwd())
