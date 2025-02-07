# -*- coding: utf-8 -*-
import numpy as np
import os
import re
import sys
from os.path import exists

from biopandas.pdb import PandasPdb as ppdb
from concurrent import futures
from log import get_logger


# 重原子のみでコンタクトプロファイルを作成する際の閾値
# CONAN, Difference contact map による解析例での設定値
ONLY_HEAVY_ATOM_CONTACT_CUTOFF = 5.0

# 水素原子 + 重原子でコンタクトプロファイルを作成する際の閾値
# Rodriguez-Bussey, 2018等 による解析例での設定値
ALL_ATOM_CONTACT_CUTOFF = 4.5

CA_ATOM_CONTACT_CUTOFF = 8.0
CB_ATOM_CONTACT_CUTOFF = 10.0

def profile_main(trajectory: str,
                 region: str = "",
                 contact_cutoff: float = None,
                 num_process: int = 1,
                 out_dist_map: bool = False,
                 target_atom: int = 0,
                 contact_file: str = '',
                 distance_file: str = '',
                 log_file_path: str = ''):
    """Trajectory を読み込み、平均距離行列と平均維持率行列を出力する

    Args:
        trajectory (str): マルチPDB / PDB
        region (str): 計算対象残基領域
        contact_cutoff (float, optional): コンタクト判定閾値. Defaults to None.
                                          target_atom=0. Defaults to ALL_ATOM_CONTACT_CUTOFF.
                                          target_atom=1. Defaults to ONLY_HEAVY_ATOM_CONTACT_CUTOFF.
                                          target_atom=2. Defaults to CA_ATOM_CONTACT_CUTOFF.
                                          target_atom=3. Defaults to CB_ATOM_CONTACT_CUTOFF.
        num_process (int, optional): マルチプロセス実行時のプロセス数. Defaults to 1.
        out_dist_map (bool, optional): 中間ファイルの出力 (距離行列). Defaults to False.
        target_atom (int, optional): データ取得対象原子タイプ. Defaults to 0.
                                    0: 全原子
                                    1: 重原子
                                    2: CA原子
                                    3: CB原子
        contact_file (str, optional): 出力ディレクトリ. Defaults to ''.
        distance_file (str, optional): 出力ディレクトリ. Defaults to ''.
        log_file_path (str, optional): 出力ディレクトリ. Defaults to ''.
    """
    if all([exists(f) for f in [contact_file, distance_file, log_file_path]]):
        return None

    if(contact_cutoff is None):
        # Setting default values for contact determination thresholds
        if(target_atom == 0):
            contact_cutoff = ALL_ATOM_CONTACT_CUTOFF

        elif(target_atom == 1):
            contact_cutoff = ONLY_HEAVY_ATOM_CONTACT_CUTOFF

        elif(target_atom == 2):
            contact_cutoff = CA_ATOM_CONTACT_CUTOFF

        elif(target_atom == 3):
            contact_cutoff = CB_ATOM_CONTACT_CUTOFF

    lg = get_logger(name=__name__, log_file_path=log_file_path, is_stdout=True)
    lg.info("START Contact_profile_calculator")
    lg.info(
        "args: \n"
        "trajectory: {t}\n"
        "region: {r}\n"
        "contact_cutoff: {c}\n"
        "num_process: {n}\n"
        "out_dist_map: {m}\n"
        "target_atom: {a}".format(
            t=trajectory, r=region, c=contact_cutoff, n=num_process,
            m=out_dist_map, a=['all', 'heavy atoms', 'CA', 'CB'][target_atom])
    )

    try:
        _validate(trajectory, region, contact_cutoff, num_process)

        total_dist_map = None
        total_binary_map = None
        map_col_name = None
        total_snapshot_num = 0
        snapshot_num_list = list()
        process_list = list()

        with futures.ProcessPoolExecutor(max_workers=num_process) as executor:
            for entry, df_atom in _generate_data_frame(trajectory=trajectory, region=region,
                                                       target_atom=target_atom,
                                                       logger=lg):
                total_snapshot_num += 1
                snapshot_num_list.append(entry)
                if(map_col_name is None):
                    map_col_name = sorted(list(df_atom['chain_resi_icode'].unique()),
                                          key=_natural_keys)
                process_list.append(
                    executor.submit(_calc_dist_map, df_atom=df_atom, snapshot=entry)
                )
            futures.wait(fs=process_list, return_when=futures.FIRST_EXCEPTION)

        for process_idx in range(len(process_list)):
            process_execption = process_list[process_idx].exception()
            if(process_execption is not None):
                raise process_execption

            dist_map = process_list[process_idx].result()
            entry = snapshot_num_list[process_idx]
            # if(out_dist_map):
            #     output = os.path.join(output_dir, "distance_map_{}.csv".format(entry))
            #     _save_matrix(dist_map, map_col_name, output)

            if(total_dist_map is None):
                total_dist_map = dist_map
            else:
                total_dist_map += dist_map

            # 距離行列を閾値で二値化した行列の総和の算出
            if(total_binary_map is None):
                total_binary_map = (dist_map < contact_cutoff) * 1
            else:
                total_binary_map += (dist_map < contact_cutoff) * 1

        # 二値化行列の対角成分をゼロにする。
        np.fill_diagonal(total_binary_map, 0)

        # 平均距離行列、維持率行列の算出
        assert total_snapshot_num > 0, "No Trajectory!!"
        mean_dist_map = total_dist_map / total_snapshot_num
        mean_binary_map = total_binary_map / total_snapshot_num

        _save_matrix(mean_dist_map, map_col_name, distance_file)
        _save_matrix(mean_binary_map, map_col_name, contact_file)
    except Exception as e:
        lg.error(e)
        sys.exit(1)
    finally:
        lg.info("END Contact_profile_calculator")
        for handle in lg.handlers[::-1]:
            handle.close()
            lg.removeHandler(handle)


def _calc_dist_map(df_atom, logger=None, snapshot=None):
    """trajectory のデータフレームから距離行列を求める関数

    Args:
        df_atom (pandas.core.frame.DataFrame): ATOM行の情報が含まれたデータフレーム
        logger (logging.Logger, optional): ロガー. Defaults to None.

    Returns:
        numpy.ndarray: 距離行列
    """
    if(logger is None):
        logger = get_logger(name=__name__, log_file_path="", is_stdout=False)

    if(snapshot is None):
        logger.info("Start distance matrix calculation")
    else:
        logger.info("Start distance matrix calculation: MODEL={}".format(snapshot))

    resi_list = sorted(list(df_atom['chain_resi_icode'].unique()), key=_natural_keys)
    logger.debug("Target Residue list: {}".format(resi_list))

    dist_map = np.zeros(len(resi_list) ** 2).reshape(len(resi_list), len(resi_list))

    for row_i in range(len(resi_list)):
        row_resi = resi_list[row_i]
        row_resi_coords = list()
        df = df_atom[df_atom['chain_resi_icode'] == row_resi]
        for atom_num in range(len(df)):
            x = df.iloc[atom_num]['x_coord']
            y = df.iloc[atom_num]['y_coord']
            z = df.iloc[atom_num]['z_coord']
            row_resi_coords.append(np.array([x, y, z]))

        for col_i in range(row_i + 1, len(resi_list)):
            col_resi = resi_list[col_i]
            col_resi_coords = list()
            df = df_atom[df_atom['chain_resi_icode'] == col_resi]

            for atom_num in range(len(df)):
                x = df.iloc[atom_num]['x_coord']
                y = df.iloc[atom_num]['y_coord']
                z = df.iloc[atom_num]['z_coord']
                col_resi_coords.append(np.array([x, y, z]))

            shortest_distance = None
            for row_coord in row_resi_coords:
                for col_coord in col_resi_coords:
                    dist = np.linalg.norm(row_coord - col_coord)

                    if(shortest_distance is None or dist < shortest_distance):
                        shortest_distance = dist

            dist_map[row_i][col_i] = shortest_distance

    assert np.array_equal(dist_map, np.triu(dist_map))
    dist_map = dist_map + dist_map.T - np.diag(dist_map.diagonal())

    if(snapshot is None):
        logger.info("End distance matrix calculation")
    else:
        logger.info("End distance matrix calculation: MODEL={}".format(snapshot))
    return dist_map


def _save_matrix(matrix: list, col_name_list: list, output: str):
    """行列をファイル出力する関数

    Args:
        matrix (list): 出力する行列
        col_name_list (list): 行名/列名のヘッダー
        output (str): 出力ファイル
    """
    assert len(matrix) == len(col_name_list),\
           "len(matrix)={} != len(col_name_list)={}".format(len(matrix), len(col_name_list))
    with open(output, "w") as w:
        w.write(",".join([""] + col_name_list) + "\n")
        for row_i in range(len(matrix)):
            str_dist_row = [str(d) for d in matrix[row_i]]
            w.write(",".join([col_name_list[row_i]] + str_dist_row) + "\n")


def _natural_keys(text: str):
    """sort 関数のkey引数への組み込み関数
    自然順に並び返す
    Args:
        text (str): 文字列

    """
    return [_atoi(c) for c in re.split(r'(-\d+|\d+)', text)]


def _atoi(text: str):
    """数値文字を数値として返す関数

    Args:
        text (str): 文字列

    Returns:
        str/int: 文字列が数値ならばint, それ以外は str
    """
    return int(text) if _is_num(text) else text


def _is_num(x: str) -> bool:
    """数値文字(整数)を判定する関数

    Args:
        x (str): 文字列

    Returns:
        bool: 文字列が数値ならば True
    """
    if(re.match(r'^(-|)[0-9]+$', x) is not None):
        return True
    else:
        return False


def _generate_data_frame(trajectory: str, region: str, target_atom: int = 0, logger=None):
    """PDBを読み込んで対象残基のデータフレームを返す関数

    Args:
        trajectory (str): multi PDB / PDB
        region (str): 計算対象残基領域
        target_atom (int, optional): データ取得対象原子タイプ. Defaults to 0.
                                    0: 全原子
                                    1: 重原子
                                    2: CA原子
                                    3: CB原子
        logger (logging.Logger, optional): ロガー. Defaults to None.

    Yields:
        Pandas: 対象残基のデータフレーム
    """
    if(logger is None):
        logger = get_logger(name=__name__, log_file_path="", is_stdout=False)

    logger.info("START READ Trajectory Data.")
    logger.debug(
        "args: \n"
        "trajectory: {t}\n"
        "region: {r}\n"
        "target atom: {a}".format(
            t=trajectory,
            r=region,
            a=['all', 'heavy atoms', 'CA', 'CB'][target_atom])
    )

    target_region = _parse_region(region)

    df_pdb = ppdb().read_pdb(trajectory)
    # ATOM, HETATM, ANISOU 以外の行からスナップショット認識行（MODEL, ENDMDL）を取得
    df_others = df_pdb.df['OTHERS'].query('record_name in ["MODEL", "ENDMDL"]')
    model_line_idx = None
    endmdl_line_idx = None
    entry = None

    # データフレームに対するクエリー文字列をスナップショット番号ごとに保持する
    atom_query_dict = dict()

    if(len(df_others) == 0):
        # マルチでないPDBを読み込んだ場合
        atom_query_dict[''] = ""

    for idx in range(len(df_others)):
        record_name = df_others.iloc[idx]['record_name']
        if(record_name == 'MODEL'):
            assert model_line_idx is None,\
                  'There is no "ENDMDL" before "MODEL": {}'.format(model_line_idx)

            model_line_idx = df_others.iloc[idx]['line_idx']
            entry = df_others.iloc[idx]['entry'].strip()

        elif(record_name == 'ENDMDL'):
            assert endmdl_line_idx is None,\
                   'There is no "MODEL" before "ENDMDL": {}'.format(endmdl_line_idx)
            endmdl_line_idx = df_others.iloc[idx]['line_idx']

        if(model_line_idx is not None and endmdl_line_idx is not None):
            atom_query_dict[entry] = '({s} < line_idx < {e})'.format(
                s=model_line_idx, e=endmdl_line_idx)

            model_line_idx = None
            endmdl_line_idx = None

    for entry, atom_query in atom_query_dict.items():
        # ATOM行情報から水素を除外するためのクエリー文字列を追加
        if(target_atom == 1):
            if(len(atom_query) == 0):
                atom_query = 'element_symbol != "H"'
            else:
                atom_query += 'and (element_symbol != "H")'

        elif(target_atom == 2):
            if(len(atom_query) == 0):
                atom_query = 'atom_name == "CA"'
            else:
                atom_query += 'and (atom_name == "CA")'

        elif(target_atom == 3):
            if(len(atom_query) == 0):
                atom_query = 'atom_name == "CB"'
            else:
                atom_query += 'and (atom_name == "CB")'

        elif(target_atom != 0):
            raise "Undifined atom type number: {}".format(target_atom)

        df_atom = df_pdb.df['ATOM'].query(atom_query)

        region_query = ""
        for chain, resi_list in target_region.items():
            if(len(region_query) != 0):
                region_query += " or "

            if(len(resi_list) == 0):
                region_query += "(chain_id == '{c}')".format(c=chain)
                continue
            else:
                region_query += "(chain_id == '{c}' and ".format(c=chain)

            for i in range(len(resi_list)):
                if(i == 0):
                    region_query += "("
                else:
                    region_query += " or "

                resi = resi_list[i]
                if("_" in resi):
                    start = resi.split("_")[0]
                    start_icode = ""
                    if(len(start) != 0 and not _is_num(start[-1])):
                        start_icode = start[-1]
                        start = start[:-1]

                    end = resi.split("_")[1]
                    end_icode = ""
                    if(len(end) != 0 and not _is_num(end[-1])):
                        end_icode = end[-1]
                        end = end[:-1]

                    if(len(start) == 0):
                        df_end_resi = None
                        if(len(end_icode) == 1):
                            df_end_resi = df_atom.query(
                                '(residue_number == {e}) and (insertion == "{i}")'.format(
                                    e=end, i=end_icode))
                        else:
                            df_end_resi = df_atom.query('residue_number == {e}'.format(e=end))

                        assert len(df_end_resi['atom_number']) != 0,\
                               "Not Found resi={}!!".format(resi)

                        end_resi_atom_index = max(df_end_resi['atom_number'])
                        region_query += "atom_number <= {e}".format(e=end_resi_atom_index)

                    elif(len(end) == 0):
                        df_start_resi = None
                        if(len(start_icode) == 1):
                            df_start_resi = df_atom.query(
                                '(residue_number == {s}) and (insertion == "{i}")'.format(
                                    s=start, i=start_icode))
                        else:
                            df_start_resi = df_atom.query('residue_number == {s}'.format(s=start))

                        assert len(df_start_resi['atom_number']) != 0,\
                               "Not Found resi={}!!".format(resi)

                        start_resi_atom_index = min(df_start_resi['atom_number'])
                        region_query += "{s} <= atom_number".format(s=start_resi_atom_index)

                    else:
                        assert int(start) <= int(end), "invalid region {} <= {} in chain {}".format(
                            start, end, chain)

                        df_end_resi = None
                        if(len(end_icode) == 1):
                            df_end_resi = df_atom.query(
                                '(residue_number == {e}) and (insertion == "{i}")'.format(
                                    e=end, i=end_icode))
                        else:
                            df_end_resi = df_atom.query('residue_number == {e}'.format(e=end))

                        assert len(df_end_resi['atom_number']) != 0,\
                               "Not Found resi={}!!".format(resi)

                        end_resi_atom_index = max(df_end_resi['atom_number'])

                        df_start_resi = None
                        if(len(start_icode) == 1):
                            df_start_resi = df_atom.query(
                                '(residue_number == {s}) and (insertion == "{i}")'.format(
                                    s=start, i=start_icode))
                        else:
                            df_start_resi = df_atom.query('residue_number == {s}'.format(s=start))

                        assert len(df_start_resi['atom_number']) != 0,\
                               "Not Found resi={}!!".format(resi)

                        start_resi_atom_index = min(df_start_resi['atom_number'])
                        region_query += "{s} <= atom_number <= {e}".format(
                            s=start_resi_atom_index, e=end_resi_atom_index)

                else:
                    region_query += "residue_number == {}".format(resi)

                if(i == len(resi_list) - 1):
                    region_query += ")"

            region_query += ")"

        if(len(region_query) != 0):
            df_atom = df_atom.query(region_query)

        df_chain_resi_icode = df_atom['chain_id'] + df_atom['residue_number'].astype(str)\
                                                  + df_atom['insertion']
        df_atom = df_atom.assign(chain_resi_icode=df_chain_resi_icode)

        if(len(df_atom) == 0):
            logger.warning(
                "No target data !!\n"
                "MODEL {}\n"
                "DataFrame Query: {}".format(entry, atom_query)
            )
        else:
            logger.debug("DataFrame Query: {}".format(atom_query))
            logger.info("MODEL {}: DataFrame size (row, col) = {}".format(entry, df_atom.shape))

        yield entry,\
            df_atom[['atom_number', 'chain_id', 'chain_resi_icode',
                     'x_coord', 'y_coord', 'z_coord']]


def _validate(trajectory: str, region: str, contact_cutoff: float, num_process: int):
    """main 引数のバリデーションチェックをする関数

    Args:
        trajectory (str): multi-PDB / PDB
        region (str): 計算対象残基領域
        contact_cutoff (float): コンタクト判定閾値
        num_process (int): プロセス数
    """

    # リージョンの下端、上端
    RESI_MIN = -99
    RESI_MAX = 999

    if(type(trajectory) is not str):
        raise TypeError('[TypeError] trajectory "{}" is not str!!'.format(trajectory))
    elif(not trajectory.endswith('.pdb')):
        raise ValueError('[ValueError] trajectory "{}" is not supported!!'
                         '(".pdb" is supported.)'.format(trajectory))

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

    if(type(contact_cutoff) is not float and type(contact_cutoff) is not int):
        raise TypeError('[TypeError] contact_cutoff "{}" is not Float(or Integer)!!'.format(
                                                                                contact_cutoff))

    elif(contact_cutoff < 0):
        raise ValueError('[ValueError] contact_cutoff "{}" less than 0!!'.format(contact_cutoff))

    if(type(num_process) is not int):
        raise TypeError('[TypeError] num_process "{}" is not Integer!!'.format(num_process))

    elif(num_process < 0):
        raise ValueError('[ValueError] num_process "{}" less than 0!!'.format(num_process))


def _parse_region(region: str) -> dict:
    """残基指定書式を構文解析してChainごとに領域をリスト化する関数

    Args:
        region (str): 残基指定書式 (<chain>):(<resnum>)
                      ex) A:300_400C,B:-1_-20

    Returns:
        dict: Chainごとの領域をリスト化した辞書
              領域は ':' で区切られる
              ex) dict['A'] = ['300_400C']
                  dict['B'] = ['-1_-20']
    """
    region_dict = dict()
    if(len(region) == 0):
        return region_dict

    for single_region in region.split(','):
        if(len(single_region) == 0):
            continue
        elif(len(single_region.split(':')) != 2):
            raise ValueError('[ValueError] "{}" is Invalid format!! in "{}"'.format(
                                                                single_region, region))

        chain = single_region.split(':')[0]

        if(chain not in region_dict.keys()):
            region_dict[chain] = list()

        resnum = single_region.split(':')[1]
        if(len(resnum) == 0):
            continue
        else:
            region_dict[chain].append(resnum)

    return region_dict
#
#
# if(__name__ == '__main__'):
#
#     parser = argparse.ArgumentParser(
#         description='Read Trajectory (multi-PDB / PDB) '
#                     'and calculate the contact profile for each residue.',
#         formatter_class=RawTextHelpFormatter
#     )
#     parser.add_argument(
#         'trajectory',
#         help='multi-PDB / PDB',
#         type=str
#     )
#     parser.add_argument(
#         '-r', '--region',
#         help='Contact profile calculation target area\n'
#              'Format: (<chain>):(<resnum>) \n'
#              'Multiple specifications can be specified by separating them with ","',
#         default=''
#     )
#     parser.add_argument(
#         '-c', '--contact_cutoff',
#         help='Contact judgment threshold.\n'
#              'Option "-t 0". Defaults to  {0}.\n'
#              'Option "-t 1". Defaults to  {1}.\n'
#              'Option "-t 2". Defaults to  {2}.\n'
#              'Option "-t 3". Defaults to  {3}.\n'.format(
#                  ALL_ATOM_CONTACT_CUTOFF,
#                  ONLY_HEAVY_ATOM_CONTACT_CUTOFF,
#                  CA_ATOM_CONTACT_CUTOFF,
#                  CB_ATOM_CONTACT_CUTOFF
#              ),
#         type=float
#     )
#     parser.add_argument(
#         '-p', '--num_process',
#         help='Number of processes.(Run on multiprocessing)',
#         default=1,
#         type=int
#     )
#     parser.add_argument(
#         '-m', '--out_dist_map',
#         help='Output an intermediate file (distance matrix).',
#         action='store_true'
#     )
#     parser.add_argument(
#         '-t', '--target_atom',
#         help='Interaction analysis target atoms.\n'
#              '0: All atoms\n'
#              '1: Heavy atoms (Atoms excluding hydrogen atoms)[Default]\n'
#              '2: CA atoms\n'
#              '3: CB atoms\n',
#         default=1,
#         choices=[0, 1, 2, 3],
#         type=int
#     )
#
#     args = parser.parse_args()
#
#     profile_main(trajectory=args.trajectory,
#          region=args.region,
#          contact_cutoff=args.contact_cutoff,
#          num_process=args.num_process,
#          out_dist_map=args.out_dist_map,
#          target_atom=args.target_atom)
