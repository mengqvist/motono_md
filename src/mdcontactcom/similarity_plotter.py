import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import re
import sys
import traceback

from argparse import RawTextHelpFormatter
from contact_profile_calculator import _parse_region, _natural_keys
from similarity_drawer import _is_include_region
from log import get_logger


def main(similarity_table_file: str, bin: int = 50, exclude_region: str = '',
         output_dir: str = None):
    """similarity table を読み込んでグラフを描画する関数

    Args:
        similarity_table_file (str): similarity table ファイル
                                     1 列目: 残基ID
                                     2 列目: tanimoto係数によるSimilarityの値
                                     3 列目: ユークリッド距離によるSimilarityの値
        output_dir (str, optional): 出力先. Defaults to None.
        bin (int): 目盛りの刻み幅. Defaults to 50.
    """
    if(output_dir is None):
        output_dir = os.path.join(os.getcwd(), 'similarity_plotter')
    os.makedirs(output_dir, exist_ok=False)

    log_file_path = os.path.join(output_dir, "similarity_plotter.log")
    lg = get_logger(name=__name__, log_file_path=log_file_path, is_stdout=True)
    lg.info('START similarity_plotter')

    lg.info(
        "args: \n"
        "similarity_table_file: {t}\n"
        "bin: {b}\n"
        "exclude_region: {e}\n"
        "output_dir: {o}".format(
            t=similarity_table_file, b=bin, e=exclude_region, o=output_dir)
    )

    try:
        _validate(bin=bin, region=exclude_region)

        residue_id_data = list()
        tanimoto_data = list()
        euclidean_data = list()

        # 処理対象から除外する残基領域
        exclude_region = _parse_region(exclude_region)

        with open(similarity_table_file, 'r') as rf:
            header = True
            # プロットするデータの取得
            for line in rf.readlines():
                items = line.strip().split(',')
                assert len(items) >= 3, 'Missing columns!!: {}'.format(line.strip())

                if(header):
                    # 各列に想定する値
                    # 1列目：残基ID
                    # 2列目：tanimoto係数によるSimilarityの値
                    # 3列目：ユークリッド距離によるSimilarityの値

                    if(items[0] != 'residue_id'):
                        lg.warning('The column name of the first column is not "residue_id"!'
                                   ': {}'.format(items[0]))

                    if(items[1] != 'Tanimoto coefficient'):
                        lg.warning('The column name of the second column is not '
                                   '"Tanimoto coefficient"!: {}'.format(items[1]))

                    if(items[2] != 'Euclidean distance'):
                        lg.warning('The column name of the third column is not '
                                   '"Euclidean distance"!: {}'.format(items[2]))

                    header = False
                    continue

                if(_is_include_region(items[0], exclude_region)):
                    # 指定残基領域に含まれる残基の除外
                    continue

                residue_id_data.append(items[0])
                tanimoto_data.append(float(items[1]))
                euclidean_data.append(float(items[2]))

        # x軸のダミーデータ
        x_datas = range(len(residue_id_data))

        # グラフのプロット
        fig = plt.figure(figsize=(9.0, 6.0))
        axis1 = fig.add_subplot(111)

        # 左軸：TANIMOTO グラフ
        axis1.plot(x_datas, tanimoto_data,
                   linewidth=1, linestyle='solid', color='blue',
                   label='TANIMOTO')
        axis1.set_ylim(0, 1.2)
        axis1.set_ylabel('TANIMOTO')

        lg.info('Max(TANIMOTO)={}'.format(max(tanimoto_data)))
        lg.info('Min(TANIMOTO)={}'.format(min(tanimoto_data)))

        axis2 = axis1.twinx()
        axis2.plot(x_datas, euclidean_data,
                   linewidth=1, linestyle='solid', color='black',
                   label='EUCLIDEAN')
        axis2.set_ylim(0, (1.2 * max(euclidean_data)) / min(tanimoto_data))
        axis2.set_ylabel('EUCLIDEAN')

        lg.info('Max(euclidean_data)={}'.format(max(euclidean_data)))
        lg.info('Min(euclidean_data)={}'.format(min(euclidean_data)))

        # 目盛りラベルの変更
        # 1. set_xticks で bin の整数倍に目盛りを配置するx座標を指定
        # 2. set_xticklabels で配置された目盛りを文字列リストで置き換え
        visible_xticks = list()
        replace_xlabels = list()
        for idx in range(len(residue_id_data)):
            if(idx % bin == 0):
                visible_xticks.append(idx)
                replace_xlabels.append(residue_id_data[idx])

        axis1.set_xticks(visible_xticks)
        axis1.set_xticklabels(replace_xlabels)

        # # 凡例位置の調整
        hans1, labs1 = axis1.get_legend_handles_labels()
        hans2, labs2 = axis2.get_legend_handles_labels()
        axis1.legend(hans1+hans2, labs1+labs2, bbox_to_anchor=(0.5, -0.1), loc='upper center',
                     borderaxespad=0, ncol=2)

        # 余白調整
        plt.subplots_adjust(bottom=0.20)

        output_file = os.path.join(output_dir, 'similarity_plot.png')
        fig.savefig(output_file)

    except Exception as e:
        lg.error('{}\n{}'.format(e, traceback.format_exc()))
        sys.exit(1)
    finally:
        lg.info("END similarity_plotter")
        for handle in lg.handlers[::-1]:
            handle.close()
            lg.removeHandler(handle)


def _validate(bin: int, region: str):
    """引数に対するバリデーション処理

    Args:
        bin (int): 目盛りの刻み幅
        region (str): 計算対象除外残基領域
    """
    if(type(bin) is not int):
        raise TypeError('[TypeError] bin "{}" is Integer!!'.format(bin))

    if(bin <= 0):
        raise ValueError('[ValueError] bin "{}" less than 1!!'.format(bin))

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


if(__name__ == '__main__'):

    parser = argparse.ArgumentParser(
        description='Read similarity table file and draw a graph.',
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        'similarity_table',
        help='similarity table file',
        type=str
    )
    parser.add_argument(
        '--bin', '-b',
        help='Display width of graph scale. (Defalt 50)',
        type=int,
        default=50
    )
    parser.add_argument(
        '--noregion', '-n',
        help='Exclude residues.\n'
             'Format: (<chain>):(<resnum>) \n'
             'Multiple specifications can be specified by separating them with ","',
        default=''
    )

    args = parser.parse_args()
    table = args.similarity_table
    bin = args.bin
    exclude_region = args.noregion
    main(similarity_table_file=table,
         bin=bin,
         exclude_region=exclude_region,
         output_dir=os.getcwd())
