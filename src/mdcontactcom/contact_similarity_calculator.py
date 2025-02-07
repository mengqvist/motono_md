from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import sys
import traceback
from os.path import join, exists
from log import get_logger


def similarity_main(cp1: str,
                    cp2: str,
                    diff_file: str = '',
                    similarity_file: str = '',
                    log_file: str = ''):
    """contact profile and perform similarity calculation

    Args:
        cp1 (str): contact profile file
        cp2 (str): contact profile file
        diff_file (str): Output destination. Defaults to ''.
        similarity_file (str): Output destination. Defaults to ''.
        log_file (str): Output destination. Defaults to ''.
    """
    # skip if files exist
    if all([exists(f) for f in [diff_file, similarity_file, log_file]]):
        return None

    log_file_path = log_file
    lg = get_logger(name=__name__, log_file_path=log_file_path, is_stdout=True)
    lg.info('START contact_similariry_calculator')

    lg.info(
        "args: \n"
        "contact_profile_1: {cp1}\n"
        "contact_profile_2: {cp2}".format(cp1=cp1,
                                            cp2=cp2)
    )

    lg.info(
        "output: \n"
        "contact profile diff: {diff}\n"
        "similarity table: {s}".format(diff=diff_file,
                                       s=similarity_file)
    )

    try:
        with open(cp1, 'r') as cp1f,\
             open(cp2, 'r') as cp2f,\
             open(diff_file, 'w') as cpdf,\
             open(similarity_file, 'w') as stf:

            header = True
            for cp1_row, cp2_row in zip(cp1f.readlines(), cp2f.readlines()):
                if(header):
                    cpdf.write(cp1_row)
                    stf.write(','.join(['residue_id', 'Tanimoto coefficient',
                                        'Euclidean distance\n']))
                    header = False

                else:
                    residue_id = cp1_row.split(',')[0]
                    cp1_matrix = np.array(cp1_row.split(',')[1:], dtype='float')
                    cp2_matrix = np.array(cp2_row.split(',')[1:], dtype='float')
                    lg.debug("residue id: {id}\n"
                             "cp1: {cp1}\n"
                             "cp2: {cp2}".format(id=residue_id, cp1=cp1_matrix, cp2=cp2_matrix))

                    diff = _calc_contact_profile_diff(list(cp1_matrix), list(cp2_matrix))

                    cpdf.write(",".join([residue_id] + diff) + '\n')

                    tanimoto = _calc_tanimoto(cp1_matrix, cp2_matrix)
                    tanimoto = Decimal(str(tanimoto)).quantize(Decimal('0.001'),
                                                               rounding=ROUND_HALF_UP)

                    euclidian = _calc_euclidian(cp1_matrix, cp2_matrix)
                    euclidian = Decimal(str(euclidian)).quantize(Decimal('0.001'),
                                                                 rounding=ROUND_HALF_UP)

                    stf.write(','.join([residue_id, str(tanimoto), str(euclidian)]) + '\n')

    except Exception as e:
        lg.error('{}\n{}'.format(e, traceback.format_exc()))
        sys.exit(1)
    finally:
        lg.info("END contact_similariry_calculator")
        for handle in lg.handlers[::-1]:
            handle.close()
            lg.removeHandler(handle)


def _calc_contact_profile_diff(x: list, y: list) -> list:
    """1次元のcontact profile map A function that computes the difference array for

    Args:
        x (list): contact profile map
        y (list): contact profile map

    contact profile map: One-dimensional matrix Each element is 'float'

    Returns:
        list: contact profile map difference matrix of
              One-dimensional matrix Each element is 'str'
    """
    assert type(x) is list, 'Not List!!: {}'.format(x)
    assert type(y) is list, 'Not List!!: {}'.format(y)
    assert len(x) == len(y)

    contact_profile_diff_map = list()
    for idx in range(len(x)):

        x_data = Decimal(str(x[idx]))
        y_data = Decimal(str(y[idx]))

        if(x_data == 0 and y_data == 0):
            # both contact profile is the value of 0 In the case of
            # The value of the difference is regarded as meaningless and is treated as an empty string.。
            contact_profile_diff_map.append('')

        else:
            contact_profile_diff_map.append(str(x_data - y_data))

    return contact_profile_diff_map


def _calc_tanimoto(x: list, y: list) -> float:
    """A function that calculates the TANIMOTO coefficient from a two-array

    Args:
        x (list/numpy.ndarray): one-dimensional array
        y (list/numpy.ndarray): one-dimensional array

    Returns:
        float: Arrangement x, y calculated from TANIMOTO coefficient
    """
    assert len(x) == len(y), 'len(x)={}, len(y)={}'.format(len(x), len(y))
    if(type(x) is list):
        x = np.array(x, dtype='float')
    if(type(y) is list):
        y = np.array(y, dtype='float')

    assert type(x) is np.ndarray, 'type(x): {}'.format(type(x))
    assert type(y) is np.ndarray, 'type(y): {}'.format(type(y))

    x_sum_square = (x ** 2).sum()
    y_sum_square = (y ** 2).sum()
    x_y_sum = (x * y).sum()

    denominator = (x_sum_square + y_sum_square - x_y_sum)
    if (denominator == 0):
        return 0
    else:
        return x_y_sum / denominator


def _calc_euclidian(x: list, y: list) -> float:
    """A function that computes the Euclidean distance between two points in n dimensions

    Args:
        x (list/numpy.ndarray): `v`
        y (list/numpy.ndarray): 1D array representing coordinates

        *x, y Each element of indicates the value of each axis

    Returns:
        float: between two points x, y Euclidean distance of
    """
    assert len(x) == len(y), 'len(x)={}, len(y)={}'.format(len(x), len(y))
    if(type(x) is list):
        x = np.array(x, dtype='float')
    if(type(y) is list):
        y = np.array(y, dtype='float')

    assert type(x) is np.ndarray, 'type(x): {}'.format(type(x))
    assert type(y) is np.ndarray, 'type(y): {}'.format(type(y))

    return np.sqrt(((x - y) ** 2).sum())

#
# if(__name__ == '__main__'):
#
#     parser = argparse.ArgumentParser(
#         description='Read contact profile files and calculate the similarity.',
#         formatter_class=RawTextHelpFormatter
#     )
#     parser.add_argument(
#         'contact_profile_1',
#         help='contact profile file',
#         type=str
#     )
#     parser.add_argument(
#         'contact_profile_2',
#         help='contact profile file',
#         type=str
#     )
#
#     args = parser.parse_args()
#     cp1 = args.contact_profile_1
#     cp2 = args.contact_profile_2
#     similarity_main(cp1, cp2)
