import logging
import os
import sys

"""ログに関するモジュール
"""


def get_logger(name: str, log_file_path: str = "", is_stdout: bool = True) -> logging.Logger:
    """ロガーを返す関数

    Args:
        name (str): ロガーの名前
        log_file_path (str, optional): ログファイル. Defaults to "".
        is_stdout (bool, optional): 標準出力フラグ. Defaults to True.

    Returns:
        logging.Logger: ログレベル、メッセージフォーマット等が適切に設定されたロガー
    """
    logger_log_level = logging.INFO
    stream_log_level = logging.WARNING
    file_log_level = logging.INFO
    if(sys.flags.debug):
        # デバッグモード(-d オプションによる実行)時のログレベル
        logger_log_level = logging.DEBUG
        stream_log_level = logging.DEBUG
        file_log_level = logging.DEBUG

    logger = logging.getLogger(name)
    logger.setLevel(logger_log_level)

    # https://docs.python.org/ja/3.6/library/logging.html#logrecord-attributes
    formatter = logging.Formatter(
        '%(asctime)s - %(filename)s - Func: %(funcName)-25s - %(levelname)-8s - %(message)s')

    # ハンドラ（ログの送信先）の設定
    if(is_stdout):
        # ストリームハンドラの追加（標準出力）
        sh = logging.StreamHandler()
        sh.setLevel(stream_log_level)
        sh.setFormatter(formatter)
        logger.addHandler(sh)

    if(len(log_file_path) != 0):
        if(len(os.path.dirname(log_file_path)) != 0):
            output_dir = os.path.dirname(log_file_path)
            assert os.path.exists(output_dir), "Not Found '{}'".format(output_dir)

        # ファイルハンドラの追加（ファイル出力）
        fh = logging.FileHandler(log_file_path, mode='w', encoding='utf-8')
        fh.setLevel(file_log_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger
