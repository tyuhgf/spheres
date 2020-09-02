import logging
import tempfile

LOG_LEVEL = logging.DEBUG
log_formatter = logging.Formatter('[%(asctime)s - %(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s')
log_handler = logging.StreamHandler()
log_handler.setLevel(LOG_LEVEL)
log_handler.setFormatter(log_formatter)

TMP_DIR = tempfile.gettempdir()
