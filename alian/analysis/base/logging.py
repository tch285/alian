import logging
import os


class ColoredFormatter(logging.Formatter):
    COLORS = {
        'WARNING': '\033[33m',
        'ERROR': '\033[31m',
        'DEBUG': '\033[34m',
        'INFO': '\033[32m',
        'CRITICAL': '\033[35m'
    }
    RESET = '\033[0m'

    def format(self, record):
        color = self.COLORS.get(record.levelname, '')
        if color:
            # Color the entire line
            formatted_msg = super().format(record)
            return f"{color}{formatted_msg}{self.RESET}"
        return super().format(record)

def setup_logger(name, handler_level = logging.INFO, logger_level = logging.DEBUG,
                 log_fmt = '%(asctime)s - %(levelname)s - %(filename)s:%(lineno)d - %(funcName)s - %(message)s',
                 date_fmt = '%Y-%m-%d %H:%M:%S'
                 ):
    logger = logging.getLogger(name)
    handler = logging.StreamHandler()

    # Create a formatter and set it for the handler
    if "SLURM_SUBMIT_DIR" in os.environ:
        # use uncolored Formatter for Slurm jobs
        formatter = logging.Formatter(log_fmt, datefmt = date_fmt)
    else:
        # use colored Formatter for local testing
        formatter = ColoredFormatter(log_fmt, datefmt = date_fmt)

    handler.setFormatter(formatter)
    handler.setLevel(handler_level)
    logger.addHandler(handler)
    logger.setLevel(logger_level)
    return logger