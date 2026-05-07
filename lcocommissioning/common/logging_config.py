"""
Centralized logging configuration for lcocommissioning package.

Eliminates repeated logging setup across all CLI modules.
"""

import logging
import sys


def setup_logging(name=None, level='INFO', format_style='detailed'):
    """
    Configure logging with standardized format and suppressed third-party loggers.

    Parameters
    ----------
    name : str, optional
        Logger name (typically __name__). If None, uses root logger.
    level : str, optional
        Logging level: 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
        Default: 'INFO'
    format_style : str, optional
        Format style: 'detailed', 'simple', or 'minimal'
        Default: 'detailed'

    Returns
    -------
    logging.Logger
        Configured logger instance

    Examples
    --------
    >>> log = setup_logging(__name__)
    >>> log.info("Starting observation submission")

    >>> log = setup_logging(__name__, level='DEBUG')
    """
    formats = {
        'detailed': '%(asctime)s.%(msecs)03d %(levelname)7s: %(module)20s: %(message)s',
        'simple': '%(levelname)s - %(name)s - %(message)s',
        'minimal': '%(levelname)s: %(message)s',
    }

    log_format = formats.get(format_style, formats['detailed'])
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format=log_format,
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Suppress verbose third-party loggers
    suppress_loggers = [
        'matplotlib',
        'matplotlib.font_manager',
        'opensearch',
        'connectionpool',
        'urllib3',
        'botocore',
        'boto3',
    ]

    for logger_name in suppress_loggers:
        logging.getLogger(logger_name).setLevel(logging.WARNING)

    return logging.getLogger(name)


def get_logger(name):
    """
    Get logger without reconfiguring (for use after initial setup).

    Parameters
    ----------
    name : str
        Logger name (typically __name__)

    Returns
    -------
    logging.Logger
        Logger instance
    """
    return logging.getLogger(name)
