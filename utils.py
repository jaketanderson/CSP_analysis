import argparse
import logging
import sys

import colorlog


def setup_logging(no_log: bool, log: str):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s- %(levelname)s - %(message)s")

    stdout_handler = colorlog.StreamHandler(sys.stdout)
    stdout_handler.setFormatter(
        colorlog.ColoredFormatter(
            "%(log_color)s%(levelname)s%(reset)s - %(message)s",
            log_colors={
                "DEBUG": "cyan",
                "INFO": "green",
                "WARNING": "orange",
                "ERROR": "red",
            },
        )
    )

    logger.addHandler(stdout_handler)

    if not no_log:
        file_handler = logging.FileHandler(log)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def check_args(args, logger) -> None:
    # Check that the plot filetype is valid
    if args.plot_format not in ("png", "svg", "pdf"):
        logger.error(f'File format "{args.plot_format}" is invalid!')
        raise AssertionError(f'File format "{args.plot_format}" is invalid!')

    # Check that badresidues is a valid list of ints
    if not all(
        [badres.strip().strip("\"'").strip().isdigit() for badres in args.bad_residues]
    ):
        logger.error(
            f'Bad residues list "{args.bad_residues}" contained empty or non-integer values!'
        )
        raise AssertionError(
            f'Bad residues list "{args.bad_residues}" contained empty or non-integer values!'
        )

    # Check that badresidues is a valid list of floats
    try:
        plot_percentiles = [
            float(plot_perc.strip().strip("\"'").strip())
            for plot_perc in args.plot_percentiles.split(",")
        ]
    except:
        logger.error(
            f'Plot percentile list "{args.plot_percentiles}" contained empty or non-float values!'
        )
        raise AssertionError(
            f'Plot percentile list "{args.plot_percentiles}" contained empty or non-float values!'
        )
    if len(plot_percentiles) > 3:
        logger.error(f"Plot percentile list is too long. Must be 3 or fewer.")
        raise AssertionError(f"Plot percentile list is too long. Must be 3 or fewer.")
