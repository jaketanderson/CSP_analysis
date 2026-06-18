import argparse
import logging
import os
import sys
from datetime import datetime

import colorlog
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Chemical shift perturbations, from Eq 8 of doi.org/10.1016/j.pnmrs.2013.02.001
def CSP(w1: float, w2: float, alpha1: float, alpha2: float):
    return np.sqrt(0.5 * ((alpha1 * w1) ** 2 + (alpha2 * w2) ** 2))


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


def get_consistent_max_residue(peaklists: List[pd.DataFrame], logger) -> (int, float):
    max_residue_indices = []
    max_vals = []
    for peaklist in peaklists[1:]:
        max_residue_indices.append(peaklist["CSP"].idxmax())
        max_vals.append(peaklist["CSP"].max())

    max_residue_indices = set(max_residue_indices)
    if len(max_residue_indices) == 1:
        logger.info(
            f"Residue {peaklist.loc[list(max_residue_indices)[0]]['Assignment']} (index {list(max_residue_indices)[0]}) has the maximum CSP for all titrants."
        )
        index = list(max_residue_indices)[0]
        value = max(max_vals)
        return index, value
    else:
        logger.info(
            "The residue with the maximum CSP is different across the titration. No autoscaling possible"
        )
        return None, None


def main():
    script_start_time_str = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    parser = argparse.ArgumentParser(
        description="CSP_analysis is a script to analyze chemical shift perturbations across a titration."
    )

    parser.add_argument(
        "--bad-residues",
        type=str,
        default="",
        help='A list of residues for which the peaks are bad and ought to be excluded from analysis. E.g. "9,15,20".',
    )
    parser.add_argument(
        "--residue-offset",
        type=int,
        default=0,
        help="The starting index will be shifted by this amount.",
    )
    parser.add_argument(
        "--log",
        type=str,
        default="CSP_analysis.log",
        help="The file to which you'd like to save logs from this script. It will live inside the directory specified by --dir. Set --nolog to disable.",
    )
    parser.add_argument("--no-log", action="store_true")
    parser.add_argument(
        "--dir",
        type=str,
        default=f"analysis_{script_start_time_str}",
        help="The directory to which you'd like to save logs and plots from this script.",
    )
    parser.add_argument(
        "--plot-format",
        type=str,
        default="pdf",
        help='The file format you\'d like the plots to be given in. Valid options are "png", "svg", and "pdf".',
    )
    args = parser.parse_args()

    # Set up logging
    if not os.path.exists(args.dir):
        os.mkdir(args.dir)
    else:
        raise AssertionError(f"The analysis directory {args.dir} already exists!")
    logfile = os.path.abspath(f"{args.dir}/{args.log}")
    logger = setup_logging(args.no_log, logfile)
    logger.info("You ran: " + str(sys.executable) + " " + os.path.abspath(__file__))
    logger.info(f"This log file is located at {logfile}")

    # Check that the filetype is valid
    if args.plot_format not in ("png", "svg", "pdf"):
        logger.error(f'File format "{args.plot_format}" is invalid!')
        raise AssertionError(f'File format "{args.plot_format}" is invalid!')

    badresidues = [
        badres for badres in args.bad_residues.split(",") if badres.isdigit()
    ]

    peakfiles = []
    ctr = 0
    while "done" not in [f.lower() for f in peakfiles]:
        desc = "first" if ctr == 0 else "next"
        peakfiles.append(
            input(f'Give the {desc} peak file. When finished, write "done": ')
        )
        ctr += 1
    # Remove the "done" since it's not a peakfile
    peakfiles.pop(-1)
    # Remove any peakfiles that were from just pressing enter
    peakfiles = [
        os.path.abspath(peakfile) for peakfile in peakfiles if peakfile.strip() != ""
    ]
    logger.info(f"Using peakfiles:\n  {'\n  '.join(peakfiles)}")

    titration_percs = []
    peaklists = []
    for peakfile in peakfiles:

        df = pd.read_csv(peakfile, sep=r"\s{2,}", engine="python")
        df["ResidueType"] = pd.Series([a[:1] for a in df["Assignment"]])
        df["ResidueIndex"] = pd.Series([a.split("_")[0][1:] for a in df["Assignment"]])
        df["Nucleus1"] = pd.Series(
            [a.split("_")[1].split("-")[0] for a in df["Assignment"]]
        )
        df["Nucleus2"] = pd.Series(
            [a.split("_")[1].split("-")[1] for a in df["Assignment"]]
        )

        # Remove bad residues
        df = df[df["ResidueIndex"].isin(badresidues) == False]

        peaklists.append(df)
        titration_percs.append(
            float(
                input(
                    f"Give the titration percentage for peak file {peakfile}: "
                ).strip()
            )
        )

    # Sort the peaklists so they are in ascending order of titration percentage
    peaklists = [
        pl
        for pl, _ in sorted(
            zip(peaklists, titration_percs), key=lambda zipped: zipped[1]
        )
    ]
    titration_percs.sort()

    if titration_percs[0] != 0:
        logger.warning(
            "The lowest titration percentage is not zero. Did you make a mistake? It will still be used as a reference... be careful!"
        )
    logger.info(f"Using titration percentages {titration_percs}")

    assert (
        len(set([len(peaklist) for peaklist in peaklists])) == 1
    ), "The number of peaks per file is different. Each file must have the same number of peaks after bad-residue exclusion."

    unique_nonH_nuclei = set(
        [x for peaklist in peaklists for x in peaklist["Nucleus1"].to_list()]
    )
    unique_nonH_nuclei |= set(
        [x for peaklist in peaklists for x in peaklist["Nucleus2"].to_list()]
    )
    unique_nonH_nuclei.discard("H")
    alphas = {"H": 1.0}
    for nucleus in unique_nonH_nuclei:
        alphas[nucleus] = float(
            input(
                f"Give the alpha for scaling perturbations in {nucleus} relative to H: "
            ).strip()
        )
    logger.info(f"Using alphas: {alphas}")

    for peaklist in peaklists:
        # We are comparing this peaklist to the first peaklist, which has the lowest titration value
        zipped = list(
            zip(peaklist.itertuples(index=True), peaklists[0].itertuples(index=True))
        )
        assert all(
            [
                (rowi.Nucleus1 == rowj.Nucleus1) and (rowi.Nucleus2 == rowj.Nucleus2)
                for rowi, rowj in zipped
            ]
        ), "There are peaks with mismatched nuclei from peaklist to peaklist! Please review the peaklists line-by-line and ensure they match. E.g. N-H with N-H, so on."

        peaklist["CSP"] = pd.Series(
            [
                CSP(
                    rowj.w1 - rowi.w1,
                    rowj.w2 - rowi.w2,
                    alphas[rowi.Nucleus1],
                    alphas[rowi.Nucleus2],
                )
                for rowi, rowj in zipped
            ]
        )

    # Perform autoscaling to a consistent residue's maximum CSP
    max_residue_index, max_value = get_consistent_max_residue(peaklists, logger)
    if max_residue_index:
        while True:
            autoscale = (
                input(
                    "Would you like to autoscale so this residue's maximum is at the top of each plot? (y/n): "
                )
                .lower()
                .strip()
            )
            if autoscale in ("y", "yes"):
                autoscale = True
                break
            elif autoscale in ("n", "no"):
                autoscale = False
                break

    logger.info(f"Creating {args.plot_format.upper()}s")
    # Include the full log in the metadata of the image
    with open(f"{args.dir}/{args.log}", "r") as f:
        full_log = f.read()
    plt.rcParams.update({"font.size": 14})
    for nucleus in unique_nonH_nuclei:
        for i, peaklist in enumerate(
            peaklists[1:]
        ):  # We do [1:] because first peaklist is baseline and the CSPs=0
            fig = plt.figure(figsize=(16, 9))
            plt.title(f"H{nucleus} at {titration_percs[i+1]}%")
            peaks_filtered = peaklist[(peaklist["Nucleus1"] == nucleus)]
            CSPs_np = peaks_filtered.CSP.to_numpy()
            mean = np.mean(CSPs_np)
            median = np.median(CSPs_np)
            stdev = np.std(CSPs_np)

            colors = []
            for row in peaks_filtered.itertuples(index=True):
                if row.CSP > mean + stdev:
                    colors.append("red")
                else:
                    colors.append("blue")

            plt.bar(
                peaks_filtered.ResidueIndex.to_numpy(dtype=int),
                CSPs_np,
                color=colors,
                width=1.0,
            )
            plt.axhline(
                y=mean,
                color="black",
                linestyle="-",
                label=r"$\overline{\Delta \delta}$",
            )
            plt.axhline(
                y=mean + stdev,
                color="black",
                linestyle="--",
                label=r"$\overline{\Delta \delta} \pm \sigma$",
            )
            plt.axhline(
                y=mean - stdev,
                color="black",
                linestyle="--",
            )

            plt.scatter(
                min(peaks_filtered.ResidueIndex.to_numpy(dtype=int)),
                0,
                label=" ",
                visible=False,
            )
            plt.scatter(
                min(peaks_filtered.ResidueIndex.to_numpy(dtype=int)),
                0,
                label=r"$\overline{\Delta \delta}$"
                + f" = {mean:0.3f} ppm,\n"
                + r"$\sigma$"
                + f"   = {stdev:0.3f} ppm",
                visible=False,
            )

            # For labeling the residues on the x axis. Maybe allow as an option later.
            # plt.xticks(
            #     peaks_filtered.ResidueIndex.to_numpy(dtype=int),
            #     labels=peaks_filtered["Assignment"].str.split("_").str[0].tolist(),
            #     fontsize=8,
            #     rotation=45,
            # )
            plt.xticks(
                np.arange(0, max(peaks_filtered.ResidueIndex.to_numpy(dtype=int)), 5),
            )
            plt.xlim(
                xmin=peaks_filtered.ResidueIndex.to_numpy(dtype=int).min() - 0.5,
                xmax=peaks_filtered.ResidueIndex.to_numpy(dtype=int).max() + 0.5,
            )
            plt.xlabel("Residues")
            plt.ylabel(r"$\Delta \delta$" + " (ppm)")

            if autoscale:
                plt.ylim(0, max_value)
            plt.legend(framealpha=0.75)
            plt.tight_layout()

            plt.savefig(
                f"{args.dir}/H{nucleus}_{titration_percs[i+1]}.{args.plot_format}",
                metadata={"Title": full_log},  # Title should work across all 3 formats
            )


if __name__ == "__main__":
    main()
