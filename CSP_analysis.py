import argparse
import logging
import os
import pprint
import sys
from datetime import datetime

import colorlog
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import utils

# Use this dictionary to overwrite the alphas from the shift database with custom values.
CUSTOM_ALPHAS = {
    #"GLYHA3": 4.3333,
}

shift_db = pd.read_csv("shift_database.csv", sep=",")

# Chemical shift perturbation calcs
# Note that this is DIFFERENT from eq 8 of doi.org/10.1016/j.pnmrs.2013.02.001
def CSP(w1: float, w2: float, alpha1: float, alpha2: float):
    return np.sqrt(0.5 * ((alpha1 * w1) ** 2 + (alpha2 * w2) ** 2))


def get_max_CSP(peaklists: List[pd.DataFrame], logger) -> float:
    max_vals = []
    for peaklist in peaklists[1:]:
        max_vals.append(peaklist["CSP"].max())

    return max(max_vals)


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
    parser.add_argument(
        "--plot-percentiles",
        type=str,
        default="75,90",
        help='A list of percentiles you\'d like to plot. The highest will be colored red. E.g. "75,90"',
    )
    parser.add_argument(
        "--scale-plots",
        action="store_true",
        default=False,
        help="This flag will scale up all plots to each take up the entire y-axis. Otherwise, by default, all plots will have the same y-axis.",
    )
    args = parser.parse_args()

    # Set up logging
    if not os.path.exists(args.dir):
        os.mkdir(args.dir)
    else:
        raise AssertionError(f"The analysis directory {args.dir} already exists!")
    logfile = os.path.abspath(f"{args.dir}/{args.log}")
    logger = utils.setup_logging(args.no_log, logfile)
    logger.info("You ran: " + str(sys.executable) + " " + os.path.abspath(__file__))
    logger.info(f"This log file is located at {logfile}")
    logger.info("Using the following args:\n" + pprint.pformat(vars(args)) + "\n")

    utils.check_args(args, logger)

    badresidues = [
        badres.strip().strip("\"'").strip() for badres in args.bad_residues.split(",")
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

    alphas = {}
    # Use the lowest stdev from all residue/nucleus combos as the reference
    reference = min(shift_db["std"])
    # Populate alphas as ratios of stdev to reference
    for row in shift_db.itertuples():
        key = row.comp_id + row.atom_id
        alphas[key] = row.std / reference

    for peaklist, peakfile in list(zip(peaklists, peakfiles)):
        # We are comparing this peaklist to the first peaklist, which has the lowest titration value
        zipped = list(
            zip(peaklist.itertuples(index=True), peaklists[0].itertuples(index=True))
        )
        assert all(
            [
                (rowi.ResidueIndex == rowj.ResidueIndex) and
                (rowi.ResidueType == rowj.ResidueType) and
                (rowi.Nucleus1 == rowj.Nucleus1) and
                (rowi.Nucleus2 == rowj.Nucleus2)
                for rowi, rowj in zipped
            ]
        ), "There are peaks with mismatched nuclei/residues from peaklist to peaklist! Please review the peaklists line-by-line and ensure they match. E.g. N-H with N-H, so on."

        peaklist["CSP"] = pd.Series(
            [
                CSP(
                    rowj.w1 - rowi.w1,
                    rowj.w2 - rowi.w2,
                    utils.get_alphas(alphas, utils.AA[rowi.ResidueType]+rowi.Nucleus1),
                    utils.get_alphas(alphas, utils.AA[rowi.ResidueType]+rowi.Nucleus2),
                )
                for rowi, rowj in zipped
            ]
        )
        peakfile_trimmed = peakfile.split('/')[-1].split(".", 1)
        peaklist.to_csv(f"{args.dir}/{peakfile_trimmed[0]}_CSPs.{peakfile_trimmed[1]}", index=False)
    used_alphas = {key: round(alphas[key], 3) for key in utils.used_alpha_keys if key in alphas}
    logger.info(f"List of alphas relevant to your residues:\n{used_alphas}")


    logger.info(f"Creating {args.plot_format.upper()}s")
    # Include the full log in the metadata of the image
    with open(f"{args.dir}/{args.log}", "r") as f:
        full_log = f.read()
    plt.rcParams.update({"font.size": 14})
    plot_percentiles = [
        float(plot_perc.strip().strip("\"'").strip())
        for plot_perc in args.plot_percentiles.split(",")
    ]
    max_CSP = get_max_CSP(peaklists, logger)
    
    # Choose subplot layout based on number of peaklists
    # This is naïve, but I think it is the most readable option
    l = len(peaklists) - 1
    if 1 <= l <= 3:
        r, c = (1, l)
    elif l == 4:
        r, c = (2, 2)
    elif 5 <= l <= 6:
        r, c = (2, 3)
    else:
        r, c = (round(l/3)+1, 3)
    
    # Iterate over nucleus types and peaklists (titration percentages)
    nucleus_pairs = list(zip(peaklists[0]["Nucleus1"], peaklists[0]["Nucleus2"]))
    for nucleus_pair in set(nucleus_pairs):
        bigfig, axes = plt.subplots(nrows=r, ncols=c, figsize=(8*r,8*c))
        for i, peaklist in enumerate(
            peaklists[1:]
        ):  # We do [1:] because first peaklist is baseline and the CSPs=0
            peaks_filtered = peaklist[
                (peaklist["Nucleus1"] == nucleus_pair[0]) & (peaklist["Nucleus2"] == nucleus_pair[1])
            ]
            if len(peaks_filtered) < 2:
                continue
            fig = plt.figure(figsize=(16, 9))
            plt.title(f"{nucleus_pair[0]}{nucleus_pair[1]} at {titration_percs[i+1]}%")
            if l == 1:
                axes.flat = [axes] # This is needed because a 1x1 Axes is not a list of Axes by default
            axes.flat[i].set_title(f"{nucleus_pair[0]}{nucleus_pair[1]} at {titration_percs[i+1]}%")
            CSPs_np = peaks_filtered.CSP.to_numpy()
            mean = np.mean(CSPs_np)
            median = np.median(CSPs_np)
            stdev = np.std(CSPs_np)
            plot_perc_values = [np.percentile(CSPs_np, q) for q in plot_percentiles]
            plot_perc_values.sort(reverse=True)
            plot_percentiles.sort(reverse=True)

            colors = []
            for row in peaks_filtered.itertuples(index=True):
                if row.CSP > plot_perc_values[0]:
                    colors.append("red")
                else:
                    colors.append("blue")

            plt.bar(
                peaks_filtered.ResidueIndex.to_numpy(dtype=int),
                CSPs_np,
                color=colors,
                edgecolor="black",
                width=1.0,
            )
            axes.flat[i].bar(
                peaks_filtered.ResidueIndex.to_numpy(dtype=int),
                CSPs_np,
                color=colors,
                edgecolor="black",
                width=1.0,
            )

            linestyles = ["dashed", "dashdot", "dotted"][: len(plot_perc_values)]
            for perc, perc_val, ls in list(
                zip(plot_percentiles, plot_perc_values, linestyles)
            ):
                plt.axhline(
                    y=perc_val,
                    color="black",
                    linestyle=ls,
                    label=f"{perc}%     = {perc_val:0.4f}",
                )
                axes.flat[i].axhline(
                    y=perc_val,
                    color="black",
                    linestyle=ls,
                    label=f"{perc}%     = {perc_val:0.4f}",
                )

            plt.axhline(
                y=median,
                color="black",
                linestyle="-",
                label=r"$\operatorname{med}(\Delta \delta)$" + f" = {median:0.4f}",
            )
            axes.flat[i].axhline(
                y=median,
                color="black",
                linestyle="-",
                label=r"$\operatorname{med}(\Delta \delta)$" + f" = {median:0.4f}",
            )

            plt.xticks(
                np.arange(0, max(peaks_filtered.ResidueIndex.to_numpy(dtype=int)), 5),
            )
            axes.flat[i].set_xticks(
                np.arange(0, max(peaks_filtered.ResidueIndex.to_numpy(dtype=int)), 5),
            )
            plt.xlim(
                xmin=peaks_filtered.ResidueIndex.to_numpy(dtype=int).min() - 0.5,
                xmax=peaks_filtered.ResidueIndex.to_numpy(dtype=int).max() + 0.5,
            )
            axes.flat[i].set_xlim(
                xmin=peaks_filtered.ResidueIndex.to_numpy(dtype=int).min() - 0.5,
                xmax=peaks_filtered.ResidueIndex.to_numpy(dtype=int).max() + 0.5,
            )
            plt.xlabel("Residues")
            axes.flat[i].set_xlabel("Residues")
            plt.ylabel(r"$\Delta \delta$" + " (ppm)")
            axes.flat[i].set_ylabel(r"$\Delta \delta$" + " (ppm)")

            if not args.scale_plots:
                plt.ylim(0, max_CSP)
                axes.flat[i].set_ylim(0, max_CSP)
            plt.legend(framealpha=0.75, loc="upper left")
            axes.flat[i].legend(framealpha=0.75, loc="upper left")
            plt.tight_layout()

            fig.savefig(
                f"{args.dir}/{nucleus_pair[0]}{nucleus_pair[1]}_{titration_percs[i+1]}.{args.plot_format}",
                metadata={"Title": full_log},  # Title should work across all 3 formats
            )
            
        bigfig.savefig(
            f"{args.dir}/{nucleus_pair[0]}{nucleus_pair[1]}.{args.plot_format}",
            metadata={"Title": full_log},
        )


if __name__ == "__main__":
    main()
