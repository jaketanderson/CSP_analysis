
## Summary
This repository contains a script for analyzing and plotting chemical shift perturbations in 2D NMR.

## Usage
### Install
This script only requires Python, Matplotlib and Pandas. Once you have Python you can install these dependencies using `pip install matplotlib pandas`.

Additionally, I recommend installing `colorlog` for prettier logging.
This script was tested on Python 3.14.4 but any modern Python version should work.

### Run
Use `python CSP_analysis.py`. For help and available arguments, you can use `python CSP_analysis --help` or you can refer to the [reference](#reference) section below.

Running the script will, by default, produce an analysis folder named `analysis_<date>_<time>`. Inside, there will be versions of the peak lists (with extra columns, including calculated CSPs), plots by residue/nucleus, plots by residue/nucleus/titration, and a log file called `CSP_analysis.log`.

Ensure the peak lists you provide match the format of the example peak lists; e.g. do not use lists with assignment columns that have different formats for residue name, residue number, and nuclei.
This script will work with or without an intensities column; it does _not_ make use of it if it is given.

### Customize
The file `CSP_analysis.py` contains a dictionary at the top called `CUSTOM_ALPHAS` that can be edited to specify alpha values for particular residues and nuclei. These values will overwrite any values automatically derived from `shift_database.csv`.

## Files
|File| Purpose |
|--|--|
| `CSP_analysis.py` | The main script. Run this directly. |
| `utils.py` | A helper script to hold some of the analysis functions. |
| `[1-5].list` | Example peak lists with short names, for testing. |
| `Multi_analysis_CSP_INT_Sparky_7_NA.m` | The original MATLAB script from DF. This repository intends to replicate and extend its functionality in Python. |
| `shift_database.csv` | A database of chemical shifts. The `std` entries are used in ratios to derive alpha values. Downloaded June 26, 2026 from [BMRB](https://bmrb.io/ref_info/csstats.php?set=full&restype=aa&output=csv). |

## Reference

    usage: CSP_analysis.py [-h] [--bad-residues BAD_RESIDUES] [--residue-offset RESIDUE_OFFSET] [--log LOG] [--no-log] [--dir DIR]
                           [--plot-format PLOT_FORMAT] [--plot-percentiles PLOT_PERCENTILES] [--scale-plots]
    
    CSP_analysis is a script to analyze chemical shift perturbations across a titration.
    
    options:
      -h, --help            show this help message and exit
      --bad-residues BAD_RESIDUES
                            A list of residues for which the peaks are bad and ought to be excluded from analysis. E.g. "9,15,20".
      --residue-offset RESIDUE_OFFSET
                            The starting index will be shifted by this amount.
      --log LOG             The file to which you'd like to save logs from this script. It will live inside the directory specified by
                            --dir. Set --nolog to disable.
      --no-log
      --dir DIR             The directory to which you'd like to save logs and plots from this script.
      --plot-format PLOT_FORMAT
                            The file format you'd like the plots to be given in. Valid options are "png", "svg", and "pdf".
      --plot-percentiles PLOT_PERCENTILES
                            A list of percentiles you'd like to plot. The highest will be colored red. E.g. "75,90"
      --scale-plots         This flag will scale up all plots to each take up the entire y-axis. Otherwise, by default, all plots will
                            have the same y-axis.

