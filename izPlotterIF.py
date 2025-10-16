#########################################################
# Script to plot averaged I(z) data. 
# Designed for datasets taken at different biases (e.g. 0.5 V and 1 V) and 
# groups of files associated with different sources (adatom, molecule, etc.).
# 
# Multiple groups of files can be defined via the CLI, each with its own 
# source label and z calibration. Data are extracted, binned in z, averaged, 
# and plotted with standard error bands. Outputs include current, conductance, 
# κappa, and dI/dz vs z, with offsets applied to separate curves clearly.
#
# Figures are saved into the `figures/` directory with the provided output stub.
# Optionally, the combined input DataFrame can also be saved as CSV.  
# 
# Command-line structure (Click CLI):
# -----------------------------------
# Top-level options (apply globally):
#   -o / --outstub <str>      : Output filename stub (required).
#   -ocsv / --outcsv          : If set, write input metadata table to figures/<outstub>.csv.
#
# Group subcommands (repeatable, chainable):
#   group -i <files...> -s <source> -c <zcalib>
#
#   -i / --infiles           : One or more input files (wildcards etc. allowed).
#   -s / --source            : Source label, must be one of [adatom | darkNoBromo | molecule].
#   -c / --zcalib <float>    : Z-calibration offset (in pm) added to raw heights.
#
# Multiple `group` commands can be chained in a single execution
# to combine different datasets into one set of averaged plots.
# -------------------------------------------------------
# Example command:
#   python3 izPlotterIFdevelop.py -o izTest group -i data/2024-04-09/multi2_Inj_240409_2_\*.txt -s adatom -c 550 group -i data/2024-04-09/multi2_Inj_240409_9_\*.txt -s adatom -c 700 group -i data/2024-10-22/multi2_Inj_241022_2_\*.txt -i data/2024-10-22/multi2_Inj_241022_\[7-8]_\*.txt -s molecule -c 650 group -i data/2024-10-22/multi2_Inj_241022_\[3-6]_\*.txt -s molecule -c 500
# -------------------------------------------------------
# TODO - refine offsets for multi-source plotting. Is hard coded best way?
# TODO - Add CLI options for binning resolution
# TODO - Option to save combined summary CSV with averaged values and SEM.
#########################################################

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import binned_statistic
import sys
import numpy as np
import matplotlib as mpl
import cycler
from matplotlib.font_manager import FontProperties
import click

import commonFunctions as cf

def exp_function(z, I_0, kappa):
    return I_0*np.exp(-2*kappa*z)

def gaussian(x, *popt):
    A, mu, sigma = popt
    return (A*np.exp(-(x-mu)**2/(2.*sigma**2)))

def offset_gaussian(x, *popt):
    A, mu, sigma, offset = popt
    return (A*np.exp(-(x-mu)**2/(2.*sigma**2))) + offset

def line(x, *popt):
    m, c = popt
    return m*x + c

def binned_mean_stdE(x, y, bins=70):
    """
    Compute the binned mean and standard error of the mean (SEM) for a quantity.
    Bins with no data points will have NaN for both mean and SEM.
    
    Parameters
    ----------
    x : array-like
        Independent variable values used to define the bins (e.g., height).
    y : array-like
        Dependent variable values for which to compute mean and SEM in each bin (e.g., current, conductance).
    bins : int or sequence of floats, optional
        Number of bins or bin edges. Default is 70.

    Returns
    -------
    bin_centers : np.ndarray
        Centers of each bin.
    mean_vals : np.ndarray
        Mean of `y` in each bin.
    sem_vals : np.ndarray
        Standard error on the mean in each bin, calculated as std / sqrt(count).
        Bins with zero counts are returned as np.nan.
    """
    mean_vals, bin_edges, _ = binned_statistic(x, y, statistic='mean', bins=bins)
    std_vals, _, _ = binned_statistic(x, y, statistic='std', bins=bin_edges)
    count_vals, _, _ = binned_statistic(x, y, statistic='count', bins=bin_edges)

    
    # Standard error on the mean
    sem_vals = std_vals / np.sqrt(count_vals)
    # Where counts are zero, set sem to nan
    sem_vals[count_vals == 0] = np.nan 

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    return bin_centers, mean_vals, sem_vals

# --------------------- Begin CLI ---------------------
# Define allowed sources
ALLOWED_SOURCES = ["adatom", "darkNoBromo", "molecule"]

@click.group(chain=True) # declares cli as a Click command group to give acess to group subcommand
                         # chain=True allows multiple group subcommands
# Top level options that apply to all groups
@click.option("--outstub", "-o", help="Out file name stub", type=str, multiple=False, is_eager=True, required=True)
@click.option("--outcsv", "-ocsv", is_flag=True, help="If set, save the combined DataFrame to figures/<outstub>.csv", is_eager=True)

@click.pass_context # Pass the Click contect obj (ctx) which maintains/carries state between the (sub)commands
def cli(ctx, outstub, outcsv):
    """Top-level CLI entry point.

    Initializes the Click context object and sets up shared storage for groups of files, sources, and z-calibration values. 
    This allows multiple `group` subcommands to be accumulated in a single execution.

    Arguments
    ---------
    ctx : click.Context
        The Click context object that persists across subcommands.
    outstub : str
        String to be used for output file names.
    outcsv : str or None
        Path to write the resulting DataFrame as a csv, or None if no file output is required.

    Returns
    -------
    None
        Function modifies the context object in-place.
    """
    ctx.ensure_object(dict) # store input here across subcommands
    ctx.obj["groups"] = [] # collect group from each group subcommand here
    ctx.obj["outstub"] = outstub
    ctx.obj["outcsv"] = outcsv

@cli.command("group") # register the group (cli) subcommand, then define its options
@click.option("--infiles", "-i", type=str, multiple=True, required=True, help="Input file(s) with z lock-in modulated.")
@click.option("--source", "-s", type=click.Choice(ALLOWED_SOURCES, case_sensitive=True), required=True, help=f"Source name (allowed: {', '.join(ALLOWED_SOURCES)}).")
@click.option("--zcalib", "-c", type=float, multiple=False, default=0, help="Calibration height to be added to the raw values, given in pm.")

@click.pass_context # Grabs the Click context object (ctx) for the command handler so it can append to ctx.obj["groups"]
def group_cmd(ctx, infiles, source, zcalib):
    """Define a group of infiles with an associated source and z calibration.

    Expands file specifications and stores the group definition in the Click context object. 
    Multiple groups can be defined in a single command-line invocation.

    Arguments
    ---------
    ctx : click.Context
        The Click context object shared across subcommands.
    infiles : (str)
        Tuple of file specifications (explicit names or wildcards) passed from the command line.
    source : str
        The source identifier. Must be one of the allowed sources.
    zcalib : float
        Z calibration value associated with the infiles, sepcified in pm.

    Returns
    -------
    None
        The function appends a dictionary describing the group to `ctx.obj["groups"]` in the top level cli function.
    """

    all_inFiles = cf.globFiles(infiles)

    if not all_inFiles:
        raise click.ClickException("No files matched the given patterns.")

    ctx.obj["groups"].append( # Adds a dict describing this group to the global groups list in ctx.obj
        {"infiles": all_inFiles, "source": source, "zcalib": zcalib})
    
def collect_groups_to_df(groups, outstub=None, outcsv=False):
    """Helper function for building the pandas df
    Flattens the group definitions into the df and optionally saves to CSV.

    Arguments
    ----------
    groups : list of dict
        Each dict should contain the keys "infiles", "source", and "zcalib".
    outstub : str, optional
        File name stub for the output CSV. Only used if `outcsv` is True.
    outcsv : bool, optional
        If True, write the resulting DataFrame to figures/<outstub>.csv.

    Returns
    -------
    df : pandas.DataFrame
        A DataFrame with one row per file and columns:
        ["inFile", "source", "zcalib"].
    """
    rows = []
    for grp in groups:
        for inF in grp["infiles"]:
            rows.append({"inFile": inF, "source": grp["source"], "zcalib": grp["zcalib"]})

    df = pd.DataFrame(rows)

    # Print nicely
    # click.echo(df.to_string(index=False))

    # Save to CSV if requested
    if outcsv and outstub:
        df.to_csv(f"figures/{outstub}.csv", index=False)
        click.echo(f"\nDataFrame written to: figures/{outstub}.csv")

    return df

# Now to use what we have collected...
@cli.result_callback() # starts w/ decorator to register a function that is called after the chain is finished
@click.pass_context # passes the cxt into the callback
def process_groups(ctx, results, **kwargs):
    """Assemble all defined groups into a df after CLI parsing is complete.

    This callback runs once after all chained `group` subcommands are processed. 
    It collects the group definitions stored in the context to make the df with the helper function, and saves the DataFrame back into the context.
    Also writes to CSV with the helper func if requested. 
    
    Arguments
    ----------
    ctx : click.Context
        The Click context object, with accumulated state from subcommands.
    results : Any
        Return values from the subcommands (unused here).
    **kwargs : dict
        Extra keyword arguments passed automatically by Click (includes outstub and outcsv from the top-level CLI).

    Returns
    -------
    df : pandas.DataFrame
        The combined DataFrame with one row per input file.

    Raises
    ------
    click.ClickException
        If no groups were defined on the command line.
    """
    groups = ctx.obj.get("groups", [])
    outstub = ctx.obj.get("outstub")
    outcsv = ctx.obj.get("outcsv", False)

    if not groups:
        raise click.ClickException("No groups defined. Use `group -i <files> -s <source> -c <zcalib>` at least once.")

    # Now get to the df from the groups w/ the helper func
    df = collect_groups_to_df(groups, outstub, outcsv)  

    # and store the df in ctx.obj so main() can grab it
    ctx.obj["df"] = df  
    return df
# --------------------- End CLI ---------------------


def izPlot(iZfilesdf, outstub):
    """Generate and save averaged i(z) analysis plots for multiple experiments.

    Extracts experimental data from input files, performs averaging with binning
    and standard error estimation, and produces plots of current, conductance,
    κappa, and dI/dz versus z.

    When multiple sources are plotted together, vertical offsets are applied: multiplicative scaling for
    log-scaled plots (current, conductance) and additive offsets for linearly-scaled plots (κappa, dI/dz).

    Arguments
    ---------
    iZfilesdf : pandas.DataFrame
        DataFrame containing file paths and metadata for i(z) experiments.
        Must include columns: 'inFile', 'source', and 'zcalib'.
    outstub : str
        Output filename stub used when saving figures.

    Returns
    -------
    None
        The function writes figures to the `figures/` directory. No objects are returned.
    """

    
    # Get intended bias of the exp and add to file info df
    expBiasList = [] # hold the bias for each experiment
    for inFileZ in iZfilesdf['inFile']:
        try:
            dfZ = cf.extract_file(inFileZ)
        except cf.UnhandledFileError as err:
            print(err)
            continue
        except ValueError as err:
            print(err)
            continue
        
        expBiasList.append(round(dfZ['bias'].mean(), 1))
    iZfilesdf['bias'] = expBiasList
    
    # Generate unique sources and voltages
    sources = sorted(set(iZfilesdf['source']))
    voltages = sorted(set(iZfilesdf['bias']))
    nCombinations = len(sources)*len(voltages)
    
    for source in sources:
        for voltage in voltages:
            mask = (iZfilesdf['source']==source) & (iZfilesdf['bias']==voltage)
            
            print(f'\n### Processing {source} data at {voltage} V. ###')
            
            outerHeight = [] # hold the relative tip height for each point in each experiment
            outerCurrent = [] # hold the current for each point in each experiment
            outerCond = [] # hold the conductance (in G_0) for each point in each experiment
            outerLIX = [] # hold the dIdz for each point in each experiment
            outerKappa = [] # hold the kappa for each point in each experiment
            
            for _, inFileZ in iZfilesdf[mask].iterrows():
                try:
                    dfZ = cf.extract_file(inFileZ['inFile'])
                except cf.UnhandledFileError as err:
                    print(err)
                    continue
                except ValueError as err:
                    print(err)
                    continue
                
                # Check got a sufficient spectrum, not tip withdrawn
                # remove NaNs from the df and check have sufficient data points left
                if len(dfZ['current'][~np.isnan(dfZ['current'])]) < 30:
                    print(f'\n### Not plotted "{inFileZ[0]}". Insufficient data points.\n')
                    continue
                
                # Calibrate z data
                dfZ['height'] = dfZ['height'] + 1e-12*inFileZ['zcalib']
                
                ########### option to switch lockin phase sign if needed
                # dfZ['lIX'] = -dfZ['lIX']
            
                # Mask data to non-noisy region
                quiet = abs(dfZ['current']) > 8e-12
                
                # Get the info for the averaged plots
                outerHeight.extend(dfZ[quiet]['height'].tolist())
                outerCurrent.extend(dfZ[quiet]['current'].to_list())
                outerLIX.extend(dfZ[quiet]['lIX'].to_list())
                # including the conductance - convert from Siemens (1/V * I) to G_0 -> 1 G_0 = 7.748x10^-5 S
                cond = (1/dfZ[quiet]['bias'])*dfZ[quiet]['current']/7.748e-5
                outerCond.extend(cond)
                # and kappa
                kappa = -dfZ[quiet]['lIX']/(2*dfZ[quiet]['current'])
                outerKappa.extend(kappa)
                
            # Then make list of lists into a df
            outerdf = pd.DataFrame(np.column_stack([outerHeight, outerCurrent, outerCond, outerLIX, outerKappa]), columns=['height', 'current', 'conductance', 'lIX', 'kappa'])
            if outerdf.empty:
                print(f'### No valid data for {source} at {voltage} V. Skipping. ###')
                continue
            
            # Compute binned mean and SEM for current, conductance, lIX, kappa
            binnedHeight, aveCurrent, semCurrent = binned_mean_stdE(outerdf['height'], outerdf['current'], bins=70)
            _, aveCond, semCond = binned_mean_stdE(outerdf['height'], outerdf['conductance'], bins=70)
            _, aveLIX, semLIX = binned_mean_stdE(outerdf['height'], outerdf['lIX'], bins=70)
            _, aveKappa, semKappa = binned_mean_stdE(outerdf['height'], outerdf['kappa'], bins=70)

            # Create DataFrame with mean and SEM
            avedf = pd.DataFrame({'binnedHeight': binnedHeight,
                'aveCurrent': aveCurrent, 'semCurrent': semCurrent,
                'aveCond': aveCond, 'semCond': semCond,
                'aveLIX': aveLIX, 'semLIX': semLIX,
                'aveKappa': aveKappa, 'semKappa': semKappa})
            
            # Add vertical offsets for each source...
            # additive for plots on linear scale...
            # and on log scales exponential multiplicative separation as preserves shapes
            source_idx = sources.index(source)
            offset_current = 1.1**source_idx # Current (log scale)
            offset_cond    = 10000** source_idx # Conductance (log scale, gentler separation)
            offset_kappa   = 6*source_idx # κappa (linear scale, mild spread)
            offset_lix     = 2e-12*source_idx # dIdz (linear scale, mild spread)
        
            # Make the averaged plots
            color = plt.cm.plasma(np.linspace(0.1, 0.9, nCombinations))
            mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)
            # starting with current (multi offset)
            plt.figure(1)
            plt.plot(1e12*avedf['binnedHeight'], 1e12*avedf['aveCurrent']*offset_current, marker='o', linestyle='none', markersize=1, label=f'{voltage} V {source}')
            plt.fill_between(1e12*avedf['binnedHeight'], 1e12*(avedf['aveCurrent']-1*avedf['semCurrent'])*offset_current, 1e12*(avedf['aveCurrent']+1*avedf['semCurrent'])*offset_current, alpha=0.2)
            
            # ... and the averaged conductance (multi offset)
            plt.figure(2)
            plt.plot(1e12*avedf['binnedHeight'], avedf['aveCond']*offset_cond, marker='o', linestyle='none', markersize=1, label=f'{voltage} V {source}')
            plt.fill_between(1e12*avedf['binnedHeight'], (avedf['aveCond']-1*avedf['semCond'])*offset_cond, (avedf['aveCond']+1*avedf['semCond'])*offset_cond, alpha=0.2)
            
            # now the averaged kappa (linear offset)
            plt.figure(3, figsize=(6, 5))
            plt.plot(1e12*avedf['binnedHeight'], 1e-9*avedf['aveKappa']+offset_kappa, marker='o', linestyle='none', markersize=1, label=f'{voltage} V {source}')
            plt.fill_between(1e12*avedf['binnedHeight'], 1e-9*(avedf['aveKappa']-1*avedf['semKappa'])+offset_kappa, 1e-9*(avedf['aveKappa']+1*avedf['semKappa'])+offset_kappa, alpha=0.2)
            
            # separating the dIdz (linear offset)
            plt.figure(4)
            plt.plot(1e12*avedf['binnedHeight'], avedf['aveLIX']+offset_lix, marker='o', linestyle='none', markersize=1, label=f'{voltage} V {source}')
            plt.fill_between(1e12*avedf['binnedHeight'], (avedf['aveLIX']-1*avedf['semLIX'])+offset_lix, (avedf['aveLIX']+1*avedf['semLIX'])+offset_lix, alpha=0.2)
            
        # end loop over voltages
    #end loop over sources  
            
    plt.figure(1) # current fig
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel('I (pA)', fontsize=20)
    plt.xlim(250, 890)
    plt.tight_layout()
    plt.yscale('log')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('small')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper right', prop=fontP)
    plt.savefig(f'figures/ave_current_fit_{outstub}.png', bbox_inches='tight')
    plt.close()
    
    plt.figure(2) # conductance fig
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel(r'Conductance ($G_0 = \frac{2e^2}{h}$)', fontsize=20)
    plt.xlim(250, 890)
    plt.tight_layout()
    plt.yscale('log')
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('small')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper right', prop=fontP)
    plt.savefig(f'figures/ave_cond_fit_{outstub}.png', bbox_inches='tight')
    plt.close()
    
    plt.figure(3) # kappa fig
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel(r'$\kappa$ (nm$^{-1}$)', fontsize=20)
    plt.xlim(250, 890)
    plt.tight_layout()
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('small')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='center left', prop=fontP)
    plt.savefig(f'figures/ave_kappa_{outstub}.png', bbox_inches='tight')
    plt.savefig(f'figures/ave_kappa_{outstub}.pdf', bbox_inches='tight')
    plt.close()
    
    plt.figure(4)
    plt.xlabel('z (pm)', fontsize=20)
    plt.ylabel(r'$\partial I / \partial z$', fontsize=20)
    plt.xlim(250, 890)
    plt.tight_layout()
    fontP = FontProperties() # Making legend smaller
    fontP.set_size('small')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='lower right', prop=fontP)
    plt.savefig(f'figures/ave_dIdz_{outstub}.png', bbox_inches='tight')
    plt.close()
    
    return

def main():
    # Prepare a dict to persist CLI context state
    obj = {}
    # Run CLI without exiting, capture context
    iZfilesdf = cli.main(standalone_mode=False, obj=obj)  # stop instant sys.exit from Click and give initial context object 
    if iZfilesdf is None:
        raise RuntimeError("No DataFrame returned from CLI.")
    # Now grab outstub from the copy of the click ctx in main
    outstub = obj.get("outstub")
    if not outstub:
        raise RuntimeError("outstub must be provided via -o/--outstub")
    
    izPlot(iZfilesdf, outstub=outstub)
    
if __name__=='__main__':
 	main()