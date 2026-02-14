#Author: Louise Dupont (with some additions by John)
import matplotlib.pyplot as plt
# Explicitly import the SVG backend so tools like Nuitka bundle it correctly
from matplotlib.backends import backend_svg  # noqa: F401
import seaborn as sns
import pandas as pd
import numpy as np
from textwrap import wrap
import warnings
# silence noisy FutureWarnings & generic matplotlib UserWarnings
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")
warnings.filterwarnings("ignore", category=UserWarning,  module="matplotlib")

##matplotlib.use('svg')

def create_sps_plot(csv_file_path=None,
                    df=None, RMSE_rank = 0.05,
                    bw_method = 0.07,
                    fig_path = 'plot.png',
                    title = 'Model Predictions',
                    neg_pheno_name = "Negative",
                    pos_pheno_name = "Positive",
                    neg_pheno_color = '#F55E54',
                    pos_pheno_color = '#2fc8cc',
                    percent_accuracy = True,
                    axes = None,
                    min_genes = 0):
    '''creates density plot of sequence prediction score for each species in each ESL run 
    plot has two lines: one for the negative phenotype and one for the positive phenotype
    can either pass csv file or dataframe but not both
    required columns in data: SPS values, true phenotype, input RMSE
    
    Args:
        param1: csv file path to create chart from
        param2: dataframe to create chart from
        param3: RMSE percentage rank to be included in chart (0.05 = lowest 5 percent)
        param4: bw_method (smoothing) used to make plot, higher numbers smooth plot more
        param5: path to save plot to, default is 'plot.png'
        param6: title of plot
        param7: name given to negative phenotype in legend
        param8: name given to positive phenotype in legend
        param9: color of negative phenotype line
        param10: color of positive phenotype line
        param11: display percent accuracy, default is true
    Saves plot to file path provided
    '''
    # reads cvs into pandas dataframe object if no dataframe provided
    if df is None:
        df = pd.read_csv(csv_file_path)

    # copies dataframe if provided to avoid mutating dataframe passed as argument
    else:
        df = df.copy()

    # select only models that have a minimum number of genes
    df =  df[df['num_genes'] > min_genes]
    
    # converts input_RMSE to percentage rank
    df['RMSE_Rank'] = df.input_RMSE.rank(pct = True)

    # selects rows based on RMSE percentage rank
    df = df[df['RMSE_Rank'] < RMSE_rank]
    # Keep only labeled binary rows (±1); drop 0/unassigned and NaNs
    tp = pd.to_numeric(df['true_phenotype'], errors='coerce')
    mask = tp.isin([-1, 1]) & df['SPS'].notna()
    df = df.loc[mask].copy()
    # Normalize phenotype labels to string "-1"/"1" for consistent pivots
    df['true_phenotype'] = tp.loc[mask].astype(int).astype(str)
    data_wide = df.pivot(columns='true_phenotype', values='SPS')

    # creates density plot using seaborn
    if '-1' in data_wide and data_wide['-1'].dropna().size > 0:
        sns.kdeplot(data_wide['-1'], color=neg_pheno_color, bw_method=bw_method, ax=axes, label=neg_pheno_name)
    if '1' in data_wide and data_wide['1'].dropna().size > 0:
        sns.kdeplot(data_wide['1'], color=pos_pheno_color, bw_method=bw_method, ax=axes, label=pos_pheno_name)
    
    # labels axes and adds title
    title = title + '\nlowest ' + "{0:.0%}".format(RMSE_rank) + ' of MFS models combined'
    
    #adds caption with percent accuracy
    if percent_accuracy:
        xlabel = 'Sequence Prediction Score\n' + "accuracy: {0:.0%}".format(calc_percent_accuracy_from_df(df))
        xlabel += "; TPR: {:.0%}; TNR: {:.0%}; balanced accuracy: {:.0%}".format(*calc_balanced_accuracy(df))
    else:
        xlabel = 'Sequence Prediction Score'
    axes.set_xlabel(xlabel)
    
    axes.set_title(title, wrap=True)
    axes.set_ylabel('Density')
    
    # style
    axes.grid(color = 'white', linestyle='-', linewidth=1)
    axes.set_facecolor('#ececec')
    sns.despine(left=True, bottom=True)

    # legend labels and title
    handles, labels = axes.get_legend_handles_labels()
    axes.legend(handles=handles, labels=labels, title='True Phenotype')
    
    # Save figure only if caller did not provide custom axes and requested a path
    if fig_path and (axes is None):
        plt.savefig(fig_path)

def create_sps_plot_violin(csv_file_path=None,
                    df=None, RMSE_rank = 0.05,
                    bw_method = 0.2,
                    fig_path = 'plot.png',
                    title = 'Model Predictions',
                    neg_pheno_name = "Negative",
                    pos_pheno_name = "Positive",
                    neg_pheno_color = '#F55E54',
                    pos_pheno_color = '#2fc8cc',
                    percent_accuracy = True,
                    axes = None,
                    min_genes = 0):
    
    # reads csv into pandas dataframe object if no dataframe provided
    if df is None:
        df = pd.read_csv(csv_file_path)

    # copies dataframe if provided to avoid mutating dataframe passed as argument
    else:
        df = df.copy()

    # select only models that have a minimum number of genes
    df =  df[df['num_genes'] > min_genes]
    
    # converts input_RMSE to percentage rank
    df['RMSE_Rank'] = df.input_RMSE.rank(pct = True)

    # selects rows based on RMSE percentage rank
    df = df[df['RMSE_Rank'] < RMSE_rank]
    # Keep only labeled binary rows (±1); drop 0/unassigned and NaNs
    tp = pd.to_numeric(df['true_phenotype'], errors='coerce')
    mask = tp.isin([-1, 1]) & df['SPS'].notna()
    df = df.loc[mask].copy()
    # Normalize phenotype labels to string "-1"/"1" for mapping
    df['true_phenotype'] = tp.loc[mask].astype(int).astype(str)

    # create the figure and subplot
    if axes is None:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8,5))
    else:
        fig = axes.get_figure()

    #adds caption with percent accuracy
    if percent_accuracy:
        xlabel = 'True Phenotype\n\n' + "accuracy: {0:.0%}".format(calc_percent_accuracy_from_df(df))
        xlabel += "; TPR: {:.0%}; TNR: {:.0%};\n balanced accuracy: {:.0%}".format(*calc_balanced_accuracy(df))
    else:
        xlabel = 'True Phenotype'

    # set SPS scores greater than 1 equal to 1 and less than -1 equal to -1
    df['SPS'] = np.clip(df['SPS'], -1, 1)

    # create a dictionary with the phenotype names
    pheno_dict = {'1': pos_pheno_name, '-1': neg_pheno_name}

    # replace the true_phenotype values in the dataframe with the corresponding names
    df['true_phenotype'] = df['true_phenotype'].map(pheno_dict)
    
    # create the violin plots for each true phenotype
    sns.violinplot(x='true_phenotype', y='SPS', data=df,
                   palette=[neg_pheno_color, pos_pheno_color],
                   split=True, inner='box', bw=bw_method, ax=axes,
                   label=[neg_pheno_name, pos_pheno_name],
                   order=[neg_pheno_name, pos_pheno_name],
                   legend=False, cut = 0,
                   linewidth = 0)

    # set the plot title and axis labels
    title_cutoff = '\n'+ '\n'.join(wrap(title, width = 33)) + '\nlowest ' + "{0:.0%}".format(RMSE_rank) + ' of MFS models combined'
    axes.set_ylabel('Sequence Prediction Score')
    axes.set_title(title_cutoff, wrap=True)

    axes.set_xlabel(xlabel)

    # style the plot
    axes.grid(color='white', linestyle='-', linewidth=.5)
    axes.set_facecolor('#ececec')
    sns.despine(left=True, bottom=True)
    
    # save the plot to a file. Let matplotlib infer format from file extension to ensure the required backend is available.
    # Save only if caller did not supply external axes and requested a path
    if fig_path and axes is None:
        fig.savefig(fig_path)


def create_continuous_plot(csv_file_path=None,
                          df=None,
                          fig_path='plot.png',
                          title='Phenotype vs SPS',
                          axes=None,
                          min_genes=0):
    """Create a 2D density plot of true phenotype vs SPS.

    Args:
        csv_file_path: Optional path to the predictions CSV file.
        df: Optional dataframe containing predictions. If provided it will
            be copied to avoid mutation.
        fig_path: Where to save the resulting figure (format inferred from
            extension).
        title: Plot title.
        axes: Optional matplotlib axes to draw on.
        min_genes: Minimum number of genes a model must have to be included.
    """
    if df is None:
        df = pd.read_csv(csv_file_path)
    else:
        df = df.copy()

    df = df[df['num_genes'] > min_genes]

    if axes is None:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 5))
    else:
        fig = axes.get_figure()

    sns.kdeplot(data=df, x='true_phenotype', y='SPS', fill=True,
                cmap='viridis', thresh=0, levels=100, ax=axes)

    axes.set_xlabel('True phenotype')
    axes.set_ylabel('SPS')
    axes.set_title(title, wrap=True)
    axes.grid(color='white', linestyle='-', linewidth=0.5)
    axes.set_facecolor('#ececec')
    sns.despine(left=True, bottom=True)

    if fig_path and axes is None:
        fig.savefig(fig_path)


# calculates percent accuracy directly from csv file- no longer use
def calc_percent_accuracy_from_csv(csv_file_path):
    csv_file = open(csv_file_path, 'r')
    num_correct = 0
    total_num = 0
    for row in csv_file.readlines()[1:]:
        row = row.strip().split(',')
        total_num += 1
        true_pheno = float(row[-1])
        sps = float(row[-2])
        if true_pheno * sps > 0:
            num_correct += 1
    csv_file.close()
    return num_correct / total_num


# calculates percent accuracy from data frame
def calc_percent_accuracy_from_df(df):
    # Normalize and filter to labeled binary rows (±1)
    tp = pd.to_numeric(df['true_phenotype'], errors='coerce')
    mask = tp.isin([-1, 1]) & df['SPS'].notna()
    sps = df.loc[mask, 'SPS']
    labels = tp.loc[mask].astype(int)
    total_num = len(sps)
    if total_num == 0:
        return 0.0
    num_correct = sum(is_row_correct(x, y) for x, y in zip(sps, labels))
    return num_correct / total_num

def calc_balanced_accuracy(df):
    # Normalize and filter to labeled binary rows (±1)
    tp = pd.to_numeric(df['true_phenotype'], errors='coerce')
    mask = tp.isin([-1, 1]) & df['SPS'].notna()
    work = df.loc[mask, ['SPS']].copy()
    labels = tp.loc[mask].astype(int)
    # Correct prediction indicator
    correct = (work['SPS'] > 0) == (labels > 0)
    # Counts for denominators
    pos_mask = labels == 1
    neg_mask = labels == -1
    pos_total = int(pos_mask.sum())
    neg_total = int(neg_mask.sum())
    TPR = (int((pos_mask & correct).sum()) / pos_total) if pos_total > 0 else 0.0
    TNR = (int((neg_mask & correct).sum()) / neg_total) if neg_total > 0 else 0.0
    balanced_acc = (TPR + TNR) / 2 if (pos_total > 0 or neg_total > 0) else 0.0
    return (TPR, TNR, balanced_acc)

def is_row_correct(SPS, true_phenotype):
    # returns 1 if correct, 0 if incorrect
    # correct sequence prediction when SPS*true_phenotype > 0
    # seq with +1 phenotype should have positive SPS, seq with -1 phenotype should have negative SPS
    if SPS*true_phenotype > 0:
        return 1
    return 0

