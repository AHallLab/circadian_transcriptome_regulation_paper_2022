import pandas as pd
import numpy as np
import tqdm

# Import triad groups and expression data
triads = pd.read_csv('grps_hc_triad_table.csv')
X_data = pd.read_csv('200721_2020_traes_circ_gene_TPM_av_exp_0_68_for_MC.csv')


# Define cross correlation function between two arrays
def cross_corr(x1, x2):
    x1 = (x1 - np.mean(x1)) / (np.std(x1) * x1.shape[0])
    x2 = (x2 - np.mean(x2)) / np.std(x2)
    corr = np.correlate(x1, x2, mode='full')  # Returns array of cross-correlation scores at different lags
    lag = corr.argmax() - (len(x1) - 1)  # The lag with the highest cross-correlation,
    # need to multiply by 4 to convert to hours
    return corr, -lag * 4


# Create results dataframe for pairs of genes in a triad
triads_uniq = list(triads['Triad'].unique())
lags = {'gene1': [], 'gene2': [], 'max_corr': [], 'max_corr_lag': [], 'triad': []}

for triad in tqdm.tqdm(triads_uniq):
    # Get genes in triad
    genes = triads.loc[triads['Triad'] == triad]
    genes = genes.iloc[0, 1:].to_list()

    # Get expression data for triad if genes are present in expression matrix
    subset = X_data.loc[X_data['transcript'].isin(genes)].T
    if subset.shape[1] < 2:
        continue
    new_header = subset.iloc[0]
    subset = subset[1:]
    subset.columns = new_header

    # If all three genes present in expression matrix
    if subset.shape[1] == 3:
        a, b = cross_corr(subset.iloc[:, 0], subset.iloc[:, 1])
        c, d = cross_corr(subset.iloc[:, 0], subset.iloc[:, 2])
        e, f = cross_corr(subset.iloc[:, 1], subset.iloc[:, 2])

        lags['gene1'].extend([subset.columns[0], subset.columns[0], subset.columns[1]])
        lags['gene2'].extend([subset.columns[1], subset.columns[2], subset.columns[2]])
        lags['max_corr'].extend([np.max(a), np.max(c), np.max(e)])
        lags['max_corr_lag'].extend([b % 24, d % 24, f % 24])
        lags['triad'].extend([triad, triad, triad])

    # If only two genes found in expression matrix
    if subset.shape[1] == 2:
        a, b = cross_corr(subset.iloc[:, 0], subset.iloc[:, 1])

        lags['gene1'].append(subset.columns[0])
        lags['gene2'].append(subset.columns[1])
        lags['max_corr'].append(np.max(a))
        lags['max_corr_lag'].append(b % 24)
        lags['triad'].append(triad)

lags_df = pd.DataFrame(lags)
