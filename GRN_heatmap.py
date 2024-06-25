import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform


# Directory where your CSV files are stored
directory = "/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/"

# Aggregate all CSV files into a single DataFrame
# Aggregate all CSV files into a single DataFrame
all_data = pd.DataFrame()
for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        file_path = os.path.join(directory, filename)
        df = pd.read_csv(file_path, encoding='latin1')

        # Extract age and inj from the filename for additional columns
        # Ensure '.csv' is not included in the extracted parts
        parts = filename.replace('.csv', '').split('_')
        df['Age'] = parts[0]
        df['inj'] = parts[1]
        all_data = pd.concat([all_data, df], ignore_index=True)

gene = "Nfia"
# Pivot the data for sources
sources_data = all_data[all_data['source'] == gene]
pivot_sources = sources_data.pivot_table(index=['target', 'inj'], columns='Age', values='coef_mean', aggfunc='mean').reset_index()

# Pivot the data for targets
targets_data = all_data[all_data['target'] == gene]
pivot_targets = targets_data.pivot_table(index=['source', 'inj'], columns='Age', values='coef_mean', aggfunc='mean').reset_index()


############Clustering##########
def plot_clustered_heatmap(pivot_data, gene_column, inj, color_scheme='coolwarm', dpi=300):
    # Filter and prepare data for heatmap
    inj_data = pivot_data[pivot_data['inj'] == inj].drop(columns=['inj']).set_index(gene_column)
    
    # Fill NaN values and calculate pairwise distances
    inj_data_filled = inj_data.fillna(0)
    pairwise_dists = pdist(inj_data_filled, metric='euclidean')
    linkage_matrix = linkage(pairwise_dists, method='ward')

    # Order the rows according to hierarchical clustering
    dendro = dendrogram(linkage_matrix, no_plot=True)
    row_order = dendro['leaves']
    clustered_data = inj_data_filled.iloc[row_order]

    # Create the heatmap without dendrogram
    plt.figure(figsize=(4, 8))
    sns.heatmap(clustered_data, cmap=color_scheme)
    plt.title(f'Clustered Heatmap of {gene_column.capitalize()} of {gene} - Predicted.id: {inj}')
    plt.xlabel('Age')
    plt.ylabel('Gene')

    # Save with high DPI, including gene_column in the filename
    plt.savefig(f'Index_regression_GRN_predicted/{gene}/heatmap_{gene_column}_{inj}_predicted.pdf', dpi=dpi)
    plt.show()
# Iterate and plot heatmaps for sources and targets
for inj in pivot_sources['inj'].unique():
    if not pivot_sources[pivot_sources['inj'] == inj].empty:
        plot_clustered_heatmap(pivot_sources, 'target', inj, dpi=600)

for inj in pivot_targets['inj'].unique():
    if not pivot_targets[pivot_targets['inj'] == inj].empty:
        plot_clustered_heatmap(pivot_targets, 'source', inj, dpi=600)
