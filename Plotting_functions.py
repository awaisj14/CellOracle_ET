#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 13:18:30 2024

@author: javed
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import celloracle as co
co.__version__
'0.10.14'
# visualization settings
%config InlineBackend.figure_format = 'retina'
%matplotlib inline

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300
save_folder = "figures"
os.makedirs(save_folder, exist_ok=True)


###### CLUSTERS
adata = sc.read_h5ad("/Users/javed/Documents/Philipp/Final_celloracle_predicted_v2.h5ad")

import matplotlib
matplotlib.use('Agg')  # Use the non-interactive 'Agg' backend
import matplotlib.pyplot as plt
import scanpy as sc
colors = {'proximal': '#E1BD6D', 'distal': '#0B775E'}

# Your plotting code here
sc.pl.draw_graph(adata, color='new_clusters', palette=colors)

# Save the plot as a PDF file
plt.savefig('new_clusters.pdf', dpi=600,bbox_inches='tight')
plt.close()  # Close the plot to free up memory

###### AGE
import matplotlib
matplotlib.use('Agg')  # Use the non-interactive 'Agg' backend
import matplotlib.pyplot as plt
import scanpy as sc
colors = {'P04': '#E3E6BB', 'P08': '#92AAA8','P12': '#6087C5','P21': '#1A4070'}
adata.obs["Age"]
# Your plotting code here
sc.pl.draw_graph(adata, color='Age', palette=colors)

# Save the plot as a PDF file
plt.savefig('nAge.pdf', dpi=600,bbox_inches='tight')
plt.close()  # Close the plot to free up memory



##### PSUEDOTIME



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scanpy as sc


# Assuming 'adata' is already computed with sc.tl.draw_graph

# Normalize 'Pseudotime' to range between 0 and 1
adata.obs['Pseudotime_normalized'] = (adata.obs['Pseudotime'] - adata.obs['Pseudotime'].min()) / (adata.obs['Pseudotime'].max() - adata.obs['Pseudotime'].min())

# Create a custom color map from #D2D4AB to #22518A
cmap = mcolors.LinearSegmentedColormap.from_list("pseudotime_gradient", ["#D2D4AB", "#22518A"])

# Get the graph layout coordinates
coords = adata.obsm['X_draw_graph_fa']  # or the key corresponding to your layout

# Apply the colormap to the normalized 'Pseudotime' values to get colors for each point
point_colors = cmap(adata.obs['Pseudotime_normalized'])

# Plotting
plt.figure(figsize=(10, 6))
plt.scatter(coords[:, 0], coords[:, 1], c=point_colors, s=20, edgecolor='none')

# Optionally, create and add a colorbar
plt.title('Graph layout colored by Pseudotime')
plt.xlabel('Component 1')
plt.ylabel('Component 2')
plt.show()
plt.savefig('Pseudotime.pdf', dpi=600,bbox_inches='tight')
plt.close()  # Close the plot to free up memory







### EXPRESSION

sc.pl.draw_graph(oracle.adata, color=["Meis2"],
                 layer="imputed_count", use_raw=False, cmap="viridis")
plt.savefig('Meis2_expression.pdf', dpi=600,bbox_inches='tight')
plt.close()  # Close the plot to free up memory


sc.pl.draw_graph(oracle.adata, color=["Nfia"],
                 layer="imputed_count", use_raw=False, cmap="viridis")
plt.savefig('Nfia_expression.pdf', dpi=600,bbox_inches='tight')
plt.close()  # Close the plot to free up memory


sc.pl.draw_graph(oracle.adata, color=["Zbtb16"],
                 layer="imputed_count", use_raw=False, cmap="viridis")
plt.savefig('Zbtb16_expression.pdf', dpi=600,bbox_inches='tight')
plt.close()  # Close the plot to free up memory










import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import os

def read_data(file_path):
    return pd.read_csv(file_path)

# Now TFF is a list of transcription factors
TFF = ['Zbtb16', 'Meis2', 'Nfia']  

def create_network(df, tfs):
    G = nx.DiGraph()

    # Filter the DataFrame to include only rows where any of the TFFs is a source or target
    df_filtered = df[df['source'].isin(tfs) | df['target'].isin(tfs)]

    # Adding nodes and edges for filtered data
    for index, row in df_filtered.iterrows():
        source = row['source']
        target = row['target']
        coef_mean = row['coef_mean']
        logp = row['-logp']

        # Add or update nodes
        if source not in G:
            G.add_node(source, coef_mean=coef_mean)
        if target not in G:
            G.add_node(target, coef_mean=coef_mean)

        # Add edges with attributes
        G.add_edge(source, target, weight=logp)

    return G
def draw_network(G, tfs, file_name, age, inj):
    plt.figure(figsize=(14, 12))
    pos = nx.spring_layout(G, scale=2)  # Initial layout for positions
    coef_mean_values = np.array([G.nodes[node]['coef_mean'] for node in G.nodes()])
    logp_values = np.array([G[u][v]['weight'] for u, v in G.edges()])
    # Identify upstream genes (nodes that have edges pointing to TFFs)
    upstream_nodes = set()
    for tf in tfs:
        upstream_nodes.update([source for source, target in G.edges() if target == tf])
    upstream_nodes.difference_update(tfs)  # Remove TFFs if they are self-regulating

    # Position the upstream genes at the very top
    upstream_x_positions = np.linspace(-1, 1, len(upstream_nodes))
    for node, x_pos in zip(sorted(upstream_nodes), upstream_x_positions):
        pos[node] = np.array([x_pos, 1])  # Set y=1 for upstream genes, positioned at the top

    # Position the TFFs just below the upstream genes
    tff_x_positions = np.linspace(-0.5, 0.5, len(tfs))
    for tf, x_pos in zip(tfs, tff_x_positions):
        pos[tf] = np.array([x_pos, 0.8])  # Set y=0.8 for TFFs

    # Identify downstream genes (non-TFF and non-upstream nodes)
    non_tff_nodes = [node for node in G.nodes() if node not in tfs and node not in upstream_nodes]

    # Evenly distribute downstream genes below the TFFs
    if non_tff_nodes:
        grid_width = np.ceil(np.sqrt(len(non_tff_nodes)))
        grid_height = np.ceil(len(non_tff_nodes) / grid_width)
        x_positions = np.linspace(-1, 1, int(grid_width))
        y_positions = np.linspace(-0.9, 0.6, int(grid_height))  # Adjust y range for downstream nodes
        positions = [(x, y) for y in y_positions for x in x_positions]
        for node, (x, y) in zip(non_tff_nodes, positions):
            pos[node] = np.array([x, y])

    # Define node colors and sizes
    node_colors = ['lightblue' if G.nodes[node]['coef_mean'] >= 0 else 'red' for node in G.nodes()]
    node_sizes = [abs(G.nodes[node]['coef_mean']) * 3000 for node in G.nodes()]

    # Drawing the network
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes)
    nx.draw_networkx_edges(G, pos, width=[G[u][v]['weight'] * 0.2 for u, v in G.edges()], edge_color='lightgrey', arrowsize=10, arrowstyle='-|>')
    nx.draw_networkx_labels(G, pos, font_color='black')
# Select representative values for node sizes and edge widths for the legend
    representative_node_sizes = np.percentile(coef_mean_values, [25, 50, 75])
    representative_edge_widths = np.percentile(logp_values, [25, 50, 75])

    # Create legends with actual values and adjusted sizes
    node_legend_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
                                  markersize=np.sqrt(abs(size)) * 40,  # Adjust multiplier for visual preference
                                  label=f'{size:.2f}') for size in representative_node_sizes]
    edge_legend_handles = [Line2D([0], [0], color='black', lw=width / max(representative_edge_widths) * 2.5,  # Normalize and scale for visibility
                                   label=f'{width:.2f}') for width in representative_edge_widths]

    # Node color legend
    color_legend_handles = [Patch(facecolor='lightblue', edgecolor='b', label='Positive coef_mean'),
                            Patch(facecolor='red', edgecolor='r', label='Negative coef_mean')]

    # Combine legend handles
    legend_handles = node_legend_handles + edge_legend_handles + color_legend_handles
    plt.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.5), title='Legend')

    plt.title(f"GRN of {', '.join(tfs)} at {age} in {inj} identity neurons")
    plt.axis('off')
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust to accommodate legend
    plt.savefig(file_name)
    plt.show()

# Directory where your CSV files are stored
directory = "/users/javed/Documents/Philipp/Index_regression_GRN_predicted/CellOracleShifts_pdfs_v2/GRNs/"

# Iterate over each file in the directory
for filename in os.listdir(directory):
    if filename.endswith(".csv"):  # Make sure to process only CSV files
        parts = filename.split('_')
        age, inj = parts[0], parts[1]

        file_path = os.path.join(directory, filename)
        df = read_data(file_path)
        G = create_network(df, TFF)

        output_file = f"network_plot_{age}_{inj}_{'_'.join(TFF)}.pdf"
        draw_network(G, TFF, output_file, age, inj)
        print(f"Generated plot for {filename} saved as {output_file}")
