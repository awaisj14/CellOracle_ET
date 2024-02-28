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
adata = sc.read_h5ad("/Users/javed/Documents/Philipp/Final_celloracle.h5ad")

adata.obs['cluster.id.final'] = adata.obs['cluster.id3_named2'].copy()

# Merge 'distal' and 'intermediate' into 'merged_cluster'
adata.obs['cluster.id3.final'] = adata.obs['cluster.id3_named2'].replace({'distal': 'distal', 'intermediate': 'distal'})
# Step 2: Now, you can use 'cluster.id3_named2_merged' for your analyses or visualizations
# For example, to visualize the graph with the new cluster assignment
import matplotlib
matplotlib.use('Agg')  # Use the non-interactive 'Agg' backend
import matplotlib.pyplot as plt
import scanpy as sc
colors = {'proximal': '#E1BD6D', 'distal': '#0B775E'}

# Your plotting code here
sc.pl.draw_graph(adata, color='cluster.id3.final', palette=colors)

# Save the plot as a PDF file
plt.savefig('cluster.id3.final.pdf', dpi=300,bbox_inches='tight')
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
norm = mcolors.Normalize(vmin=adata.obs['Pseudotime'].min(), vmax=adata.obs['Pseudotime'].max())
mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array([])
plt.colorbar(mappable, label='Pseudotime')

plt.title('Graph layout colored by Pseudotime')
plt.xlabel('Component 1')
plt.ylabel('Component 2')
plt.show()
plt.savefig('Pseudotime.pdf', dpi=300,bbox_inches='tight')
plt.close()  # Close the plot to free up memory







### EXPRESSION

sc.pl.draw_graph(oracle.adata, color=["Meis2"],
                 layer="imputed_count", use_raw=False, cmap="viridis")
plt.savefig('Meis2_expression.pdf', dpi=300,bbox_inches='tight')
plt.close()  # Close the plot to free up memory


sc.pl.draw_graph(oracle.adata, color=["Nfia"],
                 layer="imputed_count", use_raw=False, cmap="viridis")
plt.savefig('Nfia_expression.pdf', dpi=300,bbox_inches='tight')
plt.close()  # Close the plot to free up memory


sc.pl.draw_graph(oracle.adata, color=["Zbtb16"],
                 layer="imputed_count", use_raw=False, cmap="viridis")
plt.savefig('Zbtb16_expression.pdf', dpi=300,bbox_inches='tight')
plt.close()  # Close the plot to free up memory






######## NODE NETWORK

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os
import re

# Define the target genes to highlight
highlight_genes = ["Nfia", "Meis2", "Zbtb16"]

# Define a function to process the CSV files and build the network diagrams
def process_files():
    for file in os.listdir("."):
        match = re.match(r"P(\d+)_(Hpgd|Npsr1|Slco2a1)\.csv", file)
        if match:
            age, cluster = match.groups()
            # Merge Hpgd and Npsr1 clusters into ETUpper, and treat Slco2a1 as ETlower
            cluster = "ETUpper" if cluster in ["Hpgd", "Npsr1"] else "ETlower"
            
            # Load the CSV file
            df = pd.read_csv(file)
            # Create a directed graph from the dataframe
            G = nx.from_pandas_edgelist(df, source='source', target='target', create_using=nx.DiGraph())
            
            # Position nodes using the hierarchical layout
            pos = nx.spring_layout(G, seed=42)
            
            # Draw the graph
            plt.figure(figsize=(10, 8))
            nx.draw(G, pos, with_labels=True, node_color='lightblue', 
                    node_size=500, edge_color='gray', linewidths=0.5, 
                    font_size=8, arrows=True)
            
            # Highlight specific genes
            nx.draw_networkx_nodes(G, pos, nodelist=highlight_genes, node_color='red')
            nx.draw_networkx_labels(G, pos, labels={gene: gene for gene in highlight_genes}, font_color='white')
            
            # Save the plot to a PDF file
            plt.title(f"Gene Network for Age {age}, Cluster {cluster}")
            plt.savefig(f"{age}_{cluster}_network.pdf")
            plt.close()

# Run the function to process the files and generate the diagrams
process_files()





######### Sankey
import pandas as pd
import glob
import os
from matplotlib.sankey import Sankey

# Detect CSV files
files = glob.glob('*.csv')

# Define a function to preprocess and merge CSV files
def preprocess_files(files):
    all_data = []
    for file in files:
        # Extract age and cluster from filename
        age, cluster = os.path.splitext(file)[0].split('_')
        if cluster in ['Hpgd', 'Npsr1']:
            cluster = 'ETUpper'
        elif cluster == 'Slco2a1':
            cluster = 'ETlower'
        
        # Read CSV file
        data = pd.read_csv(file)
        data['Age'] = age
        data['Cluster'] = cluster
        all_data.append(data)
    
    # Combine all data into a single DataFrame
    combined_data = pd.concat(all_data, ignore_index=True)
    
    # Filter for specific genes and their targets
    filtered_data = combined_data[(combined_data['source'].isin(['Nfia', 'Meis2', 'Zbtb16'])) | 
                                  (combined_data['target'].isin(['Nfia', 'Meis2', 'Zbtb16']))]
    
    return filtered_data

# Process the CSV files
processed_data = preprocess_files(files)

# This is a simplified version to demonstrate the concept. Adjustments may be needed for your specific case.

def visualize_network(data):
    # Example of building a simple Sankey diagram - you would need to expand this
    # to handle the hierarchical aspect and multiple levels.
    Sankey(flows=[1, -1], labels=['Source', 'Target']).finish()
    # Save the figure as a PDF
    plt.savefig('network_diagram.pdf')

# Since building the actual hierarchical network requires detailed data analysis
# and could be quite complex, the above code is a placeholder to show the direction.

# You might need to develop a custom function to map out the network hierarchy
# and then adjust the Sankey diagram accordingly.
