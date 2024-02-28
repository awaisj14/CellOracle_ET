#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 16:44:52 2023

@author: javed
"""

# 0. Import

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

adata = sc.read_h5ad("Final_celloracle_predicted.h5ad")
print(f"Cell number is :{adata.shape[0]}")
print(f"Gene number is :{adata.shape[1]}")
# Random downsampling into 30K cells if the anndata object include more than 30 K cells.
n_cells_downsample = 30000
if adata.shape[0] > n_cells_downsample:
    # Let's dowmsample into 30K cells
    sc.pp.subsample(adata, n_obs=n_cells_downsample, random_state=123)
print(f"Cell number is :{adata.shape[0]}")

df = pd.read_parquet("ATAC/base_GRN_dataframe_filtered.parquet")   

base_GRN = df

# Check data
base_GRN.head()

# Instantiate Oracle object
oracle = co.Oracle()


# Check data in anndata
print("Metadata columns :", list(adata.obs.columns))
print("Dimensional reduction: ", list(adata.obsm.keys()))


# In this notebook, we use the unscaled mRNA count for the nput of Oracle object.
adata.X = adata.layers["raw_count"].copy()

# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="Age_predicted",
                                   embedding_name="X_draw_graph_fa")
# You can load TF info dataframe with the following code.
oracle.import_TF_data(TF_info_matrix=base_GRN)

# Alternatively, if you saved the informmation as a dictionary, you can use the code below.
# oracle.import_TF_data(TFdict=TFinfo_dictionary)

#Wrote this beauty myself. I feel like GOD. Database was downloaded from: https://cdn.netbiol.org/tflink/download_files/TFLink_Mus_musculus_interactions_All_GMT_proteinName_v1.0.gmt
TFlink = pd.read_csv("/Users/Javed/Documents/Humous/TFLink_Mus_musculus_interactions_All_GMT_proteinName_v1.0.csv")
tffs = TFlink.TF
TFlink.head(n=2)
del TFlink['TF']
df = TFlink
df['Target_genes'] = df[df.columns[1:]].apply(
    lambda x: ','.join(x.dropna().astype(str)),
    axis=1
)
df = df['Target_genes']
df.head(n=2)
tffs.head(n=2)

TFsmerge = pd.merge(tffs, df, right_index=True, left_index=True)
TFsmerge.to_csv('Final_tfs_links.csv')

Paul_15_data = TFsmerge

# Make dictionary: dictionary key is TF and dictionary value is list of target genes.
TF_to_TG_dictionary = {}

for TF, TGs in zip(Paul_15_data.TF, Paul_15_data.Target_genes):
    # convert target gene to list
    TG_list = TGs.replace(" ", "").split(",")
    # store target gene list in a dictionary
    TF_to_TG_dictionary[TF] = TG_list

# We invert the dictionary above using a utility function in celloracle.
TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)


# Add TF information
oracle.addTFinfo_dictionary(TG_to_TF_dictionary)
# Perform PCA
oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)


n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)

# Save oracle object.
oracle.to_hdf5("Paul_15_data_predicted.celloracle.oracle")
# Load file.
oracle = co.load_hdf5("Paul_15_data_predicted.celloracle.oracle")


# Check clustering data
sc.pl.draw_graph(oracle.adata, color="Age_predicted")


# Calculate GRN for each population in "louvain_annot" clustering unit.
# This step may take some time.(~30 minutes)
links = oracle.get_links(cluster_name_for_GRN_unit="Age_predicted", alpha=10,
                         verbose_level=10)
links.links_dict.keys()

# Save Links object.
links.to_hdf5(file_path="links_predicted.celloracle.links")



## NETWORK PRE PROCESSING

links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

plt.rcParams["figure.figsize"] = [9, 4.5]
links.plot_degree_distributions(plot_model=True,
                                               #save=f"{save_folder}/degree_distribution/",
                                               )
plt.rcParams["figure.figsize"] = [6, 4.5]
# Calculate network scores.
links.get_network_score()
links.merged_score.head()


# Save Links object.
links.to_hdf5(file_path="links_predicted.celloracle.links")
# You can load files with the following command.
links = co.load_hdf5(file_path="links_predicted.celloracle.links")

# Check cluster name
links.cluster




cluster_name = "P04_L5 PT ALM Hpgd"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P04_Hpgd.csv")

cluster_name = "P08_L5 PT ALM Hpgd"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P08_Hpgd.csv")

cluster_name = "P12_L5 PT ALM Hpgd"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P12_Hpgd.csv")

cluster_name = "P21_L5 PT ALM Hpgd"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P21_Hpgd.csv")




cluster_name = "P04_L5 PT ALM Npsr1"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P04_Npsr1.csv")

cluster_name = "P08_L5 PT ALM Npsr1"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P08_Npsr1.csv")

cluster_name = "P12_L5 PT ALM Npsr1"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P12_Npsr1.csv")

cluster_name = "P21_L5 PT ALM Npsr1"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P21_Npsr1.csv")




cluster_name = "P04_L5 PT ALM Slco2a1"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P04_Slco2a1.csv")

cluster_name = "P08_L5 PT ALM Slco2a1"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P08_Slco2a1.csv")

cluster_name = "P12_L5 PT ALM Slco2a1"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P12_Slco2a1.csv")

cluster_name = "P21_L5 PT ALM Slco2a1"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("/Users/javed/Documents/Philipp/Index_regression_GRN_predicted/P21_Slco2a1.csv")

















