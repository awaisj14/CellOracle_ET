#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 16:24:35 2023

@author: javed
"""

# Install velocyto from conda using this code if you have M1max MacOS
#conda install -c bioconda velocyto.py

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scvelo
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=600, facecolor='white')

# visualization settings
%config InlineBackend.figure_format = 'retina'
%matplotlib inline

plt.rcParams['figure.figsize'] = [6, 3]
plt.rcParams["savefig.dpi"] = 600

Phil = sc.read_h5ad("Phil.h5ad")

adata = Phil

adata = adata[adata.obs['Age'].isin(['P04', 'P08','P12','P21'])] 
adata.obs['Age_inj'] = adata.obs.apply(lambda x: f"{x['Age']}_{x['injSite']}", axis=1)

# Only consider genes with more than 1 count
sc.pp.filter_genes(adata, min_counts=1)


# Normalize gene expression matrix with total UMI count per cell
sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')


# Select top 2000 highly-variable genes
filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                              flavor='cell_ranger',
                                              n_top_genes=3000,
                                              log=False)

# Subset the genes
adata = adata[:, filter_result.gene_subset]

# Renormalize after filtering
sc.pp.normalize_per_cell(adata)


# keep raw cont data before log transformation
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()

# Log transformation and scaling
sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.pp.regress_out(adata, ['indexing'])

# PCA
sc.tl.pca(adata, svd_solver='arpack')

# Diffusion map
sc.pp.neighbors(adata)

sc.tl.diffmap(adata)
# Calculate neihbors again based on diffusionmap
sc.tl.umap(adata)

sc.pl.umap(adata, color="injSite")
sc.pl.umap(adata, color="indexing")
sc.pl.umap(adata, color="Age")
adata.obs['Age_predicted'] = adata.obs.apply(lambda x: f"{x['Age']}_{x['predicted.id']}", axis=1)

# PAGA graph construction
sc.tl.paga(adata, groups='Age_predicted')
plt.rcParams["figure.figsize"] = [6, 4.5]
sc.pl.paga(adata)

sc.tl.draw_graph(adata, init_pos='paga', random_state=123)
sc.pl.draw_graph(adata, color='Age_predicted')

adata.write_h5ad("Phil_only_CO_v3.h5ad")
adata = sc.read_h5ad("Phil_only_CO_v3.h5ad")

adata.obs['Age_predicted'] = adata.obs.apply(lambda x: f"{x['Age']}_{x['predicted.id']}", axis=1)

adata.write_h5ad("Phil_only_CO_v4.h5ad")

import copy
import glob
import time
import os
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from tqdm.auto import tqdm



import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm

#install bedtools2 in your conda environment
#git clone https://github.com/arq5x/bedtools2.git
#cd bedtools2
#make clean; make all
#bin/bedtools --version
#bedtools v2.20.1-4-gb877b35


#import time
import velocyto
import celloracle as co
from celloracle.applications import Pseudotime_calculator
co.__version__

from celloracle import motif_analysis as ma
import celloracle as co
co.__version__

from pybedtools import BedTool



#plt.rcParams["font.family"] = "arial"
plt.rcParams["figure.figsize"] = [10,5]
%config InlineBackend.figure_format = 'retina'
plt.rcParams["savefig.dpi"] = 600

%matplotlib inline



pt = Pseudotime_calculator(adata=adata,
                           obsm_key="X_draw_graph_fa", # Dimensional reduction data name
                           cluster_column_name="Age_predicted" # Clustering data name
                           )


print("Clustering name: ", pt.cluster_column_name)
print("Cluster list", pt.cluster_list)
# Check data
pt.plot_cluster(fontsize=8)



# Here, clusters can be classified into either MEP lineage or GMP lineage

clusters_in_Hpgd_lineage = ['P04_L5 PT ALM Hpgd','P08_L5 PT ALM Hpgd','P12_L5 PT ALM Hpgd','P21_L5 PT ALM Hpgd']
clusters_in_Npsr1_lineage = ['P04_L5 PT ALM Npsr1','P08_L5 PT ALM Npsr1','P12_L5 PT ALM Npsr1','P21_L5 PT ALM Npsr1']
clusters_in_Slco2a1_lineage = ['P04_L5 PT ALM Slco2a1','P08_L5 PT ALM Slco2a1','P12_L5 PT ALM Slco2a1','P21_L5 PT ALM Slco2a1']

# Make a dictionary
lineage_dictionary = {"Hpgd": clusters_in_Hpgd_lineage,"Npsr1": clusters_in_Npsr1_lineage,"Slco2a1": clusters_in_Slco2a1_lineage}

# Input lineage information into pseudotime object
pt.set_lineage(lineage_dictionary=lineage_dictionary)

# Visualize lineage information
pt.plot_lineages()


import plotly.express as px
import plotly.io as pio
pio.renderers.default = "browser"
try:
    import plotly.express as px
    def plot(adata, embedding_key, cluster_column_name):
        embedding = adata.obsm[embedding_key]
        df = pd.DataFrame(embedding, columns=["x", "y"])
        df["cluster"] = adata.obs[cluster_column_name].values
        df["label"] = adata.obs.index.values
        fig = px.scatter(df, x="x", y="y", hover_name=df["label"], color="cluster")
        fig.show()

    plot(adata=pt.adata,
         embedding_key=pt.obsm_key,
         cluster_column_name=pt.cluster_column_name)
except:
    print("Plotly not found in your environment. Did you install plotly? Please read the instruction above.")

# Estimated root cell name for each lineage
root_cells = {"Slco2a1": "P.MO4.12_AACCACAGTTAAACAG-1","Hpgd": "P.MO4.12_GCCATTCAGGTGCAGT-1","Npsr1": "L5ET_4812_TGAGTCACATCCGATA-1"}
pt.set_root_cells(root_cells=root_cells)


pt.plot_root_cells()

"umap" in pt.adata.obsm

pt.get_pseudotime_per_each_lineage()

# Check results
pt.plot_pseudotime(cmap="rainbow")


pt.adata.obs[["Pseudotime"]].head()


adata.obs = pt.adata.obs

# Save updated anndata object
adata.write_h5ad("Final_celloracle_predicted.h5ad")

adata = sc.read_h5ad("Final_celloracle_predicted.h5ad")



print(f"Cell number is :{adata.shape[0]}")
print(f"Gene number is :{adata.shape[1]}")

# Random downsampling into 30K cells if the anndata object include more than 30 K cells.
n_cells_downsample = 30000
if adata.shape[0] > n_cells_downsample:
    # Let's dowmsample into 30K cells
    sc.pp.subsample(adata, n_obs=n_cells_downsample, random_state=123)












#ALL the below part just needs peaks


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline

import seaborn as sns

import os, sys, shutil, importlib, glob
from tqdm import tqdm_notebook as tqdm

%config InlineBackend.figure_format = 'retina'

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

# Import celloracle function
from celloracle import motif_analysis as ma



# Load bed_file
file_path_of_bed_file = "L5_ET.sort.bed"
bed = ma.read_bed(file_path_of_bed_file)
print(bed.shape)
bed.head()

# Convert bed file into peak name list
peaks = ma.process_bed_file.df_to_list_peakstr(bed)
peaks

tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="mm10")

# Check results
tss_annotated.tail()


# Change format
peak_id_tss = ma.process_bed_file.df_to_list_peakstr(tss_annotated)
tss_annotated = pd.DataFrame({"peak_id": peak_id_tss,
                              "gene_short_name": tss_annotated.gene_short_name.values})
tss_annotated = tss_annotated.reset_index(drop=True)
print(tss_annotated.shape)
tss_annotated.head()


tss_annotated.to_csv("processed_peak_file.csv")


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import seaborn as sns

import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object
co.__version__
'0.13.1'
%config InlineBackend.figure_format = 'retina'
%matplotlib inline

plt.rcParams['figure.figsize'] = (15,7)
plt.rcParams["savefig.dpi"] = 600


ref_genome = "mm10"

genome_installation = ma.is_genome_installed(ref_genome=ref_genome)
print(ref_genome, "installation: ", genome_installation)
if not genome_installation:
    import genomepy
    genomepy.install_genome(name=ref_genome, provider="UCSC")
else:
    print(ref_genome, "is installed.")
    

    
# Load annotated peak data.
peaks = pd.read_csv("processed_peak_file.csv", index_col=0)
peaks.head()

def decompose_chrstr(peak_str):
    """
    Args:
        peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'

    Returns:
        tuple: chromosome name, start position, end position
    """

    *chr_, start, end = peak_str.split("_")
    chr_ = "_".join(chr_)
    return chr_, start, end

from genomepy import Genome

def check_peak_format(peaks_df, ref_genome):
    """
    Check peak format.
     (1) Check chromosome name.
     (2) Check peak size (length) and remove sort DNA sequences (<5bp)

    """

    df = peaks_df.copy()

    n_peaks_before = df.shape[0]

    # Decompose peaks and make df
    decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
    df_decomposed = pd.DataFrame(np.array(decomposed), index=peaks_df.index)
    df_decomposed.columns = ["chr", "start", "end"]
    df_decomposed["start"] = df_decomposed["start"].astype(int)
    df_decomposed["end"] = df_decomposed["end"].astype(int)

    # Load genome data
    genome_data = Genome(ref_genome)
    all_chr_list = list(genome_data.keys())


    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])


    # Filter peaks with invalid chromosome name
    n_threshold = 5
    df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]

    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])

    # Data counting
    n_invalid_length = len(lengths[lengths < n_threshold])
    n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
    n_peaks_after = df.shape[0]


    #
    print("Peaks before filtering: ", n_peaks_before)
    print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
    print("Peaks with invalid length: ", n_invalid_length)
    print("Peaks after filtering: ", n_peaks_after)

    return df


peaks = check_peak_format(peaks, ref_genome)

# Instantiate TFinfo object
tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=ref_genome)

# Scan motifs. !!CAUTION!! This step may take several hours if you have many peaks!
tfi.scan(fpr=0.02,
         motifs=None,  # If you enter None, default motifs will be loaded.
         verbose=True)

# Save tfinfo object
tfi.to_hdf5(file_path="test1.celloracle.tfinfo")
# Check motif scan results
tfi.scanned_df.head()


# Reset filtering
tfi.reset_filtering()

# Do filtering
tfi.filter_motifs_by_score(threshold=10)

# Format post-filtering results.
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)


df = tfi.to_dataframe()
df.head()


# Save result as a dataframe
df = tfi.to_dataframe()
df.to_parquet("base_GRN_dataframe_filtered.parquet")

df = pd.read_parquet("ATAC/base_GRN_dataframe_filtered.parquet")   
