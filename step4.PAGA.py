#!/usr/bin/env python3
#-*-coding:utf-8-*-

import re
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
from multiprocessing import Pool
plt.switch_backend('agg')
matplotlib.use('Agg')

#input loom file
adata = scv.read(loom, cache=False)
adata.X = adata.X.astype('float64')
sc.tl.paga(adata, groups='cluster')
sc.pl.paga(adata, color=['sample'], save=".sample.all.pdf")
sc.pl.paga(adata, color=['cluster'], save=".cluster.all.pdf")
sc.pl.paga(adata, color=['sample'], threshold = 0.03, save=".sample.pdf")
sc.pl.paga(adata, color=['cluster'], threshold = 0.03, save=".cluster.pdf")
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['sample'], legend_loc='right margin', save='.sample.pdf')
sc.pl.draw_graph(adata, color=['cluster'], legend_loc='right margin',save='.cluster.pdf')
sc.pl.paga_compare(adata, threshold=0, title='', right_margin=0.2, size=10,edge_width_scale=0.5, legend_fontsize=12, fontsize=12,frameon=False, edges=True, save=".all.pdf")
sc.pl.paga_compare(adata, threshold=0.03, title='', right_margin=0.2, size=10,edge_width_scale=0.5, legend_fontsize=12, fontsize=12,frameon=False, edges=True, save=".pdf")

#input root_cluster
adata.uns['iroot'] = np.flatnonzero(adata.obs['cluster'] == str(root_cluster))[0]
sc.tl.dpt(adata)
sc.pl.draw_graph(adata, color=['dpt_pseudotime'], legend_loc='on data', save='.pseudotime.pdf')

#save result
def write_pd_table(df, output, sep="\t", col_names=True, row_names=True, first_colname=None, columns={}):
	if not first_colname:
		first_colname = False
		row_names = False
	if columns:
		column_names = list(columns.keys())
		df = subset(df, column=column_names)
		df=df.rename(columns=columns, inplace=False)
	df.to_csv(output, sep=sep, header=col_names, index=row_names, index_label=first_colname)

write_pd_table(adata.obs,output="Trajectory.Data.xls",first_colname="Cells",columns={'cluster':'Cluster', 'sample':'Sample','dpt_pseudotime':'Pseudotime'})
adata.write("paga.h5ad", compression='gzip')
