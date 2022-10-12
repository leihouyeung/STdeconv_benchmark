# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2021-06-16 15:19:50
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2021-08-02 23:22:10


import os
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sns.set_style("ticks")
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']

DefaulfColorPalette = [
    "#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C",
    "#DAA51B", "#2F8AC4", "#764E9F", "#ED645A", "#9EB4C7", "#BCBD22",
    "#B279A2", "#EECA3B", "#17BECF", "#FF9DA6", "#778AAE", "#1B9E77",
    "#A6761D", "#526A83", "#B82E2E", "#80B1D3", "#68855C", "#D95F02",
    "#BEBADA", "#AF6458", "#D9AF6B", "#9C9C5E", "#625377", "#8C785D",
    "#88CCEE", "#E73F74", "#FFFFB3", "#CCEBC5", "#332288", "#A65628"
]

def PlotParser(subparsers):
    workflow = subparsers.add_parser("plot", 
        help = "Visualize the deconvolution result.")

    group_input = workflow.add_argument_group("Input arguments")
    group_input.add_argument("--deconv-file", dest = "deconv_res_file", required = True,
        help = "Location of the deconvolution result file (i.e., outdir/outprefix_spot_celltype_frac.txt). ")
    group_input.add_argument("--st-loc", dest = "st_loc_file", required = True,
        help = "Location of the ST spot position file. "
        "The file should be a tab-separated plain-text file with header. "
        "The first column should be the spot name, and the second and the third column should be the row and column coordinate. ")
    group_input.add_argument("--sample-id", dest = "sample_id", default = None, type = str, 
        help = "If the count matrices are merged together to perform deconvolution, please specify the sample id of interest. ")
    group_input.add_argument("--plot-type", dest = "plot_type", default = "scatterpie",
        choices = ["scatterpie", "scatter"],
        help = "Which type of plot, scatterpie or scatter. "
        "If 'scatterpie', each spot will be shown as a pie where the proportions represent different cell-type fractions. "
        "If 'scatter', each spot will be colored according to the celltype with the maximum fraction. ")
    group_input.add_argument("--pt-size", dest = "pt_size", default = None, type = float,
        help = "The size of point in the plot. "
        "If the plot type is 'scatterpie', the size will be set as 8 by default. "
        "If the plot type is 'scatter', the size will be set as 2 by default. "
        "For reference, if the plot type is 'scatterpie', the size is recommended to set from 5 to 10. "
        "If the plot type is 'scatter', the size is recommended to set from 1 to 4. ")

    group_output = workflow.add_argument_group("Output arguments")
    group_output.add_argument("--outdir", dest = "out_dir", default = ".", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: current directory. ")
    group_output.add_argument("--outprefix", dest = "out_prefix", default = "STRIDE", 
        help = "Prefix of output files. DEFAULT: STRIDE. ")


def PieMarker(loc_list, frac_list, size, color_list):
    '''
    Return marker list for a point
    '''
    frac_cumsum = np.cumsum(frac_list)
    frac_cumsum = frac_cumsum/frac_cumsum[-1]
    marker_list = []
    previous = 0
    # calculate the points of the pie pieces
    for color, frac in zip(color_list, frac_cumsum):
        curr = frac
        x  = np.cos(2 * np.pi * np.linspace(previous, curr, 50)).tolist()
        y  = np.sin(2 * np.pi * np.linspace(previous, curr, 50)).tolist()
        xy = np.row_stack([[0, 0], np.column_stack([x, y])])
        marker_list.append({'marker': xy, 's': size*np.abs(xy).max()**2, 'facecolor': color, 'edgecolor': "darkgrey", 'linewidth': 0.1})
        previous = frac
    # scatter each of the pie pieces to create pies
    point_marker_list = []
    for marker in marker_list:
        point_marker_list.append(loc_list + [marker])
    return(point_marker_list)


def ScatterpiePlot(deconv_res_file, st_loc_file, sample_id, out_dir, out_prefix, pt_size = 8):
    '''
    Draw scatter pie plot to visualize the deconvolution result
    '''
    st_deconv_df = pd.read_csv(deconv_res_file, sep = "\t", index_col = 0, header = 0)
    st_loc_df = pd.read_csv(st_loc_file, sep = "\t", index_col = 0, header = 0)
    if sample_id:
        st_loc_df.iloc[:,2] = st_loc_df.iloc[:,2].astype(str)
        st_loc_df = st_loc_df[st_loc_df.iloc[:,2] == sample_id]
        st_deconv_df = st_deconv_df.loc[st_loc_df.index,:]
    if len(st_deconv_df.columns) > 15:
        color_pal = sns.color_palette("Spectral", len(st_deconv_df.columns))
    else:
        color_pal = DefaulfColorPalette
    fig = plt.figure()
    ax = fig.add_subplot()
    for i in st_deconv_df.index:
        deconv_list = st_deconv_df.loc[i,:]
        loc_list = st_loc_df.loc[i,:].tolist()
        point_marker_list = PieMarker(loc_list[0:2], deconv_list, pt_size**2, color_pal)
        for point_marker in point_marker_list:
            ax.scatter(point_marker[0], point_marker[1], **point_marker[-1])
    # add legends
    celltypes = st_deconv_df.columns
    patch_list = []
    for i in range(len(celltypes)):
        patch_list.append(mpatches.Patch(facecolor = color_pal[i], label = celltypes[i], edgecolor = "darkgrey", linewidth=0.1))
    ax.legend(handles = patch_list, loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 'small', frameon = False,
        handlelength=1, handleheight=1)
    ax.axis('equal')
    ax.set_xlabel(st_loc_df.columns[0])
    ax.set_ylabel(st_loc_df.columns[1])
    ax.set_title(out_prefix, pad = 15)
    # save figure
    plot_file = os.path.join(out_dir, "%s_deconv_scatterpie_plot.pdf" %(out_prefix))
    fig.savefig(plot_file, bbox_inches = "tight")
    plt.close(fig)


def ScatterPlot(deconv_res_file, st_loc_file, sample_id, out_dir, out_prefix, pt_size = 2):
    '''
    Draw scatter plot to visualize the deconvolution result
    '''
    st_deconv_df = pd.read_csv(deconv_res_file, sep = "\t", index_col = 0, header = 0)
    st_loc_df = pd.read_csv(st_loc_file, sep = "\t", index_col = 0, header = 0)
    if sample_id:
        st_loc_df.iloc[:,2] = st_loc_df.iloc[:,2].astype(str)
        st_loc_df = st_loc_df[st_loc_df.iloc[:,2] == sample_id]
        st_deconv_df = st_deconv_df.loc[st_loc_df.index,:]
    st_deconv_df['Celltype'] = st_deconv_df.idxmax(axis = 1)
    fig = plt.figure()
    ax = fig.add_subplot()
    celltypes = sorted(st_deconv_df.columns[0:-1])
    if len(st_deconv_df.columns) > 15:
        color_pal = sns.color_palette("Spectral", len(celltypes))
    else:
        color_pal = DefaulfColorPalette
    cluster_list = []
    for i, label in enumerate(celltypes):
        #add data points 
        points = plt.scatter(x = st_loc_df.loc[st_deconv_df.index[st_deconv_df['Celltype']==label], st_loc_df.columns[0]], 
                    y = st_loc_df.loc[st_deconv_df.index[st_deconv_df['Celltype']==label], st_loc_df.columns[1]], 
                    color = color_pal[i], alpha=1, s = pt_size**2)
        cluster_list.append(points)
    lgd = plt.legend(cluster_list, celltypes, bbox_to_anchor=(1, 0.5), 
               loc='center left', markerscale = 1, frameon = False, handlelength=1, handleheight=1, fontsize = 'small')
    ax.axis('equal')
    ax.set_xlabel(st_loc_df.columns[0])
    ax.set_ylabel(st_loc_df.columns[1])
    ax.set_title(out_prefix, pad = 15)
    plot_file = os.path.join(out_dir, "%s_deconv_scatter_plot.pdf" %(out_prefix))
    fig.savefig(plot_file, bbox_inches = "tight")
    plt.close(fig)


def Plot(deconv_res_file, st_loc_file, plot_type, sample_id, pt_size, out_dir, out_prefix):
    sns.set_context("notebook", font_scale=1.5)
    if not pt_size:
        if plot_type == "scatterpie":
            pt_size = 8
        else:
            pt_size = 2
    if plot_type == "scatterpie":
        ScatterpiePlot(deconv_res_file, st_loc_file, sample_id, out_dir, out_prefix, pt_size)
    if plot_type == "scatter":
        ScatterPlot(deconv_res_file, st_loc_file, sample_id, out_dir, out_prefix, pt_size)
