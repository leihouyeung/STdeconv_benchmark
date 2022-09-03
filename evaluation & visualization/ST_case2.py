import numpy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

if __name__ == "__main__":

    # Plot for markers
    st_data     = "PDAC_stdata.csv"
    loc_data    = "st_location_PDAC.csv"

    df_st = pd.read_csv(st_data, header=0, index_col=0)
    df_st = pd.DataFrame.transpose(df_st)
    # df_st should be Spots x Genes

    df_loc = pd.read_csv(loc_data, header=0, index_col=0)
    # df loc should be Spots x (x, y)

    X = df_loc['x'].to_numpy()
    Y = df_loc['y'].to_numpy()
    # Coordinates to numpy array

    gene_name = ['GAPDH']
    gene_value = numpy.zeros(X.shape)
    for name in gene_name:
        gene_value += df_st[name].to_numpy()

    # Ploting in a subplot
    figure, axis = plt.subplots(4, 5, figsize=(8, 8))
    figure.tight_layout()
    for i in range(4):
        for j in range(5):
            axis[i, j].set_xticks([])
            axis[i, j].set_yticks([])
            axis[i, j].spines['top'].set_visible(False)
            axis[i, j].spines['right'].set_visible(False)
            axis[i, j].spines['bottom'].set_visible(False)
            axis[i, j].spines['left'].set_visible(False)

    # Plot for marker
    axis[0, 0].set_title("GAPDH", size=7)
    im_marker = axis[0, 0].scatter(X, Y, c=gene_value, cmap="Blues", s = 1)
    cbar = figure.colorbar(im_marker, ax=axis[0, 0], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Checkpoint
    # plt.savefig("test.tiff", dpi=600)
    # exit()


    # Plot for Celltypes

    # Berglund
    pred_data = "Berglund.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    cell_type_name = "Cancer.clone.A"
    # cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[0, 1].set_title("Berglund\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)), size=7)
    im_berglund = axis[0, 1].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[0, 1], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value)**2) * np.sum((cell_type_prop - mean_cell_prop)**2))
    print("Berglund Corr: {}".format(up/down))

    # Checkpoint
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # CARD
    pred_data = "CARD.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    # cell_type_name = "Cancer.clone.A"
    cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[0, 2].set_title("CARD\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)), size=7)
    im_berglund = axis[0, 2].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[0, 2], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value)**2) * np.sum((cell_type_prop - mean_cell_prop)**2))
    print("CARD Corr: {}".format(up/down))

    # Checkpoint
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # Cell2location
    pred_data = "Cell2location.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    cell_type_name = "Cancer.clone.A"
    # cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[0, 3].set_title("Cell2location\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)), size=7)
    im_berglund = axis[0, 3].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[0, 3], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("Cell2location Corr: {}".format(up / down))

    # Checkpoint
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # DestVI
    pred_data = "DestVI.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    # cell_type_name = "Cancer.clone.A"
    cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[0, 4].set_title("DestVI\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[0, 4].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[0, 4], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("DestVI Corr: {}".format(up / down))

    # Checkpoint
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # DTSG
    pred_data = "DSTG.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    cell_type_name = "Cancer.clone.A"
    # cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[1, 0].set_title("DSTG\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[1, 0].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[1, 0], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("DSTG Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
    # plt.savefig("test.tiff", dpi=600)
    # exit()


    # NMFreg
    pred_data = "NMFreg.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    # cell_type_name = "Cancer.clone.A"
    cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[1, 1].set_title("NMFreg\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[1, 1].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[1, 1], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("NMFreg Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # RCTD
    pred_data = "RCTD.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    # cell_type_name = "Cancer.clone.A"
    cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[1, 2].set_title("RCTD\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[1, 2].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[1, 2], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("RCTD Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # SD2
    pred_data = "SD2.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    cell_type_name = "Cancer.clone.A"
    # cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[1, 3].set_title("SD2\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[1, 3].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[1, 3], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("SD2 Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # SpatialDecon
    pred_data = "SpatialDecon.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    # cell_type_name = "Cancer.clone.A"
    cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[1, 4].set_title("SpatialDecon\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[1, 4].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[1, 4], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("SpatialDecon Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # SpatialDWLS
    pred_data = "SpatialDWLS.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    # cell_type_name = "Cancer.clone.A"
    cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[2, 0].set_title("SpatialDWLS\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[2, 0].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[2, 0], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("SpatialDWLS Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # SpiceMix
    pred_data = "SpiceMix.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    cell_type_name = "Cancer.clone.A"
    # cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[2, 1].set_title("SpiceMix\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[2, 1].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[2, 1], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("SpiceMix Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # SPOTlight
    pred_data = "SPOTlight.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    cell_type_name = "Cancer.clone.A"
    # cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[2, 2].set_title("SPOTlight\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[2, 2].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[2, 2], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("SPOTlight Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # STdeconvolve
    pred_data = "STdeconvolve.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    cell_type_name = "Cancer.clone.A"
    # cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[2, 3].set_title("STdeconvolve\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[2, 3].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[2, 3], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("STdeconvolve Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    # plt.savefig("test.tiff", dpi=600)
    # exit()


    # stereoscope
    pred_data = "stereoscope.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    cell_type_name = "Cancer.clone.A"
    # cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[2, 4].set_title("stereoscope\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[2, 4].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[2, 4], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("stereoscope Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # STRIDE
    pred_data = "STRIDE.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    cell_type_name = "Cancer.clone.A"
    # cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[3, 0].set_title("STRIDE\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[3, 0].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[3, 0], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("STRIDE Corr: {}".format(up / down))

    # Checkpoint
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    # plt.savefig("test.tiff", dpi=600)
    # exit()

    # Tangram
    pred_data = "Tangram.csv"
    df_pred = pd.read_csv(pred_data, header=0, index_col=0)
    cell_type_name = "Cancer.clone.A"
    # cell_type_name = "Cancer clone A"
    cell_type_prop = df_pred[cell_type_name].to_numpy()
    # Plot
    axis[3, 1].set_title("Tangram\nCorr: {}".format(round(stats.pearsonr(cell_type_prop, gene_value)[0], 2)),
                         size=7)
    im_berglund = axis[3, 1].scatter(X, Y, c=cell_type_prop, cmap="Blues", s=1)
    cbar = figure.colorbar(im_berglund, ax=axis[3, 1], location="bottom", orientation="horizontal", shrink=1.0)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

    # Printing out info
    mean_gene_value = np.mean(gene_value)
    mean_cell_prop = np.mean(cell_type_prop)
    up = np.sum((gene_value - mean_gene_value) * (cell_type_prop - mean_cell_prop))
    down = np.sqrt(np.sum((gene_value - mean_gene_value) ** 2) * np.sum((cell_type_prop - mean_cell_prop) ** 2))
    print("Tangram Corr: {}".format(up / down))

    # Checkpoint
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
    plt.savefig("ST_case2.png", dpi=300)
    exit()

