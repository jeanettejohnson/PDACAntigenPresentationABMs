# analyze scRNAseq from here: https://humantumoratlas.org/publications/hta12_2022_nature-genetics_daniel-cui-zhou?tab=scrna-seq
import seaborn as sns
import scanpy as sc
import anndata as ad
import pandas as pd
# Core scverse libraries
from __future__ import annotations
sc.settings.set_figure_params(dpi=50, facecolor="white")
WUSTL_HTAN_DATA = sc.read_h5ad("/Users/jeanette.johnson/HTANWUSTL/HTAN_all_adata_qc_filtered.h5ad")

sc.pp.filter_cells(WUSTL_HTAN_DATA, min_genes=200)
sc.pp.filter_genes(WUSTL_HTAN_DATA, min_cells=3)
sc.pp.filter_cells(WUSTL_HTAN_DATA, max_genes=10000)

WUSTL_HTAN_DATA.obs_names_make_unique()
sc.write("/Users/jeanette.johnson/HTANWUSTL/HTAN_all_adata_qc_filtered.h5ad", WUSTL_HTAN_DATA)


# mitochondrial genes, "MT-" for human, "Mt-" for mouse
WUSTL_HTAN_DATA.var["mt"] = WUSTL_HTAN_DATA.var_names.str.startswith("MT-")
# ribosomal genes
WUSTL_HTAN_DATA.var["ribo"] = WUSTL_HTAN_DATA.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
WUSTL_HTAN_DATA.var["hb"] = WUSTL_HTAN_DATA.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(WUSTL_HTAN_DATA, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

sc.pl.violin(
    WUSTL_HTAN_DATA,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

sc.pl.scatter(WUSTL_HTAN_DATA, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
WUSTL_HTAN_DATA = WUSTL_HTAN_DATA[WUSTL_HTAN_DATA.obs['pct_counts_mt'] < 10, :]

sns.jointplot(
    data=WUSTL_HTAN_DATA.obs,
    x="log1p_total_counts",
    y="log1p_n_genes_by_counts",
    kind="hex",
)

sns.histplot(WUSTL_HTAN_DATA.obs["pct_counts_mt"])

sc.pp.scrublet(WUSTL_HTAN_DATA, batch_key="sample")
# Saving count data
WUSTL_HTAN_DATA.layers["counts"] = WUSTL_HTAN_DATA.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(WUSTL_HTAN_DATA)
# Logarithmize the data
sc.pp.log1p(WUSTL_HTAN_DATA)

sc.pp.highly_variable_genes(WUSTL_HTAN_DATA, n_top_genes=2000, batch_key="sample")

sc.pl.highly_variable_genes(WUSTL_HTAN_DATA)
sc.tl.pca(WUSTL_HTAN_DATA)

sc.pl.pca_variance_ratio(WUSTL_HTAN_DATA, n_pcs=50, log=True)

sc.pl.pca(
    WUSTL_HTAN_DATA,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)

sc.pp.neighbors(WUSTL_HTAN_DATA, n_neighbors=30, n_pcs=30)
sc.tl.umap(WUSTL_HTAN_DATA)

sc.pl.umap(
    WUSTL_HTAN_DATA,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
)

WUSTL_HTAN_DATA.obs["treatment"]= ""
WUSTL_HTAN_DATA[WUSTL_HTAN_DATA.obs["sample"].str.startswith("HT056P1")].obs["treatment"] = "FOLFIRINOX"

samples = WUSTL_HTAN_DATA.obs["sample"]
samples = pd.DataFrame(samples)
samples['treatment'] = ""
treatments = list()

sample_names = samples["sample"].to_list()

for name in sample_names:
    if name.startswith("HT056"):
        treatments.append("FOLFIRINOX")
    elif name.startswith("HT060"):
        treatments.append("FOLFIRINOX")
    elif name.startswith("HT061"):
        treatments.append("Naive")
    elif name.startswith("HT064"):
        treatments.append("Naive")
    elif name.startswith("HT071"):
        treatments.append("Chemo-RT")
    elif name.startswith("HT085"):
        treatments.append("Naive")
    elif name.startswith("HT115"):
        treatments.append("FOLFIRINOX")
    elif name.startswith("HT121"):
        treatments.append("Naive")
    elif name.startswith("HT122"):
        treatments.append("Nab-paclitaxel")
    elif name.startswith("HT123"):
        treatments.append("Naive")    
    elif name.startswith("HT124"):
        treatments.append("Naive")   
    elif name.startswith("HT125"):
        treatments.append("Nab-paclitaxel")
    elif name.startswith("HT138"):
        treatments.append("Nab-paclitaxel")
    elif name.startswith("HT140"):
        treatments.append("Naive")  
    elif name.startswith("HT166"):
        treatments.append("FOLFIRINOX")
    elif name.startswith("HT168"):
        treatments.append("FOLFIRINOX")
    elif name.startswith("HT185"):
        treatments.append("Mixed")   
    elif name.startswith("HT190"):
        treatments.append("Nab-paclitaxel")
    elif name.startswith("HT191"):
        treatments.append("FOLFIRINOX")
    elif name.startswith("HT200"):
        treatments.append("FOLFIRINOX")
    elif name.startswith("HT204"):
        treatments.append("FOLFIRINOX")
    else:
        treatments.append("None")
        print("none!")


WUSTL_HTAN_DATA.obs["treatment"] = treatments

sc.pl.umap(
    WUSTL_HTAN_DATA,
    color="treatment",
    # Setting a smaller point size to get prevent overlap
    size=2,
)

# Using the igraph implementation and a fixed number of iterations can be significantly faster,
# especially for larger datasets
sc.tl.leiden(WUSTL_HTAN_DATA, flavor="igraph", n_iterations=2, resolution=0.5)
sc.tl.leiden(WUSTL_HTAN_DATA, key_added="leiden_res_1", n_iterations=2, resolution=1, flavor="igraph")
# sc.tl.leiden(WUSTL_HTAN_DATA)

sc.pl.umap(WUSTL_HTAN_DATA, color=["leiden"])
sc.pp.scrublet(WUSTL_HTAN_DATA, batch_key="sample")

sc.pl.umap(
    WUSTL_HTAN_DATA,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)

# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(WUSTL_HTAN_DATA, groupby="leiden", method="wilcoxon")

sc.pl.rank_genes_groups_dotplot(WUSTL_HTAN_DATA, groupby="leiden", standard_scale="var", n_genes=5)


marker_genes = {
    "Epithelial (ductal)": ["CA2", "SLC4A4", "CFTR", "KRT7", "KRT18", "KRT19", "SOX9", "HNF1B", "ONECUT1"],
    "EMT": ["KRT1", "VIM"],
    "Acinar": ["AMY2A", "CTRB2", "PRSS1"],
    "Perivascular cells ": ["ACTA2", "PDGFRB", "RGS5"],
    "Fibroblast": ["PDPN", "DCN", "LUM", "LRRC15", "SFRP2", "COL1A1", "COL1A2", "COL3A1", "PI16", "PDGFRA", "FN1", "VIM", "ENTPD1", "AIFM2", "FAP", "ACTA2"],
    "myCAF": ['TAGLN','MYL9','TPM2','MMP11','POSTN','HOPX','TWIST1','SOX4'],
    "iCAF": ['CXCL1','CXCL2','CCL2','CXCL12','PDGFRA','CFD','LMNA','DPT','HAS1','HAS2'],
    "apCAF": ['HLA-DRA','HLA-DPA1','CD74','HLA-DQA1','SLPI'],
    "Myofibroblast": ["ACTA2", "DES", "VIM"],
    # Note IGHD and IGHM are negative markers
    "B cells": [
        "MS4A1",
        "ITGB1",
        "COL4A4",
        "PRDM1",
        "IRF4",
        "PAX5",
        "BCL11A",
        "BLK",
        "IGHD",
        "IGHM",
        "CD19",
        "CD22",
        "CD40"
    ],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Lymphocytes": ["CD2","PTPRC"],
    "Plasma cell": ["SDC1", "CD27", "CD38"],
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
    "Regulatory T cell": ["IL2RA", "IL7R", "FOXP3"],
    "NK cell": ["FCGR3A", "NCAM1", "IL2RB", "CD244", "CXCR3"],
    "Dendritic cell": ["ITGAX", "S100B", "LTB", "CD1A", "CD1E", "CD1C", "PKIB", "HLA-DQA1", "CD207", "IRF8", "CLEC9A", "IDO1"],
    "Mature DC": ["LAMP3"],
    "Mast cell": ["TPSAB1"],
    "Monocyte": ["ITGAM", "VCAN", "FCN1", "SELL"],
    "Macrophage": ["CD68", "CD74", "PTPRC", "ZEB2", "HLA-DRA"],
    "Granulocytes": ["CEACAM8", "IL5RA", "FCER1A", "SIGLEC8", "ENPP3"], 
    "Tumor-infiltrating lymphocyte": ["PDCD1", "TNF", "LAG3", "TNFRSF9"],
    "CD163+ Macrophage M2 like": ["CD163", "ARG1", "TGFB1", "IL10"],
    "Inflammatory Macrophage M1 like": ["C1QC", "RGS1", "HLA-DMA", "CSF1R", "FCGR2A", "APOE", "TNF"],
    "Endocrine cells": ['GCG','LOXL4','DPP4','GC','INS','IAPP','SST','PPY']

}

sc.pl.dotplot(WUSTL_HTAN_DATA, marker_genes, groupby="leiden_res_1", standard_scale="var")


sc.pl.rank_genes_groups_dotplot(WUSTL_HTAN_DATA, groupby="leiden", standard_scale="var", n_genes=5)


sc.write("/Users/jeanette.johnson/HTANWUSTL/HTAN_all_adata_qc_filtered_with_metadata.h5ad", WUSTL_HTAN_DATA)
sc.pl.umap(
    WUSTL_HTAN_DATA,
    color=["leiden_res_1"],
    legend_loc="on data",
)