import synapseclient
syn = synapseclient.Synapse()
syn.login(authToken="eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIiwibW9kaWZ5Il0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTc2MzU4NTc3NiwiaWF0IjoxNzYzNTg1Nzc2LCJqdGkiOiIyODgxMiIsInN1YiI6IjM1NjQ5MDYifQ.FyzRhcqV1UsE0OiFp3qlrDNOYJqIG3gu9DmTFB0F5umN-slhXOYNRE-gFT__mnDYS_iYLLO6yMcTjr66xYQ_ZULVh2GDyIx0j4X3qFhBnGvu3Gk-9gOgZvQeCOQSPmzeVKje6r5yRF2Qo3Zx0YtX8tGmfJw6xnl2YbmdKn4tdDyWIvOVKpDFQm9Ae9J9B-O4eyrBobRJS54rjiRgthf9hC6KqaBawNwGEp0nt_p2jMayd5Yq4QPiFf6_Y5aX3b-xsO-tiwCs-SUXdyzyzyVLFXOnSm_Mdciby-SVqNT046YLRar14Dh0hAoDy-7S812ZqVmeGLQAXcTEvUr6RKvNMQ")
entity = syn.get("syn51117021") # matrix
import scanpy as sc

# HT056P1
# HT060P1
# HT060P1-S1R1A1G1Z1_1B2_1
syn.get("syn51117088")
syn.get("syn51117062")
syn.get("syn51116894")
HT060P1S = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT060P1S", prefix="HT060P1-S1R1A1G1Z1_1B2_1-")

# HT061P1
# HT061P1
# HT071P1
# HT071P1
# HT085P1
# HT115P1
# HT121P1

# HT122P1
# HT122P1-XB4
hta122_matrix = syn.get("syn51117047")
hta122_features = syn.get("syn51116966")
hta122_barcodes = syn.get("syn51116941")
HT122P1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT122P1-XB4", prefix="HT122P1-XB4-")

# HT123P1
# HTA12_10_8
# HT123P1-XB1 
matrix = syn.get("syn51117021") # matrix
HTA12_10_8_mat_filepath = matrix.path
features = syn.get("syn51116920") # features
HTA12_10_8_fea_filepath = features.path
barcodes = syn.get("syn51116998") # barcodes
HTA12_10_8_bar_filepath = barcodes.path
# HT123P1-XB2
ht123p1xb2_matrix = syn.get("syn51117029")
ht123p1xb2_features = syn.get("syn51116972")
ht123p1xb2_barcode = syn.get("syn51116939")
HT123P1XB2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT123P1-XB2", prefix="HT123P1-XB2-")
# HTA12_10_10
# HT123P1-XB3
hta12_10_10_matrix = syn.get("syn51117017")
hta12_10_10_features = syn.get("syn51116967")
hta12_10_10_barcodes = syn.get("syn51116896")
HT123P1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT123P1-XB3", prefix="HT123P1-XB3-")

# HT124P1
# HT124P1-XB1
hta124p1xb1_matrix = syn.get("syn51117038")
hta124p1xb1_features = syn.get("syn51116965")
hta124p1xb1_barcodes = syn.get("syn51116906")
HT124P1XB1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT124P1-XB1", prefix="HT124P1-XB1-")
# HT124P1-XB2
hta124p1xb2_matrix = syn.get("syn51117010")
hta124p1xb2_features = syn.get("syn51116981")
hta124p1xb2_barcodes = syn.get("syn51116929")
HT124P1XB2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT124P1-XB2", prefix="HT124P1-XB2-")
# HT124P1-XB3
hta124p1xb3_matrix = syn.get("syn51201375")
hta124p1xb3_features = syn.get("syn51201374")
hta124p1xb3_barcodes = syn.get("syn51201373")
HT124P1XB3 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT124P1-XB3", prefix="HT124P1-XB3-")

# HT125P1
# HT125P1-XB2
hta125p1xb2_matrix = syn.get("syn51117046")
hta125p1xb2_features = syn.get("syn51116961")
hta125p1xb2_barcodes = syn.get("syn51116917")
HT125P1XB2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT125P1-XB2", prefix="HT125P1-XB2-")
# HT125P1-XB3
hta125p1xb3_matrix = syn.get("syn51117031")
hta125p1xb3_features = syn.get("syn51116964")
hta125p1xb3_barcodes = syn.get("syn51116919")
HT125P1XB3 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT125P1-XB3", prefix="HT125P1-XB3-")
# HT125P1-XB4
ht125p1xb4_matrix = syn.get("syn51117034")
ht125p1xb4_features = syn.get("syn51116963")
ht125p1xb4_barcodes = syn.get("syn51116926")
HT125P1XB4 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT125P1-XB4", prefix="HT125P1-XB4-")
# HT125P1-XB1
ht125p1xb1_matrix = syn.get("syn51117014")
ht125p1xb1_features = syn.get("syn51116976")
ht125p1xb1_barcodes = syn.get("syn51116933")
HT125P1XB1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT125P1-XB1", prefix="HT125P1-XB1-")

# HT138P1
# HT138P1-XB1
ht138p1xb1_matrix = syn.get("syn51117020")
ht138p1xb1_features = syn.get("syn51116973")
ht138p1xb1_barcodes = syn.get("syn51116904")
HT138P1XB1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT138P1-XB1", prefix="HT138P1-XB1_1-")
# HT138P1-XB2
ht138p1xb2_matrix = syn.get("syn51117033")
ht138p1xb2_features = syn.get("syn51116986")
ht138p1xb2_barcodes = syn.get("syn51116936")
HT138P1XB2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT138P1-XB2", prefix="HT138P1-XB2_1-")
# HT138P1-XB3
ht138p1xb3_matrix = syn.get("syn51117032")
ht138p1xb3_features = syn.get("syn51116968")
ht138p1xb3_barcodes = syn.get("syn51116935")
HT138P1XB3 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT138P1-XB3", prefix="HT138P1-XB3_1-")
# HT138P1-XB4
ht138p1xb4_matrix = syn.get("syn51117042")
ht138p1xb4_features = syn.get("syn51116974")
ht138p1xb4_barcodes = syn.get("syn51116899")
HT138P1XB4 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT138P1-XB4", prefix="HT138P1-XB4_1-")

# HT140P1
# HT140P1-XB4
ht140p1xb4_matrix = syn.get("syn51117016")
ht140p1xb4_features = syn.get("syn51116978")
ht140p1xb4_barcodes = syn.get("syn51116901")
HT140P1XB4 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT140P1-XB4", prefix="HT140P1-XB4_1-")
# HT140P1-XB1
ht140p1xb1_matrix = syn.get("syn51117040")
ht140p1xb1_features = syn.get("syn51116975")
ht140p1xb1_barcodes = syn.get("syn51116913")
HT140P1XB1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT140P1-XB1", prefix="HT140P1-XB1_1-")
# HT140P1-XB2
ht140p1xb2_matrix = syn.get("syn51117049")
ht140p1xb2_features = syn.get("syn51116980")
ht140p1xb2_barcodes = syn.get("syn51116907")
HT140P1XB2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT140P1-XB2", prefix="HT140P1-XB2_1-")

# HT166P1
# HT166P1-BC1
ht166p1bc1_matrix = syn.get("syn51117043")
ht166p1bc1_features = syn.get("syn51116995")
ht166p1bc1_barcodes = syn.get("syn51116924")
HT166P1BC1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT166P1-BC1", prefix="HT166P1-BC1_1-")
# HT166P1-BC2
ht166p1bc2_matrix = syn.get("syn51117022")
ht166p1bc2_features = syn.get("syn51116983")
ht166p1bc2_barcodes = syn.get("syn51116934")
HT166P1BC2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT166P1-BC2", prefix="HT166P1-BC2_1-")
# HT166P1-BC3
ht166p1bc3_matrix = syn.get("syn51116969")
ht166p1bc3_features = syn.get("syn51117035")
ht166p1bc3_barcodes = syn.get("syn51116928")
HT166P1BC3 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT166P1-BC3", prefix="HT166P1-BC3_1-")
# HT166P1-BC4
ht166p1bc4_matrix = syn.get("syn51117058")
ht166p1bc4_features = syn.get("syn51116984")
ht166p1bc4_barcodes = syn.get("syn51116909")
HT166P1BC4 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT166P1-BC4", prefix="HT166P1-BC4_1-")

# HT168P1
# HT168P1-BC1
ht168p1bc1_matrix = syn.get("syn51117070")
ht168p1bc1_features = syn.get("syn51117008")
ht168p1bc1_barcodes = syn.get("syn51116942")
HT168P1BC1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT168P1-BC1", prefix="HT168P1-BC1_1-")
# HT168P1-BC2
ht168p1bc2_matrix = syn.get("syn51117024")
ht168p1bc2_features = syn.get("syn51116971")
ht168p1bc2_barcodes = syn.get("syn51116892")
HT168P1BC2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT168P1-BC2", prefix="HT168P1-BC2_1-")
# HT168P1-BC3
ht168p1bc3_matrix = syn.get("syn51117072")
ht168p1bc3_features = syn.get("syn51116977")
ht168p1bc3_barcodes = syn.get("syn51116944")
HT168P1BC3 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT168P1-BC3", prefix="HT168P1-BC3_1-")
# HT168P1-BC4
ht168p1bc4_matrix = syn.get("syn51117012")
ht168p1bc4_features = syn.get("syn51116960")
ht168p1bc4_barcodes = syn.get("syn51116932")
HT168P1BC4 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT168P1-BC4", prefix="HT168P1-BC4_1-")

# HT185P1
# HT185B1-S1H2A2K1G1Z1_1B1
ht185b1_matrix = syn.get("syn51117019")
ht185b1_features = syn.get("syn51116962")
ht185b1_barcodes = syn.get("syn51116911")
HT185B1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT185B1", prefix="HT185B1-S1H2A2K1G1Z1_1B1-")
# HT185P1-XBC3
ht185p1xbc3_matrix = syn.get("syn51117023")
ht185p1xbc3_features = syn.get("syn51116985")
ht185p1xbc3_barcodes = syn.get("syn51116921")
HT185P1XBC3 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT185P1XBC3", prefix="HT185P1-XBC3_1-")
# HT185P1-XBC2
ht185p1xbc2_matrix = syn.get("syn51117003")
ht185p1xbc2_features = syn.get("syn51116996")
ht185p1xbc2_barcodes = syn.get("syn51116908")
HT185P1XBC2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT185P1XBC2", prefix="HT185P1-XBC2_1-")

# HT190P1
# HT190P1-XBc4
ht190p1xbc4_matrix = syn.get("syn51117025")
ht190p1xbc4_features = syn.get("syn51116993")
ht190p1xbc4_barcodes = syn.get("syn51116945")
HT190P1XBC4 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT190P1XBC4", prefix="HT190P1-XBc4_1-")
#  HT190P1-XBc2
ht190p1xbc2_matrix = syn.get("syn51201376")
ht190p1xbc2_features = syn.get("syn51201372")
ht190p1xbc2_barcodes = syn.get("syn51201371")
HT190P1XBC2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT190P1XBC2", prefix="HT190P1-XBc2_1-")
# HT190P1-XBc3
ht190p1xbc3_matrix = syn.get("syn51117026")
ht190p1xbc3_features = syn.get("syn51116970")
ht190p1xbc3_barcodes = syn.get("syn51116914")
HT190P1XBC3 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT190P1XBC3", prefix="HT190P1-XBc3_1-")

# HT191P1
# HT191P1-XBc3
ht191p1xbc3_matrix = syn.get("syn51117002")
ht191p1xbc3_features = syn.get("syn51116997")
ht191p1xbc3_barcodes = syn.get("syn51116937")
HT191P1XBC3 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT191P1XBC3", prefix="HT191P1-XBc3_1-")
# HT191P1-XBc1
ht191p1xbc1_matrix = syn.get("syn51117039")
ht191p1xbc1_features = syn.get("syn51116979")
ht191p1xbc1_barcodes = syn.get("syn51116918")
HT191P1XBC1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT191P1XBC1", prefix="HT191P1-XBc1_1-")
# HT191P1-XBc2
ht191p1xbc2_matrix = syn.get("syn51117027")
ht191p1xbc2_features = syn.get("syn51116982")
ht191p1xbc2_barcodes = syn.get("syn51116925")
HT191P1XBC2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT191P1XBC2", prefix="HT191P1-XBc2_1-")

# HT200P1
# HT200P1-XBc1
ht200p1xbc1_matrix = syn.get("syn51117030")
ht200p1xbc1_features = syn.get("syn51116987")
ht200p1xbc1_barcodes = syn.get("syn51116927")
HT200P1XBC1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT200P1XBC1", prefix="HT200P1-XBc1_1-")
# HT200P1-XBc2
ht200p1xbc2_matrix = syn.get("syn51117028")
ht200p1xbc2_features = syn.get("syn51116991")
ht200p1xbc2_barcodes = syn.get("syn51116930")
HT200P1XBC2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT200P1XBC2", prefix="HT200P1-XBc2_1-")
# HT200P1-XBc3
ht200p1xbc3_matrix = syn.get("syn51117011")
ht200p1xbc3_features = syn.get("syn51116990")
ht200p1xbc3_barcodes = syn.get("syn51116931")
HT200P1XBC3 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT200P1XBC3", prefix="HT200P1-XBc3_1-")

# HT204P1
# HT204P1-XBc1
ht204p1xbc1_matrix = syn.get("syn51117013")
ht204p1xbc1_features = syn.get("syn51116988")
ht204p1xbc1_barcodes = syn.get("syn51116898")
HT204P1XBC1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT204P1XBC1", prefix="HT204P1-XBc1_1-")
# HT204P1-XBc2
ht204p1xbc2_matrix = syn.get("syn51117041")
ht204p1xbc2_features = syn.get("syn51116994")
ht204p1xbc2_barcodes = syn.get("syn51116897")
HT204P1XBC2 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT204P1XBC2", prefix="HT204P1-XBc2_1-")
# HT204P1-XBc3
ht204p1xbc3_matrix = syn.get("syn51116989")
ht204p1xbc3_features = syn.get("syn51116938")
ht204p1xbc3_barcodes = syn.get("syn51117018")
HT204P1XBC3 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT204P1XBC3", prefix="HT204P1-XBc3_1-")

# HT224P1
# HT224P1-S1Fc2A2N1Bmn1
syn.get("syn51116866")
syn.get("syn51116858")
syn.get("syn51116845")
HT224P1S = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT224P1S", prefix="HT224P1-S1Fc2A2N1Bmn1_1-")
# rds file: synapse get syn53214596

# HT231P1-S1H3Fc2A2N1Z1_1Bmn1
synapse get syn51116869
synapse get syn51116863
synapse get syn51116848
############################################################################
# move the barcodes, features, matrix files into a single directory

import scanpy as sc
sc.settings.set_figure_params(dpi=50, facecolor="white")

HT123P1 = sc.read_10x_mtx("/Users/jeanette.johnson/HTANWUSTL/HT123P1-XB1", prefix="HT123P1-XB1-")
HT123P1.obs_names_make_unique()
# print(HT123P1.obs["sample"].value_counts()) 
# this is where i can concat all my samples together later

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
HT123P1.var["mt"] = HT123P1.var_names.str.startswith("MT-")
# ribosomal genes
HT123P1.var["ribo"] = HT123P1.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
HT123P1.var["hb"] = HT123P1.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(HT123P1, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

sc.pl.violin(
    HT123P1,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

sc.pl.scatter(HT123P1, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

sc.pp.filter_cells(HT123P1, min_genes=100)
sc.pp.filter_genes(HT123P1, min_cells=3)

# sc.pp.scrublet(HT123P1, batch_key="sample")
sc.pp.scrublet(HT123P1)
HT123P1.layers["counts"] = HT123P1.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(HT123P1)
# Logarithmize the data
sc.pp.log1p(HT123P1)

# sc.pp.highly_variable_genes(HT123P1, n_top_genes=2000, batch_key="sample")
sc.pp.highly_variable_genes(HT123P1, n_top_genes=2000)
sc.pl.highly_variable_genes(HT123P1)

sc.tl.pca(HT123P1)
sc.pl.pca_variance_ratio(HT123P1, n_pcs=50, log=True)
sc.pl.pca(
    HT123P1,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)

sc.pp.neighbors(HT123P1)
sc.tl.umap(HT123P1)
sc.pl.umap(
    HT123P1,
    color="pct_counts_mt",
    # Setting a smaller point size to get prevent overlap
    size=2,
)

# Using the igraph implementation and a fixed number of iterations can be significantly faster,
# especially for larger datasets
sc.tl.leiden(HT123P1, n_iterations=2, flavor="igraph")
sc.pl.umap(HT123P1, color=["leiden"])

sc.pl.umap(
    HT123P1,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)
