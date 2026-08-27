"""
analyze_imc_features.py
-----------------------
Sanity check on a trained TissueMosaic checkpoint: embed patches from every IMC
ROI, then look at whether the embedding organises by anything biological.

    python tissuemosaic/analyze_imc_features.py --ckpt tissuemosaic/runs/dino_imc/ckpt_last.pt

Produces, in --out_dir:
    patch_umap.png        UMAP of patch embeddings, panelled by annotation
    patch_features.npz    embeddings + composition + Moran + ROI/patient labels

What to look for: a well-trained model gives a UMAP whose structure tracks
cell-type composition -- tumour-rich, stroma-rich and T-cell-rich patches
separating -- rather than tracking patient identity, which would mean the model
has latched onto batch. After a short run expect neither; the point of running
it then is only to prove the featurisation path works.

Shared helpers (paths, model loading, palette, to_numpy) come from imc_tm.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import imc_tm  # applies tm_compat on import

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import torch  # noqa: E402
from anndata import read_h5ad  # noqa: E402
from tissuemosaic.data.dataset import CropperSparseTensor  # noqa: E402
from tissuemosaic.models.patch_analyzer import Composition, SpatialAutocorrelation  # noqa: E402
from tissuemosaic.utils import SmartPca, SmartUmap  # noqa: E402


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ckpt", required=True, help="checkpoint written by run_train_imc.py")
    p.add_argument("--anndata_dir", default=str(imc_tm.H5AD_DIR))
    p.add_argument("--out_dir", default=str(imc_tm.REPO / "analysis/tissuemosaic"))
    p.add_argument("--n_patches", type=int, default=200, help="patches sampled per ROI")
    p.add_argument("--batch_size", type=int, default=64)
    p.add_argument("--feature_key", default="dino")
    p.add_argument("--umap_neighbors", type=int, default=25)
    return p.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(imc_tm.compat_summary())

    model = imc_tm.load_model(args.ckpt)
    from tissuemosaic.data import AnndataFolderDM
    dm = AnndataFolderDM(**model._hparams)
    channels = sorted(dm._categories_to_channels, key=dm._categories_to_channels.get)
    print("{} channels: {}".format(len(channels), ", ".join(channels)))

    use_gpu = torch.cuda.is_available()
    if use_gpu:
        model = model.cuda()

    files = sorted(f for f in Path(args.anndata_dir).iterdir() if f.suffix == ".h5ad")
    print("featurizing {} ROIs, up to {} patches each".format(len(files), args.n_patches))

    feats, comps, morans, roi_labels = [], [], [], []
    for i, path in enumerate(files):
        roi = path.stem
        sp_img = dm.anndata_to_sparseimage(read_h5ad(path))
        if use_gpu:
            sp_img = sp_img.cuda()

        sp_img.compute_patch_features(
            feature_name=args.feature_key, datamodule=dm, model=model,
            batch_size=args.batch_size, strategy="random", remove_overlap=False,
            n_patches_max=args.n_patches, overwrite=True)

        values, xywh = sp_img.read_from_patch_dictionary(key=args.feature_key)

        # composition and spatial autocorrelation of the SAME patches, for annotation
        tensors, _, _ = CropperSparseTensor.reapply_crops(sp_img.data, xywh)
        comp = torch.stack(Composition(return_fraction=True)(tensors), dim=0).cpu().numpy()
        moran = torch.stack(
            SpatialAutocorrelation(modality="moran", n_neighbours=6, neigh_correct=True)(tensors),
            dim=0).cpu().numpy()

        feats.append(imc_tm.to_numpy(values))
        comps.append(comp)
        morans.append(moran.max(axis=-1))
        roi_labels += [roi] * len(values)
        print("  [{:2d}/{}] {:<14s} {:4d} patches".format(i + 1, len(files), roi, len(values)))
        del sp_img

    feats = np.concatenate(feats, axis=0)
    comps = np.concatenate(comps, axis=0)
    morans = np.concatenate(morans, axis=0)
    roi_labels = np.asarray(roi_labels)
    patients = np.asarray([imc_tm.roi_to_patient(r) for r in roi_labels])
    print("\nembedding matrix {} from {} ROIs / {} patients".format(
        feats.shape, len(set(roi_labels)), len(set(patients))))

    pca = SmartPca(preprocess_strategy="z_score")
    umap = SmartUmap(n_neighbors=args.umap_neighbors, preprocess_strategy="raw",
                     n_components=2, min_dist=0.5, metric="euclidean")
    emb_pca = pca.fit_transform(torch.tensor(feats), n_components=0.95)
    emb_umap = imc_tm.to_numpy(umap.fit_transform(emb_pca))
    print("PCA -> {} components; UMAP -> {}".format(imc_tm.to_numpy(emb_pca).shape[-1], emb_umap.shape))

    np.savez_compressed(
        out_dir / "patch_features.npz",
        features=feats, umap=emb_umap, composition=comps, moran=morans,
        roi=roi_labels, patient=patients, channels=np.asarray(channels))

    panels = [("max Morans I", morans, "viridis"), ("patient", None, None)]
    panels += [(c + " fraction", comps[:, k], "viridis") for k, c in enumerate(channels)]

    ncols = 4
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4.2 * ncols, 4.0 * nrows))
    axes = np.asarray(axes).reshape(-1)
    for ax, (title, values, cmap) in zip(axes, panels):
        if values is None:  # categorical: patient
            codes = {p: j for j, p in enumerate(sorted(set(patients)))}
            ax.scatter(emb_umap[:, 0], emb_umap[:, 1],
                       c=[codes[p] for p in patients], s=2, cmap="tab20", linewidths=0)
        else:
            sc = ax.scatter(emb_umap[:, 0], emb_umap[:, 1], c=values, s=2, cmap=cmap, linewidths=0)
            fig.colorbar(sc, ax=ax, fraction=0.046)
        ax.set_title(title, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes[len(panels):]:
        ax.axis("off")
    fig.suptitle("TissueMosaic patch embeddings, JHH IMC cohort "
                 "({} patches, {} ROIs)".format(len(feats), len(set(roi_labels))), fontsize=13)
    fig.tight_layout()
    fig.savefig(out_dir / "patch_umap.png", dpi=130)
    print("\nwrote {} and {}".format(out_dir / "patch_umap.png", out_dir / "patch_features.npz"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
