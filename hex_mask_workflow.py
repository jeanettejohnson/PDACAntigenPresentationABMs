import argparse
import json
import math
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, shape


def load_annotations_gdf(geojson_path: Path) -> gpd.GeoDataFrame:
    """Load a GeoJSON into one row per annotation feature."""
    with geojson_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    features = data.get("features", [])
    records = []
    for feat in features:
        geom = shape(feat["geometry"])
        props = feat.get("properties", {}).copy()
        props["annotation_id"] = feat.get("id")
        props["geometry"] = geom
        records.append(props)

    gdf = gpd.GeoDataFrame(records, geometry="geometry")
    gdf["annotation_area"] = gdf.geometry.area
    gdf["annotation_perimeter"] = gdf.geometry.length
    return gdf


def make_hexagon(cx: float, cy: float, radius: float) -> Polygon:
    """Create a pointy-top hexagon centered at (cx, cy)."""
    pts = []
    for k in range(6):
        angle = math.radians(60 * k + 30)
        x = cx + radius * math.cos(angle)
        y = cy + radius * math.sin(angle)
        pts.append((x, y))
    return Polygon(pts)


def build_hex_grid(bounds, radius: float) -> gpd.GeoDataFrame:
    """Create a planar hex grid covering bounds=(minx, miny, maxx, maxy)."""
    minx, miny, maxx, maxy = bounds

    dx = math.sqrt(3) * radius
    dy = 1.5 * radius

    hexes = []
    row = 0
    y = miny - 2 * radius
    while y <= maxy + 2 * radius:
        x_offset = dx / 2 if row % 2 else 0
        x = minx - 2 * dx + x_offset
        while x <= maxx + 2 * dx:
            hexes.append(make_hexagon(x, y, radius))
            x += dx
        y += dy
        row += 1

    hex_gdf = gpd.GeoDataFrame({"geometry": hexes}, geometry="geometry")
    hex_gdf["hex_id"] = np.arange(len(hex_gdf))
    hex_gdf["hex_area"] = hex_gdf.geometry.area
    return hex_gdf


def mask_hexes_by_annotations(
    hex_gdf: gpd.GeoDataFrame,
    ann_gdf: gpd.GeoDataFrame,
    predicate: str = "intersects",
) -> gpd.GeoDataFrame:
    """Keep hexes that satisfy a spatial predicate against annotations."""
    joined = gpd.sjoin(
        hex_gdf,
        ann_gdf[["annotation_id", "geometry"]],
        how="inner",
        predicate=predicate,
    )
    masked = joined.drop_duplicates(subset=["hex_id"]).copy()
    return masked[["hex_id", "hex_area", "geometry"]]


def plot_workflow(
    ann_gdf: gpd.GeoDataFrame,
    hex_gdf: gpd.GeoDataFrame,
    masked_hexes: gpd.GeoDataFrame,
    output_png: Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)

    ann_gdf.plot(ax=axes[0], color="#2a9d8f", edgecolor="black", linewidth=0.3)
    axes[0].set_title("Annotations")
    axes[0].set_aspect("equal")

    hex_gdf.boundary.plot(ax=axes[1], color="#999999", linewidth=0.2)
    ann_gdf.boundary.plot(ax=axes[1], color="#d62828", linewidth=0.6)
    axes[1].set_title("Hex Grid Overlay")
    axes[1].set_aspect("equal")

    hex_gdf.boundary.plot(ax=axes[2], color="#dddddd", linewidth=0.2)
    masked_hexes.plot(ax=axes[2], color="#264653", edgecolor="white", linewidth=0.2)
    ann_gdf.boundary.plot(ax=axes[2], color="#e76f51", linewidth=0.5)
    axes[2].set_title("Masked Hexes")
    axes[2].set_aspect("equal")

    for ax in axes:
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.tick_params(axis="both", labelsize=8)
        ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.6)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="GeoJSON annotation -> hex mask workflow")
    parser.add_argument(
        "--geojson",
        type=Path,
        required=True,
        help="Path to annotation GeoJSON",
    )
    parser.add_argument(
        "--hex-radius",
        type=float,
        default=8.0,
        help="Hex radius in the same units as annotation coordinates",
    )
    parser.add_argument(
        "--predicate",
        choices=["intersects", "within", "contains", "overlaps", "touches"],
        default="intersects",
        help="Spatial predicate used to keep hexes",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("outputs"),
        help="Directory for outputs",
    )
    args = parser.parse_args()

    ann_gdf = load_annotations_gdf(args.geojson)
    hex_gdf = build_hex_grid(ann_gdf.total_bounds, radius=args.hex_radius)
    masked_hexes = mask_hexes_by_annotations(hex_gdf, ann_gdf, predicate=args.predicate)

    args.outdir.mkdir(parents=True, exist_ok=True)

    annotations_csv = args.outdir / "annotations_table.csv"
    masked_hexes_csv = args.outdir / "masked_hexes.csv"
    final_points_csv = args.outdir / "hex_centers_duct_filler.csv"
    figure_png = args.outdir / "annotation_hex_workflow.png"

    ann_df = pd.DataFrame(ann_gdf.drop(columns="geometry"))
    ann_df.to_csv(annotations_csv, index=False)
    masked_hexes.drop(columns="geometry").to_csv(masked_hexes_csv, index=False)

    centers = masked_hexes.geometry.centroid
    output_df = pd.DataFrame(
        {
            "x": centers.x,
            "y": centers.y,
            "z": 0,
            "type": "duct_filler",
        }
    )
    output_df.to_csv(final_points_csv, index=False)

    plot_workflow(ann_gdf, hex_gdf, masked_hexes, figure_png)

    print(f"Annotations: {len(ann_gdf)}")
    print(f"Hexes total: {len(hex_gdf)}")
    print(f"Hexes after mask ({args.predicate}): {len(masked_hexes)}")
    print(f"Saved annotation table: {annotations_csv}")
    print(f"Saved masked hex table: {masked_hexes_csv}")
    print(f"Saved final x,y,z,type table: {final_points_csv}")
    print(f"Saved visualization: {figure_png}")


if __name__ == "__main__":
    main()
