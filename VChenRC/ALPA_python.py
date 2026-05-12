"""Usage: python ALPA_python.py --batch  or  python ALPA_python.py --dem <dem.tif> --landslides <landslide.shp> --output-dir <output_dir>."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import math
import time
import warnings

import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from shapely.errors import GEOSException
from shapely.geometry import GeometryCollection, LineString, MultiLineString, Point, Polygon, mapping
from shapely.ops import split as shapely_split, triangulate as shapely_triangulate

matplotlib.use('Agg')


PNT_COLUMNS = [
    "ID",
    "IDB",
    "Cx",
    "Cy",
    "Cz",
    "grp",
    "Nx",
    "Ny",
    "Nz",
    "zip",
    "Dx",
    "Dy",
    "Dz",
    "SMG",
    "SMS",
    "_row",
    "_col",
    "_a",
    "_b",
    "_c",
]


@dataclass
class PolygonParameters:
    name: str
    min_group_point_count: int = 3
    min_group_anchor_distance: float = 0.0
    min_strip_horizontal_length: float = 30.0
    min_end_ratio_initial: float = 0.0
    min_end_ratio_distal: float = 0.0


def ALPA(
    InDEM: str | Path,
    InLandslides: str | Path,
    InParameters: str | Path | int = 0,
    InCfgEndAnchorsInt: Any = 1,
    InCfgEndAnchorsDit: Any = 1,
    InCfgEndStripsInt: Any = 1,
    InCfgEndStripsDit: Any = 1,
    InOutputHierarchical: bool = False,
    output_dir: str | Path | None = None,
) -> dict[str, Any]:
    dem_path = Path(InDEM)
    landslide_path = Path(InLandslides)
    output_root = Path(output_dir) if output_dir else landslide_path.parent
    output_root.mkdir(parents=True, exist_ok=True)

    landslides = gpd.read_file(landslide_path)
    landslides = landslides.explode(index_parts=False).reset_index(drop=True)
    landslides = landslides[landslides.geometry.notna()].copy()
    landslides = landslides[landslides.geometry.geom_type.isin(["Polygon", "MultiPolygon"])]
    landslides = landslides.reset_index(drop=True)

    parameters = _load_parameters(InParameters, len(landslides))
    anchor_int_cfg = _load_config_features(InCfgEndAnchorsInt)
    anchor_dit_cfg = _load_config_features(InCfgEndAnchorsDit)
    strip_int_cfg = _load_config_features(InCfgEndStripsInt)
    strip_dit_cfg = _load_config_features(InCfgEndStripsDit)
    failed_cases_csv = output_root / f"{landslide_path.stem}_failed_cases.csv"

    results: list[dict[str, Any]] = []
    combined_paths: list[gpd.GeoDataFrame] = []
    combined_path_records: list[pd.DataFrame] = []
    failed_cases: list[dict[str, Any]] = []
    total_started_at = time.perf_counter()

    print("Running ALPA with:")
    print(f"  DEM: {dem_path}")
    print(f"  Landslides: {landslide_path}")
    print(f"  Output dir: {output_root}")
    print(f"  Case count: {len(landslides)}\n")

    with rasterio.open(dem_path) as src:
        raster_crs = src.crs
        if landslides.crs and raster_crs and landslides.crs != raster_crs:
            landslides = landslides.to_crs(raster_crs)

        for index, row in landslides.iterrows():
            case_started_at = time.perf_counter()
            params = _parameters_for_index(parameters, index)
            case_name = params.name
            polygon = _ensure_polygon(row.geometry)
            if polygon is None or polygon.is_empty:
                failed_cases.append(
                    {
                        "case_index": int(index + 1),
                        "case_id": str(case_name),
                        "reason": "Invalid or empty polygon geometry",
                    }
                )
                elapsed_case = time.perf_counter() - case_started_at
                print(f"Process the {index + 1}/{len(landslides)} landslide failed.")
                print(f"{elapsed_case} seconds used.")
                continue

            p_end_anchor_int = _config_for_index(InCfgEndAnchorsInt, anchor_int_cfg, index)
            p_end_anchor_dit = _config_for_index(InCfgEndAnchorsDit, anchor_dit_cfg, index)
            p_end_strip_int = _config_for_index(InCfgEndStripsInt, strip_int_cfg, index)
            p_end_strip_dit = _config_for_index(InCfgEndStripsDit, strip_dit_cfg, index)
            try:
                result = f_lasld_split(
                    src=src,
                    polygon=polygon,
                    name=params.name,
                    min_group_point_count=params.min_group_point_count,
                    min_group_anchor_distance=params.min_group_anchor_distance,
                    min_strip_horizontal_length=params.min_strip_horizontal_length,
                    min_end_ratio_initial=params.min_end_ratio_initial,
                    min_end_ratio_distal=params.min_end_ratio_distal,
                    end_anchors_initial=p_end_anchor_int,
                    end_anchors_distal=p_end_anchor_dit,
                    end_strips_initial=p_end_strip_int,
                    end_strips_distal=p_end_strip_dit,
                    output_hierarchical=InOutputHierarchical,
                    output_dir=output_root,
                )
            except Exception as exc:
                failed_cases.append(
                    {
                        "case_index": int(index + 1),
                        "case_id": str(case_name),
                        "reason": str(exc),
                    }
                )
                elapsed_case = time.perf_counter() - case_started_at
                print(f"Process the {index + 1}/{len(landslides)} landslide failed.")
                print(f"{elapsed_case} seconds used.")
                continue

            if not _has_valid_path_gdf(result.get("path_gdf")):
                failed_cases.append(
                    {
                        "case_index": int(index + 1),
                        "case_id": str(case_name),
                        "reason": "Processing completed but no valid final profile line was produced",
                    }
                )
                elapsed_case = time.perf_counter() - case_started_at
                print(f"Process the {index + 1}/{len(landslides)} landslide failed (no valid final profile line).")
                print(f"{elapsed_case} seconds used.")
                continue
            results.append(result)
            combined_paths.append(result["path_gdf"])
            combined_path_records.append(pd.DataFrame(result["path_gdf"].drop(columns="geometry")))

            elapsed_case = time.perf_counter() - case_started_at
            print(f"Process the {index + 1}/{len(landslides)} landslide done.")
            print(f"{elapsed_case} seconds used.")

    merged_path_gdf = None
    merged_path_df = None
    if combined_paths:
        merged_path_gdf = pd.concat(combined_paths, ignore_index=True)
        merged_path_df = pd.concat(combined_path_records, ignore_index=True)
        base_name = landslide_path.stem
        path_prefix = output_root / f"{base_name}_paths2d"
        _safe_write_gdf(merged_path_gdf, path_prefix)
        _write_excel(merged_path_df, path_prefix.with_suffix(".xlsx"))

    failed_df = pd.DataFrame(failed_cases, columns=["case_index", "case_id", "reason"])
    failed_df.to_csv(failed_cases_csv, index=False, encoding="utf-8-sig")

    elapsed_total = time.perf_counter() - total_started_at
    print("\nRun Summary:")
    print(f"  Total cases: {len(landslides)}")
    print(f"  Success cases: {len(results)}")
    print(f"  Failed cases: {len(failed_cases)}")
    print(f"  Elapsed time (s): {elapsed_total:.3f}")
    print(f"  Failed cases CSV: {failed_cases_csv}")
    if failed_cases:
        print("  Failed case IDs: " + ", ".join(str(item["case_id"]) for item in failed_cases))

    return {
        "results": results,
        "paths": merged_path_gdf,
        "path_table": merged_path_df,
        "total_cases": len(landslides),
        "success_cases": len(results),
        "failed_cases": failed_df,
        "failed_cases_csv": failed_cases_csv,
        "elapsed_seconds": elapsed_total,
    }


def f_lasld_split(
    src: rasterio.io.DatasetReader,
    polygon: Polygon,
    name: str,
    min_group_point_count: int,
    min_group_anchor_distance: float,
    min_strip_horizontal_length: float,
    min_end_ratio_initial: float,
    min_end_ratio_distal: float,
    end_anchors_initial: Any,
    end_anchors_distal: Any,
    end_strips_initial: Any,
    end_strips_distal: Any,
    output_hierarchical: bool,
    output_dir: Path,
) -> dict[str, Any]:
    pnts = f_pnts_read(src, polygon)
    if pnts.empty or int(pnts["IDB"].max()) == 0:
        pnts.loc[:, "SMS"] = "RPF"
        _write_points_products(pnts, src.crs, output_dir / f"{name}_pnts")
        raise RuntimeError("No valid boundary points found in DEM raster (RPF)")

    _write_points_products(pnts, src.crs, output_dir / f"{name}_pnts")

    if _is_horizontal_plane(pnts):
        pnts.loc[:, "SMS"] = "LSH"
        _write_points_products(pnts, src.crs, output_dir / f"{name}_pnts_grp")
        raise RuntimeError("Landslide polygon appears to be on a horizontal plane (LSH)")

    total_count = len(pnts)
    idb_max = int(pnts["IDB"].max())
    upstream_point = pnts.loc[[pnts["zip"].idxmax()]].copy()
    groups = [pnts.copy()]
    upstream_points = [upstream_point]
    group_markers = [0]
    hierarchy_snapshots: list[list[pd.DataFrame]] = []

    while True:
        hierarchy_snapshots.append([group.copy() for group in groups])
        if min(group_markers) == 1:
            break

        new_groups: list[pd.DataFrame] = []
        new_upstream_points: list[pd.DataFrame] = []
        new_markers: list[int] = []
        grp_id = 1
        grp_max = int(max(group["grp"].iloc[0] for group in groups))

        for idx, group in enumerate(groups):
            pntu = upstream_points[idx]
            marker = group_markers[idx]
            if marker == 1:
                group = group.copy()
                group.loc[:, "grp"] = grp_id
                new_groups.append(group)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            bool_initial = int(group["grp"].iloc[0]) == 1
            bool_distal = int(group["grp"].iloc[0]) == grp_max
            end_ratio = len(group) / total_count
            min_end_ratio = min_end_ratio_initial if bool_initial else min_end_ratio_distal
            if (bool_initial or bool_distal) and end_ratio < min_end_ratio:
                stopped = _mark_stopped_group(group, grp_id, "SGR")
                new_groups.append(stopped)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            intercept = f_pnts_intercept(group)
            upper, lower = f_pnts_split(group, intercept)
            if upper.empty or lower.empty:
                stopped = _mark_stopped_group(group, grp_id, "SGP")
                new_groups.append(stopped)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            if int(pntu["ID"].iloc[0]) not in set(upper["ID"]):
                upper, lower = lower, upper

            upper, lower, topology_ok = _reassign_clusters(group, upper, lower, pntu)
            if not topology_ok:
                stopped = _mark_stopped_group(group, grp_id, "SGT")
                new_groups.append(stopped)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            if len(upper) < min_group_point_count or len(lower) < min_group_point_count:
                stopped = _mark_stopped_group(group, grp_id, "SGP")
                new_groups.append(stopped)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            upper = f_pnts_plane(upper)
            lower = f_pnts_plane(lower)
            if upper is None or lower is None:
                stopped = _mark_stopped_group(group, grp_id, "SGL")
                new_groups.append(stopped)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            if _is_horizontal_plane(upper) or _is_horizontal_plane(lower):
                stopped = _mark_stopped_group(group, grp_id, "SGH")
                new_groups.append(stopped)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            if not _validate_boundary_side_count(upper, idb_max, initial=bool_initial, distal=False):
                stopped = _mark_stopped_group(group, grp_id, "SGB")
                new_groups.append(stopped)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            if not _validate_boundary_side_count(lower, idb_max, initial=False, distal=bool_distal):
                stopped = _mark_stopped_group(group, grp_id, "SGB")
                new_groups.append(stopped)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            if bool_initial:
                upper_bnd = upper[upper["IDB"] != 0].copy()
                upper_bnd_ascending = f_pnts_bnd_ascending(upper_bnd)
                if upper_bnd_ascending.empty:
                    stopped = _mark_stopped_group(group, grp_id, "SGE")
                    new_groups.append(stopped)
                    new_upstream_points.append(pntu)
                    new_markers.append(1)
                    grp_id += 1
                    continue
                pnts0_upper = upper_bnd_ascending.iloc[[0, len(upper_bnd_ascending) - 1]][["Cx", "Cy"]].to_numpy(dtype=float)
                end_upper = _pnts_bnd_end_even(upper, pnts0_upper)
                if end_upper[0] is None or end_upper[1] is None or end_upper[2] is None:
                    stopped = _mark_stopped_group(group, grp_id, "SGE")
                    new_groups.append(stopped)
                    new_upstream_points.append(pntu)
                    new_markers.append(1)
                    grp_id += 1
                    continue

            if bool_distal:
                lower_bnd = lower[lower["IDB"] != 0].copy()
                lower_bnd_ascending = f_pnts_bnd_ascending(lower_bnd)
                if lower_bnd_ascending.empty:
                    stopped = _mark_stopped_group(group, grp_id, "SGE")
                    new_groups.append(stopped)
                    new_upstream_points.append(pntu)
                    new_markers.append(1)
                    grp_id += 1
                    continue
                pnts0_lower = lower_bnd_ascending.iloc[[len(lower_bnd_ascending) - 1, 0]][["Cx", "Cy"]].to_numpy(dtype=float)
                end_lower = _pnts_bnd_end_even(lower, pnts0_lower)
                if end_lower[0] is None or end_lower[1] is None or end_lower[2] is None:
                    stopped = _mark_stopped_group(group, grp_id, "SGE")
                    new_groups.append(stopped)
                    new_upstream_points.append(pntu)
                    new_markers.append(1)
                    grp_id += 1
                    continue

            candidate_groups = new_groups + [upper.copy(), lower.copy()] + [item.copy() for item in groups[idx + 1 :]]
            candidate_groups = _renumber_groups(candidate_groups)
            path_candidate = f_pnts_path(
                candidate_groups,
                idb_max,
                polygon,
                end_anchors_initial,
                end_anchors_distal,
                end_strips_initial,
                end_strips_distal,
                False,
            )
            if path_candidate is None:
                stopped = _mark_stopped_group(group, grp_id, "SGC")
                new_groups.append(stopped)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            anchors = path_candidate[0]
            distances = f_pnts_distances(anchors[["Cx", "Cy", "Cz"]])
            distances = distances[(distances > 0) & np.isfinite(distances)]
            if distances.size and float(distances.min()) <= float(min_group_anchor_distance):
                stopped = _mark_stopped_group(group, grp_id, "SGD")
                new_groups.append(stopped)
                new_upstream_points.append(pntu)
                new_markers.append(1)
                grp_id += 1
                continue

            parent_dip = group[["Dx", "Dy", "Dz"]].iloc[0].to_numpy(dtype=float)
            upper.loc[:, "grp"] = grp_id
            lower.loc[:, "grp"] = grp_id + 1
            new_groups.extend([upper, lower])
            new_upstream_points.append(pntu)
            new_upstream_points.append(f_pnts_MinDist_intercept(lower, parent_dip, intercept))
            new_markers.extend([0, 0])
            grp_id += 2

        groups = _renumber_groups(new_groups)
        upstream_points = new_upstream_points
        group_markers = new_markers

    groups = _merge_small_end_groups(groups, total_count, min_end_ratio_initial, min_end_ratio_distal)
    groups = _renumber_groups(groups)

    path_result = f_pnts_path(
        groups,
        idb_max,
        polygon,
        end_anchors_initial,
        end_anchors_distal,
        end_strips_initial,
        end_strips_distal,
        True,
    )
    if path_result is None:
        pnts_split = f_list2pnts(groups)
        pnts_split.loc[:, "SMS"] = "SMG"
        pnts_split.loc[:, "SMS"] = "SMO"
        _write_points_products(pnts_split, src.crs, output_dir / f"{name}_pnts_grps")
        raise RuntimeError("Path optimization returned no result after group partitioning (SMO)")

    groups = path_result[1].get("ordered_groups", groups)
    groups = _renumber_groups([group.copy() for group in groups])
    pnts_split = f_list2pnts(groups)
    pnts_split.loc[:, "SMS"] = "SMG"
    anchors = path_result[0]
    strips_package = f_strips_package(
        groups,
        polygon,
        anchors,
        path_result[1],
        end_strips_initial,
        end_strips_distal,
        min_strip_horizontal_length,
        src,
    )

    _write_points_products(pnts_split, src.crs, output_dir / f"{name}_pnts_grps")
    _write_points_products(anchors, src.crs, output_dir / f"{name}_pnts_grps_anchors")

    strip_gdf = strips_package.get("strip_gdf")
    strip_df = strips_package.get("strip_df")
    strip_nodes = strips_package.get("strip_nodes")
    triangle_gdf = strips_package.get("triangle_gdf")
    triangle_df = strips_package.get("triangle_df")
    path_gdf = f_pnts_save_path(
        anchors,
        strips_package.get("paras_strip"),
        name,
        min_group_point_count,
        min_group_anchor_distance,
        min_strip_horizontal_length,
        min_end_ratio_initial,
        min_end_ratio_distal,
        end_anchors_initial,
        end_anchors_distal,
        end_strips_initial,
        end_strips_distal,
        src.crs,
        output_dir / f"{name}_pnts_grps_strips_nodes_pln2d",
    )

    if strip_df is not None and not strip_df.empty:
        _write_excel(strip_df, output_dir / f"{name}_pnts_grps_strips.xlsx")
    if strip_gdf is not None and not strip_gdf.empty:
        _safe_write_gdf(strip_gdf, output_dir / f"{name}_pnts_grps_strips")
    if strip_nodes is not None and not strip_nodes.empty:
        _write_excel(strip_nodes, output_dir / f"{name}_pnts_grps_strips_nodes.xlsx")
        strip_nodes_gdf = gpd.GeoDataFrame(
            strip_nodes.copy(),
            geometry=gpd.points_from_xy(strip_nodes["Cx"], strip_nodes["Cy"]),
            crs=src.crs,
        )
        _safe_write_gdf(strip_nodes_gdf, output_dir / f"{name}_pnts_grps_strips_nodes")
    if triangle_df is not None and not triangle_df.empty:
        _write_excel(triangle_df, output_dir / f"{name}_pnts_grps_strips_triangles.xlsx")
    if triangle_gdf is not None and not triangle_gdf.empty:
        _safe_write_gdf(triangle_gdf, output_dir / f"{name}_pnts_grps_strips_triangles")

    if output_hierarchical:
        for level, snapshot in enumerate(hierarchy_snapshots, start=1):
            snapshot_pnts = f_list2pnts(_renumber_groups([item.copy() for item in snapshot]))
            _write_points_products(snapshot_pnts, src.crs, output_dir / f"{name}_hier_{level:03d}")

    return {
        "name": name,
        "points": pnts_split,
        "anchors": anchors,
        "path_gdf": path_gdf,
        "strip_gdf": strip_gdf,
        "strip_df": strip_df,
        "strip_nodes": strip_nodes,
        "triangle_gdf": triangle_gdf,
        "triangle_df": triangle_df,
    }


def f_polygon_rasterize(src: rasterio.io.DatasetReader, polygon: Polygon) -> tuple[Polygon, np.ma.MaskedArray, Any]:
    data, transform = mask(src, [mapping(polygon)], crop=True, filled=False)
    return polygon, data[0], transform


def f_pnts_read(src: rasterio.io.DatasetReader, polygon: Polygon) -> pd.DataFrame:
    _, raster_masked, transform = f_polygon_rasterize(src, polygon)
    valid = ~np.ma.getmaskarray(raster_masked)
    rows, cols = np.where(valid)
    if rows.size == 0:
        return pd.DataFrame(columns=PNT_COLUMNS)

    xs, ys = rasterio.transform.xy(transform, rows, cols, offset="center")
    zs = raster_masked.data[rows, cols].astype(float)

    records = []
    for idx, (x, y, z, row, col) in enumerate(zip(xs, ys, zs, rows, cols), start=1):
        records.append(
            {
                "ID": idx,
                "IDB": 0,
                "Cx": float(x),
                "Cy": float(y),
                "Cz": float(z),
                "grp": 1,
                "Nx": 0.0,
                "Ny": 0.0,
                "Nz": 0.0,
                "zip": 0.0,
                "Dx": 0.0,
                "Dy": 0.0,
                "Dz": 0.0,
                "SMG": "SMS",
                "SMS": "NA",
                "_row": int(row),
                "_col": int(col),
                "_a": 0.0,
                "_b": 0.0,
                "_c": 0.0,
            }
        )

    pnts = pd.DataFrame.from_records(records, columns=PNT_COLUMNS)
    boundary_mask = np.zeros_like(valid, dtype=bool)
    for row, col in zip(rows, cols):
        for drow, dcol in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nrow = row + drow
            ncol = col + dcol
            if nrow < 0 or nrow >= valid.shape[0] or ncol < 0 or ncol >= valid.shape[1] or not valid[nrow, ncol]:
                boundary_mask[row, col] = True
                break

    boundary_rows, boundary_cols = np.where(boundary_mask)
    if boundary_rows.size:
        boundary_xs, boundary_ys = rasterio.transform.xy(transform, boundary_rows, boundary_cols, offset="center")
        boundary_projects = [
            polygon.boundary.project(Point(float(coord_x), float(coord_y)))
            for coord_x, coord_y in zip(boundary_xs, boundary_ys)
        ]
        boundary_order = np.argsort(np.array(boundary_projects, dtype=float))
        for sequence, ordered_idx in enumerate(boundary_order, start=1):
            row = int(boundary_rows[int(ordered_idx)])
            col = int(boundary_cols[int(ordered_idx)])
            match = pnts.index[(pnts["_row"] == row) & (pnts["_col"] == col)]
            if len(match) == 1:
                pnts.loc[match[0], "IDB"] = sequence

    plane = f_pnts_plane(pnts)
    if plane is None:
        return pnts
    idsb = plane.index[plane["IDB"] != 0]
    if len(idsb) > 0:
        max_idb = int(plane.loc[idsb, "IDB"].max())
        plane.loc[idsb, "IDB"] = max_idb + 1 - plane.loc[idsb, "IDB"]
    return plane


def f_pnts_plane(pnts: pd.DataFrame) -> pd.DataFrame | None:
    if len(pnts) < 3:
        return None
    x = pnts["Cx"].to_numpy(dtype=float)
    y = pnts["Cy"].to_numpy(dtype=float)
    z = pnts["Cz"].to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(pnts)), x, y])
    rank = np.linalg.matrix_rank(design)
    if rank < 3:
        return None
    coeffs, _, _, _ = np.linalg.lstsq(design, z, rcond=None)
    c0, ax, by = coeffs
    normal = f_vector_unit(np.array([-ax, -by, 1.0], dtype=float))
    if np.isnan(normal).any():
        return None
    pnts = pnts.copy()
    pnts.loc[:, "Nx"] = normal[0]
    pnts.loc[:, "Ny"] = normal[1]
    pnts.loc[:, "Nz"] = normal[2]
    pnts.loc[:, "zip"] = c0 + ax * x + by * y
    dip = f_pnts_dip(pnts)
    pnts.loc[:, "Dx"] = dip[0]
    pnts.loc[:, "Dy"] = dip[1]
    pnts.loc[:, "Dz"] = dip[2]
    pnts.loc[:, "_a"] = ax
    pnts.loc[:, "_b"] = by
    pnts.loc[:, "_c"] = c0
    return pnts


def f_pnts_dip(pnts: pd.DataFrame) -> np.ndarray:
    normal = np.array([pnts.iloc[0]["Nx"], pnts.iloc[0]["Ny"], pnts.iloc[0]["Nz"]], dtype=float)
    if np.isnan(normal).any():
        return np.array([np.nan, np.nan, np.nan], dtype=float)
    if np.isclose(normal[0], 0.0) and np.isclose(normal[1], 0.0):
        existing = pnts[["Dx", "Dy", "Dz"]].iloc[0].to_numpy(dtype=float)
        return existing if np.any(existing) else np.zeros(3, dtype=float)
    plumb = np.array([0.0, 0.0, -1.0], dtype=float)
    strike = f_vector_unit(np.cross(normal, plumb))
    if np.allclose(strike, 0.0):
        existing = pnts[["Dx", "Dy", "Dz"]].iloc[0].to_numpy(dtype=float)
        return existing if np.any(existing) else np.zeros(3, dtype=float)
    dip = f_vector_unit(np.cross(strike, normal))
    return dip


def f_vector_unit(vector: np.ndarray) -> np.ndarray:
    magnitude = float(np.linalg.norm(vector))
    if magnitude == 0.0:
        return np.zeros_like(vector, dtype=float)
    return vector / magnitude


def f_pnts_intercept(pnts: pd.DataFrame) -> float:
    dip = pnts[["Dx", "Dy", "Dz"]].iloc[0].to_numpy(dtype=float)
    coords = pnts[["Cx", "Cy", "Cz"]].to_numpy(dtype=float)
    intercepts = -coords @ dip
    return float((intercepts.min() + intercepts.max()) / 2.0)


def f_pnts_split(pnts: pd.DataFrame, intercept: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    dip = pnts[["Dx", "Dy", "Dz"]].iloc[0].to_numpy(dtype=float)
    coords = pnts[["Cx", "Cy", "Cz"]].to_numpy(dtype=float)
    direction = coords @ dip + intercept
    upper = pnts.loc[direction <= 0].copy()
    lower = pnts.loc[direction > 0].copy()
    return upper, lower


def f_pnts_clustering(pnts: pd.DataFrame) -> list[pd.DataFrame]:
    if pnts.empty:
        return []
    index_by_cell = {
        (int(row["_row"]), int(row["_col"])): int(row["ID"])
        for _, row in pnts.iterrows()
    }
    remaining = set(index_by_cell)
    clusters: list[pd.DataFrame] = []
    offsets = [
        (-1, 0),
        (0, -1),
        (0, 1),
        (1, 0),
    ]
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        stack = [start]
        cluster_ids: list[int] = []
        while stack:
            current = stack.pop()
            cluster_ids.append(index_by_cell[current])
            for dr, dc in offsets:
                neighbor = (current[0] + dr, current[1] + dc)
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
        clusters.append(pnts[pnts["ID"].isin(cluster_ids)].copy())
    return clusters


def f_pnts_bnd_ascending(pnts_bnd: pd.DataFrame) -> pd.DataFrame:
    if pnts_bnd.empty:
        return pnts_bnd.copy()
    ordered = pnts_bnd.sort_values("IDB").reset_index(drop=True)
    start_idx = 0
    for idx in range(1, len(ordered)):
        if int(ordered.loc[idx, "IDB"]) - int(ordered.loc[idx - 1, "IDB"]) > 1:
            start_idx = idx
            break
    if start_idx == 0:
        return ordered
    return pd.concat([ordered.iloc[start_idx:], ordered.iloc[:start_idx]], ignore_index=True)


def f_pnts_bnd_sides(pnts_bnd: pd.DataFrame, idb_max: int) -> list[pd.DataFrame]:
    if pnts_bnd.empty:
        return []
    if len(pnts_bnd) == 1:
        return [pnts_bnd.copy()]
    ordered = pnts_bnd.sort_values("IDB").reset_index(drop=True)
    sides: list[pd.DataFrame] = []
    current = [ordered.iloc[[0]].copy()]
    for idx in range(1, len(ordered)):
        step = int(ordered.loc[idx, "IDB"]) - int(ordered.loc[idx - 1, "IDB"])
        if step == 1:
            current.append(ordered.iloc[[idx]].copy())
        else:
            sides.append(pd.concat(current, ignore_index=True))
            current = [ordered.iloc[[idx]].copy()]
    last_side = pd.concat(current, ignore_index=True)
    if not sides:
        sides.append(last_side)
    elif int(ordered.loc[0, "IDB"]) == 1 and int(ordered.loc[len(ordered) - 1, "IDB"]) == idb_max:
        sides[0] = pd.concat([last_side, sides[0]], ignore_index=True)
    else:
        sides.append(last_side)
    return sides


def _config_xy(config: Any) -> np.ndarray | None:
    if isinstance(config, np.ndarray) and config.size >= 2:
        return config[:2].astype(float)
    if isinstance(config, (list, tuple, np.ndarray)) and len(config) == 2:
        return np.array(config, dtype=float)
    return None


def _boundary_side_corners(list_bnd_sides: list[pd.DataFrame]) -> np.ndarray | None:
    if not list_bnd_sides:
        return None
    corners: list[np.ndarray] = []
    for side in list_bnd_sides:
        if side is None or side.empty:
            return None
        corners.append(side.iloc[0][["Cx", "Cy"]].to_numpy(dtype=float))
        corners.append(side.iloc[-1][["Cx", "Cy"]].to_numpy(dtype=float))
    return np.vstack(corners) if corners else None


def _boundary_side_corners_check(list_bnd_sides: list[pd.DataFrame]) -> tuple[bool, bool]:
    corners = _boundary_side_corners(list_bnd_sides)
    if corners is None or len(corners) < 4:
        return False, False
    diff41 = not np.allclose(corners[3], corners[0], atol=1e-12)
    diff32 = not np.allclose(corners[2], corners[1], atol=1e-12)
    return diff41, diff32


def _boundary_arc_between(ordered: pd.DataFrame, start_idx: int, end_idx: int) -> pd.DataFrame | None:
    if start_idx < 0 or end_idx < 0 or ordered.empty:
        return None
    if start_idx <= end_idx:
        return ordered.iloc[start_idx : end_idx + 1].copy().reset_index(drop=True)
    return pd.concat([ordered.iloc[start_idx:], ordered.iloc[: end_idx + 1]], ignore_index=True)


def f_pnts_bnd_sides_centerLs(list_bnd_sides: list[pd.DataFrame]) -> np.ndarray | None:
    if not list_bnd_sides:
        return None
    corners = _boundary_side_corners(list_bnd_sides)
    if corners is None:
        return None
    if len(corners) < 4:
        return None
    center1 = (corners[3] + corners[0]) / 2.0
    center2 = (corners[2] + corners[1]) / 2.0
    center_mid = (center1 + center2) / 2.0
    return np.vstack([center1, center2, center_mid])


def f_pnts_path(
    list_pnts: list[pd.DataFrame],
    idb_max: int,
    polygon: Polygon,
    end_anchors_initial: Any,
    end_anchors_distal: Any,
    end_strips_initial: Any,
    end_strips_distal: Any,
    bool_direction: bool,
) -> tuple[pd.DataFrame, dict[str, Any]] | None:
    if not list_pnts:
        return None
    groups = _renumber_groups([group.copy() for group in list_pnts])
    count_pnts = len(groups)

    if count_pnts == 1:
        group = groups[0]
        pnts_bnd = group[group["IDB"] != 0].copy()
        pnts_bnd_ascending = f_pnts_bnd_ascending(pnts_bnd)
        if pnts_bnd_ascending.empty:
            return None

        start_anchor_xy = _config_xy(end_anchors_initial)
        end_anchor_xy = _config_xy(end_anchors_distal)
        if start_anchor_xy is not None and end_anchor_xy is not None:
            anchor_mid_xy = (start_anchor_xy + end_anchor_xy) / 2.0
            anchors = pd.DataFrame(
                [
                    {
                        "Cx": float(start_anchor_xy[0]),
                        "Cy": float(start_anchor_xy[1]),
                        "Cz": f_pnts_zip(group, float(start_anchor_xy[0]), float(start_anchor_xy[1])),
                        "type": "Anchor: Initial group anchor",
                        "grp": 0,
                    },
                    {
                        "Cx": float(anchor_mid_xy[0]),
                        "Cy": float(anchor_mid_xy[1]),
                        "Cz": f_pnts_zip(group, float(anchor_mid_xy[0]), float(anchor_mid_xy[1])),
                        "type": "Anchor: Group center",
                        "grp": 1,
                    },
                    {
                        "Cx": float(end_anchor_xy[0]),
                        "Cy": float(end_anchor_xy[1]),
                        "Cz": f_pnts_zip(group, float(end_anchor_xy[0]), float(end_anchor_xy[1])),
                        "type": "Anchor: Distal group anchor",
                        "grp": 0,
                    },
                ]
            )

            coords = pnts_bnd_ascending[["Cx", "Cy"]].to_numpy(dtype=float)
            dist_initial = np.linalg.norm(coords - start_anchor_xy, axis=1)
            dist_distal = np.linalg.norm(coords - end_anchor_xy, axis=1)

            right_idx_initial = -1
            right_idx_distal = -1
            left_idx_initial = -1
            left_idx_distal = -1

            path_vector = np.array([end_anchor_xy[0] - start_anchor_xy[0], end_anchor_xy[1] - start_anchor_xy[1], 0.0], dtype=float)
            for idx1 in np.argsort(dist_initial):
                idx2 = idx1 + 1
                if idx2 >= len(pnts_bnd_ascending):
                    idx2 = 0
                vec1 = np.array([
                    coords[idx1, 0] - start_anchor_xy[0],
                    coords[idx1, 1] - start_anchor_xy[1],
                    0.0,
                ])
                vec2 = np.array([
                    coords[idx2, 0] - start_anchor_xy[0],
                    coords[idx2, 1] - start_anchor_xy[1],
                    0.0,
                ])
                cp1 = np.cross(path_vector, vec1)[2]
                cp2 = np.cross(path_vector, vec2)[2]
                if cp1 > 0 and cp2 < 0:
                    left_idx_distal = int(idx1)
                    right_idx_initial = int(idx2)
                    break

            path_vector = np.array([start_anchor_xy[0] - end_anchor_xy[0], start_anchor_xy[1] - end_anchor_xy[1], 0.0], dtype=float)
            for idx1 in np.argsort(dist_distal):
                idx2 = idx1 + 1
                if idx2 >= len(pnts_bnd_ascending):
                    idx2 = 0
                vec1 = np.array([
                    coords[idx1, 0] - end_anchor_xy[0],
                    coords[idx1, 1] - end_anchor_xy[1],
                    0.0,
                ])
                vec2 = np.array([
                    coords[idx2, 0] - end_anchor_xy[0],
                    coords[idx2, 1] - end_anchor_xy[1],
                    0.0,
                ])
                cp1 = np.cross(path_vector, vec1)[2]
                cp2 = np.cross(path_vector, vec2)[2]
                if cp1 > 0 and cp2 < 0:
                    right_idx_distal = int(idx1)
                    left_idx_initial = int(idx2)
                    break

            bnd_sides_right = _boundary_arc_between(pnts_bnd_ascending, right_idx_initial, right_idx_distal)
            bnd_sides_left = _boundary_arc_between(pnts_bnd_ascending, left_idx_initial, left_idx_distal)
            if bnd_sides_right is not None and bnd_sides_left is not None:
                corner_check = _boundary_side_corners_check([bnd_sides_right, bnd_sides_left])
                if all(corner_check):
                    pnts0_initial = pnts_bnd_ascending.iloc[[left_idx_distal, right_idx_initial]][["Cx", "Cy"]].to_numpy(dtype=float)
                    pnts0_distal = pnts_bnd_ascending.iloc[[left_idx_initial, right_idx_distal]][["Cx", "Cy"]].to_numpy(dtype=float)
                    start_strip = _resolve_end_strip_constraint_onr(
                        group,
                        polygon,
                        start_anchor_xy,
                        anchor_mid_xy,
                        end_strips_initial,
                        True,
                        pnts0_initial,
                    )
                    end_strip = _resolve_end_strip_constraint_onr(
                        group,
                        polygon,
                        end_anchor_xy,
                        anchor_mid_xy,
                        end_strips_distal,
                        False,
                        pnts0_distal,
                    )
                    group_boundary_angles = np.array([
                        float(start_strip["angle"]) if start_strip is not None else np.nan,
                        float(end_strip["angle"]) if end_strip is not None else np.nan,
                    ], dtype=float)
                    return anchors, {
                        "polygon": polygon,
                        "bool_direction": bool_direction,
                        "ordered_groups": groups,
                        "end_strips_initial": end_strips_initial,
                        "end_strips_distal": end_strips_distal,
                        "start_strip": start_strip,
                        "end_strip": end_strip,
                        "group_boundary_angles": group_boundary_angles,
                        "list_bnd_sides_right": [bnd_sides_right],
                        "list_bnd_sides_left": [bnd_sides_left],
                    }

        pnts0_initial = pnts_bnd_ascending.iloc[[0, len(pnts_bnd_ascending) - 1]][["Cx", "Cy"]].to_numpy(dtype=float)
        pnts0_distal = pnts_bnd_ascending.iloc[[len(pnts_bnd_ascending) - 1, 0]][["Cx", "Cy"]].to_numpy(dtype=float)
        start_anchor_xy = _resolve_end_anchor_onr(group, polygon, end_anchors_initial, True, pnts0_initial)
        end_anchor_xy = _resolve_end_anchor_onr(group, polygon, end_anchors_distal, False, pnts0_distal)
        if start_anchor_xy is None or end_anchor_xy is None:
            return None
        anchor_mid_xy = (start_anchor_xy + end_anchor_xy) / 2.0
        anchors = pd.DataFrame(
            [
                {
                    "Cx": float(start_anchor_xy[0]),
                    "Cy": float(start_anchor_xy[1]),
                    "Cz": f_pnts_zip(group, float(start_anchor_xy[0]), float(start_anchor_xy[1])),
                    "type": "Anchor: Initial group anchor",
                    "grp": 0,
                },
                {
                    "Cx": float(anchor_mid_xy[0]),
                    "Cy": float(anchor_mid_xy[1]),
                    "Cz": f_pnts_zip(group, float(anchor_mid_xy[0]), float(anchor_mid_xy[1])),
                    "type": "Anchor: Group center",
                    "grp": 1,
                },
                {
                    "Cx": float(end_anchor_xy[0]),
                    "Cy": float(end_anchor_xy[1]),
                    "Cz": f_pnts_zip(group, float(end_anchor_xy[0]), float(end_anchor_xy[1])),
                    "type": "Anchor: Distal group anchor",
                    "grp": 0,
                },
            ]
        )
        start_strip = _resolve_end_strip_constraint_onr(group, polygon, start_anchor_xy, anchors.iloc[1][["Cx", "Cy"]].to_numpy(dtype=float), end_strips_initial, True, pnts0_initial)
        end_strip = _resolve_end_strip_constraint_onr(group, polygon, end_anchor_xy, anchors.iloc[1][["Cx", "Cy"]].to_numpy(dtype=float), end_strips_distal, False, pnts0_distal)
        return anchors, {
            "polygon": polygon,
            "bool_direction": bool_direction,
            "ordered_groups": groups,
            "end_strips_initial": end_strips_initial,
            "end_strips_distal": end_strips_distal,
            "start_strip": start_strip,
            "end_strip": end_strip,
            "list_bnd_sides_right": None,
            "list_bnd_sides_left": None,
        }

    launch_anchor = None
    list_bnd_sides_right: list[pd.DataFrame] = []
    list_bnd_sides_left: list[pd.DataFrame] = []
    center_lines: list[np.ndarray] = []
    grpbnd_agl_initial = None
    grpbnd_agl_distal = None

    for i, group in enumerate(groups, start=1):
        pnts_i_bnd = group[group["IDB"] != 0].copy()
        if pnts_i_bnd.empty:
            return None

        if i == 1 or i == count_pnts:
            bool_initial = i == 1
            end_anchor_cfg = end_anchors_initial if bool_initial else end_anchors_distal
            end_strip_cfg = end_strips_initial if bool_initial else end_strips_distal
            pnts_i_bnd_ascending = f_pnts_bnd_ascending(pnts_i_bnd)
            if pnts_i_bnd_ascending.empty:
                return None

            if bool_initial:
                pnts0 = pnts_i_bnd_ascending.iloc[[0, len(pnts_i_bnd_ascending) - 1]][["Cx", "Cy"]].to_numpy(dtype=float)
            else:
                pnts0 = pnts_i_bnd_ascending.iloc[[len(pnts_i_bnd_ascending) - 1, 0]][["Cx", "Cy"]].to_numpy(dtype=float)

            anchor_xy = _resolve_end_anchor_onr(group, polygon, end_anchor_cfg, bool_initial, pnts0)
            if anchor_xy is None:
                return None
            bnd_sides_right, bnd_sides_left = _split_end_boundary_sides(group, pnts_i_bnd, pnts0, end_anchor_cfg, bool_initial)
            if bnd_sides_right is None or bnd_sides_left is None:
                return None

            if bool_direction:
                end_strip_meta = _resolve_end_strip_constraint_onr(group, polygon, anchor_xy, pnts0.mean(axis=0), end_strip_cfg, bool_initial, pnts0)
                if end_strip_meta is not None:
                    if bool_initial:
                        grpbnd_agl_initial = float(end_strip_meta["angle"])
                    else:
                        grpbnd_agl_distal = float(end_strip_meta["angle"])

            if bool_initial:
                if int(pnts_i_bnd_ascending.iloc[0]["IDB"]) == int(bnd_sides_right.iloc[0]["IDB"]):
                    bnd_sides_right, bnd_sides_left = bnd_sides_left, bnd_sides_right
                launch_idx = int(round(len(pnts_i_bnd) / 2.0)) - 1
                launch_idx = min(max(launch_idx, 0), len(pnts_i_bnd) - 1)
                launch_anchor = pnts_i_bnd.iloc[[launch_idx]].copy()
            else:
                if int(pnts_i_bnd_ascending.iloc[0]["IDB"]) == int(bnd_sides_left.iloc[0]["IDB"]):
                    bnd_sides_right, bnd_sides_left = bnd_sides_left, bnd_sides_right

            list_pnts_i_bnd_sides = [bnd_sides_right, bnd_sides_left]
        else:
            list_pnts_i_bnd_sides = f_pnts_bnd_sides(pnts_i_bnd, idb_max)
            list_pnts_i_bnd_sides = _f_pnts_bnd_sides_switch2(list_pnts_i_bnd_sides, launch_anchor, idb_max)

        if len(list_pnts_i_bnd_sides) != 2:
            return None

        corner_check = _boundary_side_corners_check(list_pnts_i_bnd_sides)
        if i == 1:
            if not corner_check[1]:
                return None
        elif i == count_pnts:
            if not corner_check[0]:
                return None
        else:
            if not all(corner_check):
                return None

        center_ls = f_pnts_bnd_sides_centerLs(list_pnts_i_bnd_sides)
        if center_ls is None:
            return None

        center1 = np.array(
            [
                center_ls[0, 0],
                center_ls[0, 1],
                f_pnts_zip(group, float(center_ls[0, 0]), float(center_ls[0, 1])),
            ],
            dtype=float,
        )
        center2 = np.array(
            [
                center_ls[1, 0],
                center_ls[1, 1],
                f_pnts_zip(group, float(center_ls[1, 0]), float(center_ls[1, 1])),
            ],
            dtype=float,
        )
        center_lines.extend([center1, center2])
        list_bnd_sides_right.append(list_pnts_i_bnd_sides[0])
        list_bnd_sides_left.append(list_pnts_i_bnd_sides[1])

    anchor_rows: list[dict[str, Any]] = []
    anchor_rows.append(
        {
            "Cx": float(center_lines[0][0]),
            "Cy": float(center_lines[0][1]),
            "Cz": float(center_lines[0][2]),
            "type": "Anchor: Initial group anchor",
            "grp": 0,
        }
    )

    for idx in range(len(center_lines) - 1):
        anchor_xyz = (center_lines[idx] + center_lines[idx + 1]) / 2.0
        if idx % 2 == 0:
            anchor_type = "Anchor: Group center"
            grp = idx // 2 + 1
        else:
            anchor_type = "Anchor: Inter-group center"
            grp = 0
        anchor_rows.append(
            {
                "Cx": float(anchor_xyz[0]),
                "Cy": float(anchor_xyz[1]),
                "Cz": float(anchor_xyz[2]),
                "type": anchor_type,
                "grp": grp,
            }
        )

    anchor_rows.append(
        {
            "Cx": float(center_lines[-1][0]),
            "Cy": float(center_lines[-1][1]),
            "Cz": float(center_lines[-1][2]),
            "type": "Anchor: Distal group anchor",
            "grp": 0,
        }
    )

    anchors = pd.DataFrame(anchor_rows)
    count_anchor = len(anchors)

    if _is_scalar_r_config(end_anchors_initial):
        anchor_seed_initial = anchors.iloc[[1, 0]][["Cx", "Cy"]].to_numpy(dtype=float)
        anchor_xy_initial = _nearest_boundary_intersection_along_ray(
            polygon,
            np.array(anchor_seed_initial[0], dtype=float),
            np.array(anchor_seed_initial[1] - anchor_seed_initial[0], dtype=float),
        )
    else:
        anchor_xy_initial = np.array(end_anchors_initial[:2], dtype=float)

    if _is_scalar_r_config(end_anchors_distal):
        anchor_seed_distal = anchors.iloc[[count_anchor - 2, count_anchor - 1]][["Cx", "Cy"]].to_numpy(dtype=float)
        anchor_xy_distal = _nearest_boundary_intersection_along_ray(
            polygon,
            np.array(anchor_seed_distal[0], dtype=float),
            np.array(anchor_seed_distal[1] - anchor_seed_distal[0], dtype=float),
        )
    else:
        anchor_xy_distal = np.array(end_anchors_distal[:2], dtype=float)

    if anchor_xy_initial is None or anchor_xy_distal is None:
        return None

    anchors.loc[0, ["Cx", "Cy", "Cz"]] = [
        float(anchor_xy_initial[0]),
        float(anchor_xy_initial[1]),
        f_pnts_zip(groups[0], float(anchor_xy_initial[0]), float(anchor_xy_initial[1])),
    ]
    anchors.loc[count_anchor - 1, ["Cx", "Cy", "Cz"]] = [
        float(anchor_xy_distal[0]),
        float(anchor_xy_distal[1]),
        f_pnts_zip(groups[-1], float(anchor_xy_distal[0]), float(anchor_xy_distal[1])),
    ]

    start_reference = anchors.iloc[min(2, len(anchors) - 1)][["Cx", "Cy"]].to_numpy(dtype=float)
    end_reference = anchors.iloc[max(len(anchors) - 3, 0)][["Cx", "Cy"]].to_numpy(dtype=float)
    pnts0_initial = groups[0][groups[0]["IDB"] != 0].sort_values("IDB").iloc[[0, -1]][["Cx", "Cy"]].to_numpy(dtype=float)
    pnts0_distal_sorted = groups[-1][groups[-1]["IDB"] != 0].sort_values("IDB")
    pnts0_distal = pnts0_distal_sorted.iloc[[-1, 0]][["Cx", "Cy"]].to_numpy(dtype=float)
    start_strip = _resolve_end_strip_constraint_onr(groups[0], polygon, anchor_xy_initial, start_reference, end_strips_initial, True, pnts0_initial)
    end_strip = _resolve_end_strip_constraint_onr(groups[-1], polygon, anchor_xy_distal, end_reference, end_strips_distal, False, pnts0_distal)
    return anchors, {
        "polygon": polygon,
        "bool_direction": bool_direction,
        "ordered_groups": groups,
        "end_strips_initial": end_strips_initial,
        "end_strips_distal": end_strips_distal,
        "start_strip": start_strip,
        "end_strip": end_strip,
        "group_boundary_angles": np.array([grpbnd_agl_initial, grpbnd_agl_distal], dtype=float),
        "list_bnd_sides_right": list_bnd_sides_right,
        "list_bnd_sides_left": list_bnd_sides_left,
    }


def f_strips_package(
    list_pnts: list[pd.DataFrame],
    polygon: Polygon,
    anchors: pd.DataFrame,
    path_meta: dict[str, Any],
    end_strips_initial: Any,
    end_strips_distal: Any,
    min_strip_horizontal_length: float,
    src: rasterio.io.DatasetReader,
) -> dict[str, Any]:
    path_line = LineString(anchors[["Cx", "Cy"]].to_numpy(dtype=float))
    if path_line.length == 0:
        return {
            "strip_gdf": None,
            "strip_df": None,
            "strip_nodes": None,
            "triangle_gdf": None,
            "triangle_df": None,
            "paras_strip": None,
        }

    idb_max = int(max(group["IDB"].max() for group in list_pnts if not group.empty))

    start_constraint = path_meta.get("start_strip") if path_meta else None
    end_constraint = path_meta.get("end_strip") if path_meta else None
    start_cut = _path_cut_from_constraint(path_line, polygon, start_constraint)
    end_cut = _path_cut_from_constraint(path_line, polygon, end_constraint)

    station_info = _build_station_layout_onr(
        path_line,
        anchors,
        list_pnts,
        path_meta,
        idb_max,
        start_cut,
        end_cut,
        min_strip_horizontal_length,
    )
    if station_info is None:
        return {
            "strip_gdf": None,
            "strip_df": None,
            "strip_nodes": None,
            "triangle_gdf": None,
            "triangle_df": None,
            "paras_strip": None,
        }

    cutlines = _build_station_cutlines(polygon, path_line, station_info)
    cutlines = _filter_cutlines_inside_polygon(polygon, cutlines)
    cutlines = _remove_intersecting_cutlines(cutlines)
    split_package = _split_polygon_into_strips(polygon, path_line, anchors, cutlines)
    if split_package is None:
        return {
            "strip_gdf": None,
            "strip_df": None,
            "strip_nodes": None,
            "triangle_gdf": None,
            "triangle_df": None,
            "paras_strip": None,
        }

    strip_polygons, node_events = split_package
    triangle_package = _build_triangle_surface_package(strip_polygons, path_line, node_events, anchors, src)
    strip_nodes = triangle_package["strip_nodes"]
    strip_df = triangle_package["strip_df"]
    triangle_df = triangle_package["triangle_df"]
    paras_strip = triangle_package["paras_strip"]

    strip_gdf = None
    if strip_polygons:
        strip_gdf = gpd.GeoDataFrame(strip_df.copy(), geometry=strip_polygons, crs=src.crs)

    triangle_gdf = None
    if triangle_package["triangle_geometries"]:
        triangle_gdf = gpd.GeoDataFrame(triangle_df.copy(), geometry=triangle_package["triangle_geometries"], crs=src.crs)

    return {
        "strip_gdf": strip_gdf,
        "strip_df": strip_df,
        "strip_nodes": strip_nodes,
        "triangle_gdf": triangle_gdf,
        "triangle_df": triangle_df,
        "paras_strip": paras_strip,
    }


def _group_boundary_events_onr(
    list_pnts: list[pd.DataFrame],
    anchors: pd.DataFrame,
    path_line: LineString,
    path_meta: dict[str, Any] | None,
    idb_max: int,
) -> list[dict[str, Any]]:
    sides_right = path_meta.get("list_bnd_sides_right") if path_meta else None
    sides_left = path_meta.get("list_bnd_sides_left") if path_meta else None
    if not sides_right or not sides_left or len(sides_right) != len(list_pnts) or len(sides_left) != len(list_pnts):
        return _build_group_boundary_events(list_pnts, anchors, path_line, idb_max)

    anchor_coords = anchors[["Cx", "Cy"]].to_numpy(dtype=float)
    anchor_grps = anchors["grp"].to_numpy(dtype=int)
    events: list[dict[str, Any]] = []
    count_grp = len(list_pnts)

    for idx, (side_right, side_left) in enumerate(zip(sides_right, sides_left), start=1):
        if side_right is None or side_left is None or side_right.empty or side_left.empty:
            continue
        matched = np.where(anchor_grps == idx)[0]
        if matched.size:
            center_idx = int(matched[0])
        else:
            center_idx = min(max(idx, 1), len(anchor_coords) - 2)
        prev_idx = max(center_idx - 1, 0)
        next_idx = min(center_idx + 1, len(anchor_coords) - 1)
        direction = anchor_coords[next_idx] - anchor_coords[prev_idx]
        norm = float(np.linalg.norm(direction))
        if norm == 0:
            direction = anchor_coords[-1] - anchor_coords[0]
            norm = float(np.linalg.norm(direction))
        if norm == 0:
            direction = np.array([1.0, 0.0], dtype=float)
        else:
            direction = direction / norm

        right_coords = side_right[["Cx", "Cy"]].to_numpy(dtype=float)
        left_coords = side_left[["Cx", "Cy"]].to_numpy(dtype=float)
        if len(right_coords) < 2 or len(left_coords) < 2:
            continue
        corners = np.vstack([
            right_coords[0],
            right_coords[-1],
            left_coords[0],
            left_coords[-1],
        ])
        segment_map = {
            "upper": np.vstack([corners[3], corners[0]]),
            "lower": np.vstack([corners[2], corners[1]]),
        }

        for key, label in (("upper", "Intersection: Group boundary upper"), ("lower", "Intersection: Group boundary lower")):
            if key == "upper" and idx == 1:
                continue
            if key == "lower" and idx == count_grp:
                continue
            segment = segment_map[key]
            intersection = _polyline_segment_intersection(path_line, segment)
            if intersection is None:
                continue
            events.append(
                {
                    "Cx": float(intersection[0]),
                    "Cy": float(intersection[1]),
                    "Czp": 0.0,
                    "type": label,
                    "grp": idx,
                    "station": 0,
                    "chainage": float(path_line.project(Point(float(intersection[0]), float(intersection[1])))),
                    "angle": _segment_angle(segment),
                }
            )

    events.sort(key=lambda item: item["chainage"])
    return events


def _station_chainages_onr(
    path_line: LineString,
    start_cut: dict[str, Any] | None,
    end_cut: dict[str, Any] | None,
    min_strip_horizontal_length: float,
) -> list[float] | None:
    start_chainage = float(start_cut["chainage"]) if start_cut is not None else 0.0
    end_chainage = float(end_cut["chainage"]) if end_cut is not None else float(path_line.length)
    if end_chainage <= start_chainage:
        return None

    length_path = end_chainage - start_chainage
    if not math.isfinite(length_path) or not math.isfinite(min_strip_horizontal_length):
        return None

    count_strip = 0
    if min_strip_horizontal_length <= 0:
        count_strip = 2
    else:
        while length_path / (count_strip + 1) >= min_strip_horizontal_length:
            count_strip += 1
        count_strip = max(count_strip, 2)

    inner = np.linspace(start_chainage, end_chainage, count_strip + 1)[1:-1].tolist()
    chainages = inner
    if start_cut is not None:
        chainages = [start_chainage] + chainages
    if end_cut is not None:
        chainages = chainages + [end_chainage]
    return [float(item) for item in chainages if 0.0 < float(item) < float(path_line.length)]


def _interpolate_station_angle_onr(
    chainage: float,
    anchors: pd.DataFrame,
    boundary_events: list[dict[str, Any]],
    end_angle_initial: float,
    end_angle_distal: float,
    start_cut: dict[str, Any] | None,
    end_cut: dict[str, Any] | None,
) -> float:
    start_chainage = float(start_cut["chainage"]) if start_cut is not None else 0.0
    end_chainage = float(end_cut["chainage"]) if end_cut is not None else float(anchors.iloc[-1]["Lhrz"])

    if len(boundary_events) >= 2:
        upper_idx = -1
        lower_idx = len(boundary_events)
        for idx, event in enumerate(boundary_events):
            if float(event["chainage"]) > chainage:
                lower_idx = idx
                upper_idx = idx - 1
                break

        if lower_idx == len(boundary_events):
            lhrz_diff = end_chainage - float(boundary_events[upper_idx]["chainage"])
            lhrz_diff_station = chainage - float(boundary_events[upper_idx]["chainage"])
            angle_upper = float(boundary_events[upper_idx]["angle"])
            angle_lower = end_angle_distal
        elif upper_idx < 0:
            lhrz_diff = float(boundary_events[lower_idx]["chainage"]) - start_chainage
            lhrz_diff_station = chainage - start_chainage
            angle_upper = end_angle_initial
            angle_lower = float(boundary_events[lower_idx]["angle"])
        else:
            lhrz_diff = float(boundary_events[lower_idx]["chainage"]) - float(boundary_events[upper_idx]["chainage"])
            lhrz_diff_station = chainage - float(boundary_events[upper_idx]["chainage"])
            angle_upper = float(boundary_events[upper_idx]["angle"])
            angle_lower = float(boundary_events[lower_idx]["angle"])

        if start_cut is not None:
            if upper_idx < 0 or (0 <= upper_idx < len(boundary_events) and float(boundary_events[upper_idx]["chainage"]) <= start_chainage):
                angle_upper = float(start_cut["angle"])
                lhrz_diff = float(boundary_events[lower_idx]["chainage"]) - start_chainage if lower_idx < len(boundary_events) else end_chainage - start_chainage
                lhrz_diff_station = chainage - start_chainage

        if end_cut is not None:
            if lower_idx == len(boundary_events) or (0 <= lower_idx < len(boundary_events) and float(boundary_events[lower_idx]["chainage"]) >= end_chainage):
                angle_lower = float(end_cut["angle"])
                lhrz_diff = end_chainage - (float(boundary_events[upper_idx]["chainage"]) if upper_idx >= 0 else start_chainage)
                lhrz_diff_station = chainage - (float(boundary_events[upper_idx]["chainage"]) if upper_idx >= 0 else start_chainage)

        if lhrz_diff <= 0:
            return _normalize_angle(angle_upper)
        angle_diff = _normalize_angle(angle_lower - angle_upper)
        return _normalize_angle(angle_upper + angle_diff / lhrz_diff * lhrz_diff_station)

    one_group_angle = _normalize_angle(
        math.atan2(float(anchors.iloc[-1]["Cy"] - anchors.iloc[0]["Cy"]), float(anchors.iloc[-1]["Cx"] - anchors.iloc[0]["Cx"])) - math.pi / 2.0
    )
    pivot_chainage = float(anchors.iloc[min(1, len(anchors) - 1)]["Lhrz"])
    if chainage <= pivot_chainage:
        angle_upper = float(start_cut["angle"]) if start_cut is not None else end_angle_initial
        angle_lower = one_group_angle
        base_chainage = start_chainage
        limit_chainage = pivot_chainage
    else:
        angle_upper = one_group_angle
        angle_lower = float(end_cut["angle"]) if end_cut is not None else end_angle_distal
        base_chainage = pivot_chainage
        limit_chainage = end_chainage
    lhrz_diff = max(limit_chainage - base_chainage, 0.0)
    if lhrz_diff == 0:
        return _normalize_angle(angle_upper)
    angle_diff = _normalize_angle(angle_lower - angle_upper)
    return _normalize_angle(angle_upper + angle_diff / lhrz_diff * (chainage - base_chainage))


def _build_station_layout_onr(
    path_line: LineString,
    anchors: pd.DataFrame,
    list_pnts: list[pd.DataFrame],
    path_meta: dict[str, Any] | None,
    idb_max: int,
    start_cut: dict[str, Any] | None,
    end_cut: dict[str, Any] | None,
    min_strip_horizontal_length: float,
) -> dict[str, Any] | None:
    anchors_work = anchors.copy()
    anchor_chainages = [float(path_line.project(Point(float(row["Cx"]), float(row["Cy"])))) for _, row in anchors_work.iterrows()]
    anchors_work.loc[:, "Lhrz"] = anchor_chainages

    boundary_events = _group_boundary_events_onr(list_pnts, anchors_work, path_line, path_meta, idb_max)
    station_chainages = _station_chainages_onr(path_line, start_cut, end_cut, min_strip_horizontal_length)
    if station_chainages is None:
        return None

    end_angles = path_meta.get("group_boundary_angles") if path_meta else None
    if end_angles is not None and len(end_angles) == 2 and np.all(np.isfinite(np.asarray(end_angles, dtype=float))):
        end_angle_initial = float(end_angles[0])
        end_angle_distal = float(end_angles[1])
    else:
        end_angle_initial = float(start_cut["angle"]) if start_cut is not None else _path_normal_angle(path_line, 0.0)
        end_angle_distal = float(end_cut["angle"]) if end_cut is not None else _path_normal_angle(path_line, path_line.length)

    station_events: list[dict[str, Any]] = []
    for idx, chainage in enumerate(station_chainages, start=1):
        point = path_line.interpolate(chainage)
        station_events.append(
            {
                "Cx": float(point.x),
                "Cy": float(point.y),
                "Czp": 0.0,
                "type": "Station: Strip boundary",
                "grp": 0,
                "station": idx,
                "chainage": float(chainage),
                "angle": float(
                    _interpolate_station_angle_onr(
                        float(chainage),
                        anchors_work,
                        boundary_events,
                        end_angle_initial,
                        end_angle_distal,
                        start_cut,
                        end_cut,
                    )
                ),
            }
        )

    anchor_events = []
    for _, row in anchors.iterrows():
        point = Point(float(row["Cx"]), float(row["Cy"]))
        anchor_events.append(
            {
                "Cx": float(row["Cx"]),
                "Cy": float(row["Cy"]),
                "Czp": float(row["Cz"]),
                "type": row["type"],
                "grp": int(row["grp"]),
                "station": 0,
                "chainage": float(path_line.project(point)),
            }
        )

    return {
        "anchor_events": anchor_events,
        "boundary_events": boundary_events,
        "station_events": station_events,
        "start_cut": start_cut,
        "end_cut": end_cut,
    }


def f_pnts_MinDist_intercept(pnts: pd.DataFrame, dip: np.ndarray, intercept: float) -> pd.DataFrame:
    coords = pnts[["Cx", "Cy", "Cz"]].to_numpy(dtype=float)
    distances = np.abs(coords @ dip + intercept)
    return pnts.loc[[pnts.index[int(np.argmin(distances))]]].copy()


def f_pnts_zip(pnts: pd.DataFrame, coord_x: float, coord_y: float) -> float:
    row = pnts.iloc[0]
    return float(row["_c"] + row["_a"] * coord_x + row["_b"] * coord_y)


def f_pnts_distances(column_anchors: pd.DataFrame) -> np.ndarray:
    values = column_anchors.to_numpy(dtype=float)
    distances: list[float] = []
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            distances.append(float(np.linalg.norm(values[i] - values[j])))
    return np.array(distances, dtype=float)


def f_list2pnts(list_pnts: list[pd.DataFrame]) -> pd.DataFrame:
    if not list_pnts:
        return pd.DataFrame(columns=PNT_COLUMNS)
    return pd.concat(list_pnts, ignore_index=True)


def f_pnts2paras_path(pnts: pd.DataFrame, paras_strip: tuple[float, float, float, float] | None) -> pd.DataFrame:
    coords = pnts[["Cx", "Cy", "Cz"]].to_numpy(dtype=float)
    if len(coords) < 2:
        data = {"Dall": 0.0, "Lall": 0.0, "Dhrz": 0.0, "Lhrz": 0.0, "Dlng": 0.0, "Hfnl": 0.0, "Hmax": 0.0, "Hint": 0.0, "Gmax": 0.0, "Gint": 0.0, "Tall": 0.0, "Thrz": 0.0, "Tlng": 0.0, "miuapp": 0.0}
        return pd.DataFrame([data])

    horizontal_steps = np.linalg.norm(np.diff(coords[:, :2], axis=0), axis=1)
    vertical_steps = np.diff(coords[:, 2])
    spatial_steps = np.linalg.norm(np.diff(coords, axis=0), axis=1)

    l_path = float(horizontal_steps.sum())
    h_path = float(vertical_steps.sum())
    s_path = float(spatial_steps.sum())

    l_series = np.concatenate([[0.0], np.cumsum(horizontal_steps)])
    h_series = np.concatenate([[0.0], np.cumsum(vertical_steps)])
    miuapp_path = float(-h_path / l_path) if l_path else 0.0
    g_series = (-h_series) - miuapp_path * l_series

    h_integral = 0.0
    g_integral = 0.0
    for i in range(1, len(l_series)):
        l = l_series[i] - l_series[i - 1]
        h_integral += l * ((-h_series[i]) + (-h_series[i - 1])) / 2.0
        g_integral += l * (g_series[i] + g_series[i - 1]) / 2.0

    d_overall = float(np.linalg.norm(coords[-1] - coords[0]))
    d_horizontal = float(np.linalg.norm(coords[-1, :2] - coords[0, :2]))
    d_longitudinal = float(math.sqrt(l_path**2 + h_path**2))

    data = {
        "Dall": d_overall,
        "Lall": s_path,
        "Dhrz": d_horizontal,
        "Lhrz": l_path,
        "Dlng": d_longitudinal,
        "Hfnl": -h_path,
        "Hmax": float(-h_series.min()),
        "Hint": h_integral,
        "Gmax": float(g_series.max()),
        "Gint": g_integral,
        "Tall": float(s_path / d_overall) if d_overall else 0.0,
        "Thrz": float(l_path / d_horizontal) if d_horizontal else 0.0,
        "Tlng": float(s_path / d_longitudinal) if d_longitudinal else 0.0,
        "miuapp": miuapp_path,
    }
    if paras_strip is not None:
        l_hrz2, a_hrz, l_all2, a_all = paras_strip
        data.update(
            {
                "Aall": a_all,
                "Wall": float(a_all / l_all2) if l_all2 else 0.0,
                "epsall": float((l_all2**2) / a_all) if a_all else 0.0,
                "Ahrz": a_hrz,
                "Whrz": float(a_hrz / l_hrz2) if l_hrz2 else 0.0,
                "epshrz": float((l_hrz2**2) / a_hrz) if a_hrz else 0.0,
            }
        )
    return pd.DataFrame([data])


def f_pnts_save_path(
    pnts_path: pd.DataFrame,
    paras_strip: tuple[float, float, float, float] | None,
    name: str,
    min_group_point_count: int,
    min_group_anchor_distance: float,
    min_strip_horizontal_length: float,
    min_end_ratio_initial: float,
    min_end_ratio_distal: float,
    end_anchors_initial: Any,
    end_anchors_distal: Any,
    end_strips_initial: Any,
    end_strips_distal: Any,
    crs_out: Any,
    file_out: Path,
    line_id: int = 0,
) -> gpd.GeoDataFrame:
    paras = f_pnts2paras_path(pnts_path, paras_strip)
    input_paras = {
        "ID": line_id,
        "Name": name,
        "iMGPC": min_group_point_count,
        "iMGAD": min_group_anchor_distance,
        "iMSHL": min_strip_horizontal_length,
        "iMERI": min_end_ratio_initial,
        "iMERD": min_end_ratio_distal,
        "iEAI": _config_to_string(end_anchors_initial),
        "iEAD": _config_to_string(end_anchors_distal),
        "iESI": _config_to_string(end_strips_initial),
        "iESD": _config_to_string(end_strips_distal),
    }
    data_df = pd.concat([pd.DataFrame([input_paras]), paras], axis=1)
    coords = pnts_path[["Cx", "Cy"]].to_numpy(dtype=float)
    line = LineString(coords)
    if coords.shape[0] < 2 or line.is_empty or float(line.length) <= 0:
        return gpd.GeoDataFrame(columns=[*data_df.columns, "geometry"], geometry="geometry", crs=crs_out)
    geometry = [line]
    gdf = gpd.GeoDataFrame(data_df, geometry=geometry, crs=crs_out)
    _safe_write_gdf(gdf, file_out)
    _write_excel(pd.concat([data_df, pnts_path.reset_index(drop=True)], axis=1), file_out.with_suffix(".xlsx"))
    return gdf


def _has_valid_path_gdf(path_gdf: gpd.GeoDataFrame | None) -> bool:
    if path_gdf is None or path_gdf.empty or "geometry" not in path_gdf:
        return False
    geometries = path_gdf.geometry.dropna()
    if geometries.empty:
        return False
    return any((not geometry.is_empty) and float(geometry.length) > 0 for geometry in geometries)


def _load_parameters(parameters: str | Path | int, polygon_count: int) -> pd.DataFrame | None:
    if parameters == 0:
        return None
    df = pd.read_excel(parameters)
    df = df.reset_index(drop=True)
    if len(df) < polygon_count:
        warnings.warn("Parameter rows are fewer than polygon count; defaults will be used for remaining polygons.")
    return df


def _parameters_for_index(parameters: pd.DataFrame | None, index: int) -> PolygonParameters:
    if parameters is None or index >= len(parameters):
        return PolygonParameters(name=f"lsd{index + 1:06d}")
    row = parameters.iloc[index]
    return PolygonParameters(
        name=str(row.iloc[0]),
        min_group_point_count=int(row.iloc[1]),
        min_group_anchor_distance=float(row.iloc[2]),
        min_strip_horizontal_length=float(row.iloc[3]),
        min_end_ratio_initial=float(row.iloc[4]),
        min_end_ratio_distal=float(row.iloc[5]),
    )


def _load_config_features(config: Any) -> list[np.ndarray] | None:
    if isinstance(config, (str, Path)):
        gdf = gpd.read_file(config)
        features: list[np.ndarray] = []
        for geometry in gdf.geometry:
            if geometry is None or geometry.is_empty:
                features.append(np.array([]))
            elif geometry.geom_type == "Point":
                features.append(np.array([geometry.x, geometry.y], dtype=float))
            elif geometry.geom_type == "LineString":
                coords = list(geometry.coords)
                if len(coords) >= 2:
                    direction = math.atan2(coords[-1][1] - coords[0][1], coords[-1][0] - coords[0][0])
                    if direction < 0:
                        direction += math.pi * 2.0
                    features.append(np.array([coords[0][0], coords[0][1], coords[-1][0], coords[-1][1], direction], dtype=float))
                else:
                    features.append(np.array([]))
            else:
                centroid = geometry.centroid
                features.append(np.array([centroid.x, centroid.y], dtype=float))
        return features
    return None


def _config_for_index(raw_config: Any, loaded: list[np.ndarray] | None, index: int) -> Any:
    if loaded is not None:
        if index >= len(loaded):
            raise IndexError("Configuration feature count is smaller than polygon count.")
        return loaded[index]
    if isinstance(raw_config, (int, np.integer, float, np.floating)):
        return -abs(int(raw_config))
    return raw_config


def _is_scalar_r_config(value: Any) -> bool:
    return isinstance(value, (int, np.integer, float, np.floating))


def _angle_included(angle1: float, angle2: float) -> float:
    angle_included = angle2 - angle1
    while angle_included > math.pi:
        angle_included -= math.pi * 2.0
    while angle_included < -math.pi:
        angle_included += math.pi * 2.0
    return abs(angle_included)


def _orient_end_strip_angle(angle: float, pnts0: np.ndarray) -> float:
    angle = _normalize_angle(float(angle))
    agl0 = math.atan2(float(pnts0[1, 1] - pnts0[0, 1]), float(pnts0[1, 0] - pnts0[0, 0]))
    if _angle_included(agl0, angle) > math.pi / 2.0:
        angle = _normalize_angle(angle - math.pi)
    return angle


def _resolve_end_anchor_onr(group: pd.DataFrame, polygon: Polygon, config: Any, initial: bool, pnts0: np.ndarray) -> np.ndarray | None:
    if isinstance(config, np.ndarray) and config.size >= 2:
        return config[:2].astype(float)
    if isinstance(config, (list, tuple, np.ndarray)) and len(config) == 2:
        return np.array(config, dtype=float)

    code = int(config)
    bnd = group[group["IDB"] != 0].copy()
    if bnd.empty:
        return group[["Cx", "Cy"]].mean().to_numpy(dtype=float)

    if code == -1:
        _, _, _, anchor = _pnts_bnd_end_mbb(bnd, pnts0)
        return None if anchor is None else np.array(anchor, dtype=float)
    if code == -2:
        _, _, _, anchor = _pnts_bnd_end_quad(bnd, pnts0)
        return None if anchor is None else np.array(anchor, dtype=float)
    if code == -3:
        _, _, anchor = _pnts_bnd_end_even(group, pnts0, initial)
        return None if anchor is None else np.array(anchor, dtype=float)
    return None


def _split_end_boundary_sides(
    group: pd.DataFrame,
    pnts_bnd: pd.DataFrame,
    pnts0: np.ndarray,
    config: Any,
    initial: bool,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    if isinstance(config, np.ndarray) and config.size >= 2:
        return _pnts_bnd_end_coords(pnts_bnd, pnts0, config[:2])
    if isinstance(config, (list, tuple, np.ndarray)) and len(config) == 2:
        return _pnts_bnd_end_coords(pnts_bnd, pnts0, np.array(config, dtype=float))

    code = int(config)
    if code == -1:
        right, left, _, _ = _pnts_bnd_end_mbb(pnts_bnd, pnts0)
        return right, left
    if code == -2:
        right, left, _, _ = _pnts_bnd_end_quad(pnts_bnd, pnts0)
        return right, left
    if code == -3:
        right, left, _ = _pnts_bnd_end_even(group, pnts0, initial)
        return right, left
    return None, None


def _f_pnts_bnd_sides_switch2(list_bnd_sides: list[pd.DataFrame], launch_anchor: pd.DataFrame | None, idb_max: int) -> list[pd.DataFrame]:
    if launch_anchor is None or len(list_bnd_sides) != 2:
        return list_bnd_sides
    bnd_sides_right = list_bnd_sides[0]
    bnd_sides_left = list_bnd_sides[1]
    idb_start_right = int(bnd_sides_right.iloc[0]["IDB"])
    idb_start_left = int(bnd_sides_left.iloc[0]["IDB"])
    idb_launch_anchor = int(launch_anchor.iloc[0]["IDB"])

    idb_diff_right = idb_start_right - idb_launch_anchor
    idb_diff_left = idb_start_left - idb_launch_anchor
    if idb_diff_right < 0:
        idb_diff_right += idb_max
    if idb_diff_left < 0:
        idb_diff_left += idb_max

    if idb_diff_right <= idb_diff_left:
        return list_bnd_sides
    return [bnd_sides_left, bnd_sides_right]


def _pnts_path_extending(pnts: np.ndarray, polygon: Polygon) -> np.ndarray | None:
    return _extend_path_anchor_to_boundary(pnts, polygon)


def _extend_path_anchor_to_boundary(pnts: np.ndarray, polygon: Polygon) -> np.ndarray | None:
    pnts = np.asarray(pnts, dtype=float)
    if pnts.shape[0] < 2 or pnts.shape[1] < 2:
        return None
    pnt1 = np.array(pnts[0, :2], dtype=float)
    pnt2 = np.array(pnts[1, :2], dtype=float)
    if np.allclose(pnt1, pnt2, atol=1e-12):
        return None

    direction = pnt2 - pnt1
    return _nearest_boundary_intersection_along_ray(polygon, pnt1, direction)


def _resolve_end_strip_constraint_onr(
    group: pd.DataFrame,
    polygon: Polygon,
    anchor_xy: np.ndarray,
    reference_xy: np.ndarray,
    config: Any,
    initial: bool,
    pnts0: np.ndarray,
) -> dict[str, Any] | None:
    if isinstance(config, np.ndarray) and config.size >= 5:
        segment = config[:4].astype(float).reshape(2, 2)
        angle = _orient_end_strip_angle(float(config[4]), pnts0)
        return {"mode": "segment", "segment": segment, "angle": angle}
    if isinstance(config, (list, tuple, np.ndarray)) and len(config) == 5:
        array = np.array(config, dtype=float)
        segment = array[:4].reshape(2, 2)
        angle = _orient_end_strip_angle(float(array[4]), pnts0)
        return {"mode": "segment", "segment": segment, "angle": angle}

    code = int(config)
    bnd = group[group["IDB"] != 0].copy()
    if bnd.empty:
        return None

    if code == -1:
        _, _, angle, _ = _pnts_bnd_end_mbb(bnd, pnts0)
    elif code == -2:
        _, _, angle, _ = _pnts_bnd_end_quad(bnd, pnts0)
    elif code == -3:
        angle = math.atan2(float(pnts0[1, 1] - pnts0[0, 1]), float(pnts0[1, 0] - pnts0[0, 0]))
    else:
        angle = float(code)

    angle = _orient_end_strip_angle(angle, pnts0)
    return {
        "mode": "angle",
        "angle": angle,
        "anchor": np.array(anchor_xy, dtype=float),
    }


def _pnts_bnd_end_coords(
    pnts_bnd: pd.DataFrame,
    pnts0: np.ndarray,
    coords: np.ndarray,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    pnt0 = pd.DataFrame([[float((pnts0[0, 0] + pnts0[1, 0]) / 2.0), float((pnts0[0, 1] + pnts0[1, 1]) / 2.0)]], columns=["Cx", "Cy"])
    end_anchor0 = pd.DataFrame([[float(coords[0]), float(coords[1])]], columns=["Cx", "Cy"])
    split = _pnts_group_end_split(pnts_bnd, pnt0, end_anchor0)
    if split is None:
        return None, None
    return split[0], split[1]


def _pnts_bnd_end_mbb(
    pnts_bnd: pd.DataFrame,
    pnts0: np.ndarray,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None, float, np.ndarray | None]:
    pnts_bnd = f_pnts_bnd_ascending(pnts_bnd)
    pnt0 = pd.DataFrame([[float((pnts0[0, 0] + pnts0[1, 0]) / 2.0), float((pnts0[0, 1] + pnts0[1, 1]) / 2.0)]], columns=["Cx", "Cy"])
    if len(pnts_bnd) == 3:
        angle = _normalize_angle(math.atan2(float(pnts0[1, 1] - pnts0[0, 1]), float(pnts0[1, 0] - pnts0[0, 0])))
        end_anchor0 = pnts_bnd.iloc[1][["Cx", "Cy"]].to_numpy(dtype=float)
        split = _pnts_group_end_split(pnts_bnd, pnt0, pd.DataFrame([end_anchor0], columns=["Cx", "Cy"]))
        if split is None:
            return None, None, angle, end_anchor0
        return split[0], split[1], angle, end_anchor0

    hull = Polygon(pnts_bnd[["Cx", "Cy"]].to_numpy(dtype=float)).convex_hull
    rect = hull.minimum_rotated_rectangle
    rect_coords = np.array(rect.exterior.coords[:-1], dtype=float)
    agl0 = math.atan2(float(pnts0[1, 1] - pnts0[0, 1]), float(pnts0[1, 0] - pnts0[0, 0]))
    pnt0_xy = pnt0.iloc[0][["Cx", "Cy"]].to_numpy(dtype=float)
    edge_angles: list[float] = []
    edge_angle_diff: list[float] = []
    edge_midpoints: list[np.ndarray] = []
    for i in range(4):
        p0 = rect_coords[i]
        p1 = rect_coords[(i + 1) % 4]
        angle = math.atan2(float(p1[1] - p0[1]), float(p1[0] - p0[0]))
        angle_diff = _angle_included(agl0, angle)
        if angle_diff > math.pi / 2.0:
            angle -= math.pi
            angle_diff = math.pi - angle_diff
        edge_angles.append(_normalize_angle(angle))
        edge_angle_diff.append(angle_diff)
        edge_midpoints.append((p0 + p1) / 2.0)

    min_diff = min(edge_angle_diff)
    candidate_idxs = [idx for idx, value in enumerate(edge_angle_diff) if math.isclose(value, min_diff, rel_tol=1e-9, abs_tol=1e-9)]
    if len(candidate_idxs) == 1:
        idx1 = candidate_idxs[0]
        rest = [idx for idx in range(4) if idx != idx1]
        idx2 = min(rest, key=lambda idx: edge_angle_diff[idx])
    else:
        idx1 = candidate_idxs[0]
        idx2 = candidate_idxs[1]
    dist1 = float(np.linalg.norm(edge_midpoints[idx1] - pnt0_xy))
    dist2 = float(np.linalg.norm(edge_midpoints[idx2] - pnt0_xy))
    idx = idx1 if dist1 >= dist2 else idx2
    end_anchor0 = edge_midpoints[idx]
    angle = edge_angles[idx]
    split = _pnts_group_end_split(pnts_bnd, pnt0, pd.DataFrame([end_anchor0], columns=["Cx", "Cy"]))
    if split is None:
        return None, None, angle, end_anchor0
    return split[0], split[1], angle, end_anchor0


def _pnts_bnd_end_quad(
    pnts_bnd: pd.DataFrame,
    pnts0: np.ndarray,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None, float, np.ndarray | None]:
    pnts_bnd = f_pnts_bnd_ascending(pnts_bnd)
    pnt0 = pd.DataFrame([[float((pnts0[0, 0] + pnts0[1, 0]) / 2.0), float((pnts0[0, 1] + pnts0[1, 1]) / 2.0)]], columns=["Cx", "Cy"])
    if len(pnts_bnd) == 3:
        angle = _normalize_angle(math.atan2(float(pnts0[1, 1] - pnts0[0, 1]), float(pnts0[1, 0] - pnts0[0, 0])))
        end_anchor0 = pnts_bnd.iloc[1][["Cx", "Cy"]].to_numpy(dtype=float)
        split = _pnts_group_end_split(pnts_bnd, pnt0, pd.DataFrame([end_anchor0], columns=["Cx", "Cy"]))
        if split is None:
            return None, None, angle, end_anchor0
        return split[0], split[1], angle, end_anchor0

    hull_coords = np.array(Polygon(pnts_bnd[["Cx", "Cy"]].to_numpy(dtype=float)).convex_hull.exterior.coords[:-1], dtype=float)
    if len(hull_coords) < 3:
        return _pnts_bnd_end_mbb(pnts_bnd, pnts0)
    agl0_original = math.atan2(float(pnts0[1, 1] - pnts0[0, 1]), float(pnts0[1, 0] - pnts0[0, 0]))
    agl0 = agl0_original
    perp = np.array([math.cos(agl0 + math.pi / 2.0), math.sin(agl0 + math.pi / 2.0)], dtype=float)
    pnt0_xy = pnt0.iloc[0][["Cx", "Cy"]].to_numpy(dtype=float)
    if float(((hull_coords - pnt0_xy) @ perp).sum()) < 0.0:
        agl0 += math.pi
    agl0 = _normalize_angle(agl0)
    best_value = float("inf")
    best_par: tuple[float, float, float] | None = None
    bounds = np.linspace(agl0 + math.pi / 2.0, agl0 + math.pi / 2.0 + math.pi, 11)
    for i in range(10):
        for a2 in np.linspace(bounds[i], bounds[i + 1], 9):
            for a1_ratio in np.linspace(0.01, 0.99, 9):
                for a3_ratio in np.linspace(0.01, 0.99, 9):
                    value = _pnts_bqa(hull_coords, agl0, pnt0_xy, (float(a2), float(a1_ratio), float(a3_ratio)))
                    if value < best_value:
                        best_value = value
                        best_par = (float(a2), float(a1_ratio), float(a3_ratio))
    if best_par is None or not np.isfinite(best_value):
        return _pnts_bnd_end_mbb(pnts_bnd, pnts0)
    a2, a1_ratio, a3_ratio = best_par
    a1_min = max(a2 - math.pi, agl0)
    a3_min = max(agl0 + math.pi * 2.0 - math.pi, a2)
    a1_max = min(agl0 + math.pi, a2)
    a3_max = min(a2 + math.pi, agl0 + math.pi * 2.0)
    a1 = a1_min + (a1_max - a1_min) * a1_ratio
    a3 = a3_min + (a3_max - a3_min) * a3_ratio
    ints = _pnts_bqc(hull_coords, (agl0, a1, a2, a3), pnt0_xy)
    if ints is None:
        return _pnts_bnd_end_mbb(pnts_bnd, pnts0)
    _, int12, int23, _ = ints
    angle = a2
    if _angle_included(agl0_original, angle) > math.pi / 2.0:
        angle -= math.pi
    angle = _normalize_angle(angle)
    end_anchor0 = (int12 + int23) / 2.0
    split = _pnts_group_end_split(pnts_bnd, pnt0, pd.DataFrame([end_anchor0], columns=["Cx", "Cy"]))
    if split is None:
        return None, None, angle, end_anchor0
    return split[0], split[1], angle, end_anchor0


def _pnts_bnd_end_even(
    pnts: pd.DataFrame,
    pnts0: np.ndarray,
    bool_initial: bool = True,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None, np.ndarray | None]:
    pnts_bnd = pnts[pnts["IDB"] != 0].copy()
    if len(pnts_bnd) == 3:
        ordered = f_pnts_bnd_ascending(pnts_bnd)
        if bool_initial:
            return ordered.iloc[1:3].reset_index(drop=True), ordered.iloc[0:2].reset_index(drop=True), ordered.iloc[1][["Cx", "Cy"]].to_numpy(dtype=float)
        return ordered.iloc[0:2].reset_index(drop=True), ordered.iloc[1:3].reset_index(drop=True), ordered.iloc[1][["Cx", "Cy"]].to_numpy(dtype=float)
    pnt0 = np.array([(pnts0[0, 0] + pnts0[1, 0]) / 2.0, (pnts0[0, 1] + pnts0[1, 1]) / 2.0], dtype=float)
    nonzero = pnts.loc[((pnts["Cx"] - pnt0[0]) != 0) | ((pnts["Cy"] - pnt0[1]) != 0)].copy().reset_index(drop=True)
    bnd = f_pnts_bnd_ascending(nonzero[nonzero["IDB"] != 0].copy()).reset_index(drop=True)
    if bnd.empty:
        return None, None, None
    all_angles = [math.atan2(float(row["Cy"] - pnt0[1]), float(row["Cx"] - pnt0[0])) for _, row in nonzero.iterrows()]
    bnd_angles = [math.atan2(float(row["Cy"] - pnt0[1]), float(row["Cx"] - pnt0[0])) for _, row in bnd.iterrows()]
    count_diff: list[int] = []
    for bnd_angle in bnd_angles:
        signs = []
        for angle in all_angles:
            delta = angle - bnd_angle
            if delta > math.pi:
                delta -= math.pi * 2.0
            if delta < -math.pi:
                delta += math.pi * 2.0
            if math.isclose(delta, math.pi) or math.isclose(delta, -math.pi):
                delta = 0.0
            signs.append(delta)
        count_diff.append(sum(1 for value in signs if value < 0) - sum(1 for value in signs if value > 0))
    pnt0_df = pd.DataFrame([pnt0], columns=["Cx", "Cy"])
    bnd_curve = bnd[["Cx", "Cy"]].to_numpy(dtype=float)
    unblocked = []
    for _, row in bnd.iterrows():
        ints = _pnts_bnd_end_even_ints(bnd_curve, pnt0_df, pd.DataFrame([row[["Cx", "Cy"]]], columns=["Cx", "Cy"]))
        unblocked.append(1 if ints is not None and len(ints) == 1 else 0)
    valid = [abs(diff) for diff, flag in zip(count_diff, unblocked) if flag == 1]
    if not valid:
        return None, None, None
    diff_abs_min = min(valid)
    diff_ids = [idx for idx, (diff, flag) in enumerate(zip(count_diff, unblocked)) if flag == 1 and abs(diff) == diff_abs_min]
    if count_diff[diff_ids[0]] == 0:
        if len(diff_ids) > 1:
            return None, None, None
        idx_ahead = diff_ids[0]
        idx_behind = diff_ids[0]
    else:
        idx_ahead = None
        idx_behind = None
        for i in range(len(bnd) - 1):
            if unblocked[i] == 1 and unblocked[i + 1] == 1 and count_diff[i] * count_diff[i + 1] < 0:
                idx_ahead = i + 1
                idx_behind = i
        if idx_ahead is None or idx_behind is None:
            idx_ahead = diff_ids[0]
            idx_behind = diff_ids[0]
    if bool_initial:
        right = bnd.iloc[idx_ahead:].reset_index(drop=True)
        left = bnd.iloc[: idx_behind + 1].reset_index(drop=True)
    else:
        right = bnd.iloc[: idx_behind + 1].reset_index(drop=True)
        left = bnd.iloc[idx_ahead:].reset_index(drop=True)
    return right, left, bnd.iloc[idx_ahead][["Cx", "Cy"]].to_numpy(dtype=float)


def _pnts_bnd_end_even_ints(
    pnts_curve_bnd: np.ndarray,
    pnt0: pd.DataFrame,
    pnt_bnd: pd.DataFrame,
) -> np.ndarray | None:
    segment = np.vstack([
        pnt0.iloc[0][["Cx", "Cy"]].to_numpy(dtype=float),
        pnt_bnd.iloc[0][["Cx", "Cy"]].to_numpy(dtype=float),
    ])
    return _intersection_curves_polyline(np.asarray(pnts_curve_bnd, dtype=float), segment)


def _pnts_group_end_split(
    pnts_bnd: pd.DataFrame,
    pnt0: pd.DataFrame,
    end_anchor0: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    x0 = float(pnt0.iloc[0]["Cx"])
    y0 = float(pnt0.iloc[0]["Cy"])
    lx_path = float(end_anchor0.iloc[0]["Cx"] - x0)
    ly_path = float(end_anchor0.iloc[0]["Cy"] - y0)
    ln_path = math.hypot(lx_path, ly_path)
    if ln_path == 0.0:
        return None
    ordered = f_pnts_bnd_ascending(pnts_bnd)
    reduced_rows = []
    agls = []
    for _, row in ordered.iterrows():
        lx = float(row["Cx"] - x0)
        ly = float(row["Cy"] - y0)
        ln = math.hypot(lx, ly)
        if lx != 0 or ly != 0:
            reduced_rows.append(row)
            cosine = max(-1.0, min(1.0, (lx * lx_path + ly * ly_path) / (ln * ln_path)))
            agls.append(math.acos(cosine))
    if not reduced_rows:
        return None
    reduced = pd.DataFrame(reduced_rows).reset_index(drop=True)
    idxs = [idx for idx, value in enumerate(agls) if math.isclose(value, min(agls), rel_tol=1e-9, abs_tol=1e-9)]
    if len(idxs) <= 1:
        idx = idxs[0]
    else:
        dists = [math.hypot(float(reduced.iloc[i]["Cx"] - x0), float(reduced.iloc[i]["Cy"] - y0)) for i in idxs]
        idx = idxs[int(np.argmin(dists))]
    return reduced.iloc[: idx + 1].reset_index(drop=True), reduced.iloc[idx:].reset_index(drop=True)


def _intersection_curves_polyline(curve1: np.ndarray, curve2: np.ndarray) -> np.ndarray | None:
    inter = LineString(np.asarray(curve1, dtype=float).tolist()).intersection(LineString(np.asarray(curve2, dtype=float).tolist()))
    points = _extract_intersection_points(inter)
    if not points:
        return None
    coords: list[np.ndarray] = []
    for point in points:
        xy = np.array([point.x, point.y], dtype=float)
        if not any(np.allclose(xy, existing, atol=1e-8) for existing in coords):
            coords.append(xy)
    return np.vstack(coords) if coords else None


def _intersection_curves2(pnts1: np.ndarray, pnts2: np.ndarray) -> np.ndarray | None:
    x1, y1 = float(pnts1[0, 0]), float(pnts1[0, 1])
    x2, y2 = float(pnts1[1, 0]), float(pnts1[1, 1])
    x3, y3 = float(pnts2[0, 0]), float(pnts2[0, 1])
    x4, y4 = float(pnts2[1, 0]), float(pnts2[1, 1])
    denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if math.isclose(denominator, 0.0, abs_tol=1e-12):
        return None
    intx = ((x3 - x4) * (x1 * y2 - y1 * x2) - (x1 - x2) * (x3 * y4 - y3 * x4)) / denominator
    inty = ((y3 - y4) * (x1 * y2 - y1 * x2) - (y1 - y2) * (x3 * y4 - y3 * x4)) / denominator
    return np.array([intx, inty], dtype=float)


def _pnts_bqa(pnts: np.ndarray, agl0: float, pnt0: np.ndarray, aglx: tuple[float, float, float]) -> float:
    a2, a1_ratio, a3_ratio = aglx
    a1_min = max(a2 - math.pi, agl0)
    a3_min = max(agl0 + math.pi * 2.0 - math.pi, a2)
    a1_max = min(agl0 + math.pi, a2)
    a3_max = min(a2 + math.pi, agl0 + math.pi * 2.0)
    a1 = a1_min + (a1_max - a1_min) * a1_ratio
    a3 = a3_min + (a3_max - a3_min) * a3_ratio
    ints = _pnts_bqc(pnts, (agl0, a1, a2, a3), pnt0)
    if ints is None:
        return float("inf")
    return _pnts_bqca(ints)


def _pnts_bqc(pnts: np.ndarray, agls: tuple[float, float, float, float], pnt0: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    agl0, agl1, agl2, agl3 = agls
    pnt1 = _pnts_tangency(pnts, agl0, agl1)
    pnt2 = _pnts_tangency(pnts, agl1, agl2)
    pnt3 = _pnts_tangency(pnts, agl2, agl3)
    if pnt1 is None or pnt2 is None or pnt3 is None:
        return None
    pnt0b = _pnts_pntb(pnt0, agl0, 100.0)
    pnt1b = _pnts_pntb(pnt1, agl1, 100.0)
    pnt2b = _pnts_pntb(pnt2, agl2, 100.0)
    pnt3b = _pnts_pntb(pnt3, agl3, 100.0)
    int01 = _intersection_curves2(np.vstack([pnt0, pnt0b]), np.vstack([pnt1, pnt1b]))
    int12 = _intersection_curves2(np.vstack([pnt1, pnt1b]), np.vstack([pnt2, pnt2b]))
    int23 = _intersection_curves2(np.vstack([pnt2, pnt2b]), np.vstack([pnt3, pnt3b]))
    int30 = _intersection_curves2(np.vstack([pnt3, pnt3b]), np.vstack([pnt0, pnt0b]))
    if int01 is None or int12 is None or int23 is None or int30 is None:
        return None
    return int01, int12, int23, int30


def _pnts_bqca(ints: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]) -> float:
    int01, int12, int23, int30 = ints
    a = float(np.linalg.norm(int12 - int01))
    b = float(np.linalg.norm(int23 - int12))
    c = float(np.linalg.norm(int30 - int23))
    d = float(np.linalg.norm(int01 - int30))
    p = float(np.linalg.norm(int23 - int01))
    s1 = (a + b + p) / 2.0
    s2 = (c + d + p) / 2.0
    return math.sqrt(max(s1 * (s1 - a) * (s1 - b) * (s1 - p), 0.0)) + math.sqrt(max(s2 * (s2 - c) * (s2 - d) * (s2 - p), 0.0))


def _pnts_tangency(pnts: np.ndarray, agl0: float, agl: float) -> np.ndarray | None:
    perp = np.array([math.cos(agl - math.pi / 2.0), math.sin(agl - math.pi / 2.0)], dtype=float)
    direction = np.array([math.cos(agl0), math.sin(agl0)], dtype=float)
    dot = float(perp @ direction)
    if math.isclose(dot, 0.0, abs_tol=1e-12):
        return None
    if dot < 0.0:
        perp = -perp
    dots = np.asarray(pnts, dtype=float) @ perp
    return np.asarray(pnts, dtype=float)[int(np.argmax(dots))]


def _pnts_pntb(pnt: np.ndarray, agl: float, step: float) -> np.ndarray:
    pnt = np.asarray(pnt, dtype=float)
    dir_x = math.cos(agl)
    dir_y = math.sin(agl)
    if math.isclose(dir_x, 0.0, abs_tol=1e-12):
        return np.array([pnt[0], pnt[1] + step], dtype=float)
    pnt2_x = pnt[0] + step
    pnt2_y = dir_y / dir_x * (pnt2_x - pnt[0]) + pnt[1]
    return np.array([pnt2_x, pnt2_y], dtype=float)


def _mark_stopped_group(group: pd.DataFrame, grp_id: int, marker: str) -> pd.DataFrame:
    stopped = group.copy()
    stopped.loc[:, "SMG"] = marker
    stopped.loc[:, "grp"] = grp_id
    return stopped


def _reassign_clusters(
    original: pd.DataFrame,
    upper: pd.DataFrame,
    lower: pd.DataFrame,
    upstream_point: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    upper_clusters = f_pnts_clustering(upper)
    lower_clusters = f_pnts_clustering(lower)
    if not upper_clusters or not lower_clusters:
        return upper, lower, False
    if len(upper_clusters) >= 2 or len(lower_clusters) >= 2:
        upstream_id = int(upstream_point["ID"].iloc[0])
        chosen_upper = None
        overflow_to_lower: list[pd.DataFrame] = []
        for cluster in upper_clusters:
            if upstream_id in set(cluster["ID"]):
                chosen_upper = cluster
            else:
                overflow_to_lower.append(cluster)
        if chosen_upper is None:
            return upper, lower, False
        lower_clusters = sorted(lower_clusters, key=len, reverse=True)
        chosen_lower = lower_clusters[0]
        overflow_to_upper = lower_clusters[1:]
        upper = pd.concat([chosen_upper] + overflow_to_upper, ignore_index=True)
        lower = pd.concat([chosen_lower] + overflow_to_lower, ignore_index=True)

    if len(upper) + len(lower) != len(original):
        return upper, lower, False
    if len(f_pnts_clustering(upper)) != 1 or len(f_pnts_clustering(lower)) != 1:
        return upper, lower, False
    return upper, lower, True


def _validate_boundary_side_count(group: pd.DataFrame, idb_max: int, initial: bool, distal: bool) -> bool:
    bnd = group[group["IDB"] != 0].copy()
    sides = f_pnts_bnd_sides(bnd, idb_max)
    if initial or distal:
        return len(sides) == 1
    return len(sides) == 2


def _merge_small_end_groups(
    groups: list[pd.DataFrame],
    total_count: int,
    min_end_ratio_initial: float,
    min_end_ratio_distal: float,
) -> list[pd.DataFrame]:
    groups = [group.copy() for group in groups]
    if len(groups) >= 3 and min_end_ratio_initial > 0:
        end_count = len(groups[0])
        merge_until = 0
        for idx in range(1, len(groups)):
            end_count += len(groups[idx])
            if end_count / total_count > min_end_ratio_initial:
                merge_until = idx
                break
        if merge_until > 1:
            merged = pd.concat(groups[:merge_until], ignore_index=True)
            merged.loc[:, "SMG"] = "SGR"
            groups = [merged] + groups[merge_until:]
    if len(groups) >= 3 and min_end_ratio_distal > 0:
        end_count = len(groups[-1])
        merge_from = len(groups) - 1
        for idx in range(len(groups) - 2, -1, -1):
            end_count += len(groups[idx])
            if end_count / total_count > min_end_ratio_distal:
                merge_from = idx + 1
                break
        if merge_from < len(groups) - 1:
            merged = pd.concat(groups[merge_from:], ignore_index=True)
            merged.loc[:, "SMG"] = "SGR"
            groups = groups[:merge_from] + [merged]
    return groups


def _renumber_groups(groups: list[pd.DataFrame]) -> list[pd.DataFrame]:
    renumbered = []
    for idx, group in enumerate(groups, start=1):
        updated = group.copy()
        updated.loc[:, "grp"] = idx
        renumbered.append(updated)
    return renumbered


def _resolve_end_anchor(group: pd.DataFrame, polygon: Polygon, config: Any, initial: bool) -> np.ndarray | None:
    if isinstance(config, np.ndarray) and config.size >= 2:
        return config[:2].astype(float)
    if isinstance(config, (list, tuple, np.ndarray)) and len(config) == 2:
        return np.array(config, dtype=float)
    return None


def _resolve_end_anchor_seed(
    group: pd.DataFrame,
    polygon: Polygon,
    config: Any,
    initial: bool,
    idb_max: int,
) -> np.ndarray | None:
    explicit_anchor = _resolve_end_anchor(group, polygon, config, initial)
    if explicit_anchor is not None:
        return explicit_anchor
    return _group_anchor_center(group, idb_max)


def _resolve_final_end_anchor(
    group: pd.DataFrame,
    polygon: Polygon,
    anchors: pd.DataFrame,
    config: Any,
    initial: bool,
) -> np.ndarray | None:
    explicit_anchor = _resolve_end_anchor(group, polygon, config, initial)
    if explicit_anchor is not None:
        return explicit_anchor

    anchor_coords = anchors[["Cx", "Cy"]].to_numpy(dtype=float)
    if len(anchor_coords) == 0:
        return None
    if initial:
        anchor_xy = anchor_coords[0]
        inward_xy = anchor_coords[1] if len(anchor_coords) >= 2 else group[["Cx", "Cy"]].mean().to_numpy(dtype=float)
    else:
        anchor_xy = anchor_coords[-1]
        inward_xy = anchor_coords[-2] if len(anchor_coords) >= 2 else group[["Cx", "Cy"]].mean().to_numpy(dtype=float)

    extended_anchor = _extend_anchor_to_polygon(anchor_xy, inward_xy, polygon)
    if extended_anchor is not None:
        return extended_anchor
    return np.array(anchor_xy, dtype=float)


def _extend_anchor_to_polygon(anchor_xy: np.ndarray, inward_xy: np.ndarray, polygon: Polygon) -> np.ndarray | None:
    anchor_xy = np.array(anchor_xy, dtype=float)
    inward_xy = np.array(inward_xy, dtype=float)
    direction = anchor_xy - inward_xy
    norm = float(np.linalg.norm(direction))
    if norm == 0.0:
        center = np.array([polygon.centroid.x, polygon.centroid.y], dtype=float)
        direction = anchor_xy - center
        norm = float(np.linalg.norm(direction))
    if norm == 0.0:
        return None
    direction = direction / norm
    segment = _line_segment_through_point(polygon, anchor_xy, direction)
    if segment is None:
        return None
    projections = (segment - anchor_xy) @ direction
    positive = np.where(projections > 1e-8)[0]
    if len(positive) == 0:
        return None
    nearest_idx = positive[int(np.argmin(projections[positive]))]
    return np.array(segment[nearest_idx], dtype=float)


def _legacy_end_anchor_config(group: pd.DataFrame, polygon: Polygon, config: Any, initial: bool) -> np.ndarray | None:
    code = int(config)
    bnd = group[group["IDB"] != 0].copy()
    if bnd.empty:
        centroid = group[["Cx", "Cy"]].mean().to_numpy(dtype=float)
        return centroid
    if code == 1:
        return group.sort_values("zip", ascending=not initial).iloc[0][["Cx", "Cy"]].to_numpy(dtype=float)
    if code == 2:
        ordered = f_pnts_bnd_ascending(bnd)
        edge = ordered.iloc[[0, -1]][["Cx", "Cy"]].to_numpy(dtype=float)
        return edge.mean(axis=0)
    if code == 3:
        center = np.array([polygon.centroid.x, polygon.centroid.y], dtype=float)
        coords = bnd[["Cx", "Cy"]].to_numpy(dtype=float)
        dist = np.linalg.norm(coords - center, axis=1)
        return coords[int(np.argmax(dist))]
    raise ValueError("Unsupported end-anchor configuration. Use 1, 2, 3, or explicit coordinates.")


def _resolve_end_strip_constraint(
    group: pd.DataFrame,
    polygon: Polygon,
    anchor_xy: np.ndarray,
    reference_xy: np.ndarray,
    config: Any,
    initial: bool,
) -> dict[str, Any] | None:
    if isinstance(config, np.ndarray) and config.size >= 5:
        return {
            "mode": "segment",
            "segment": config[:4].astype(float).reshape(2, 2),
            "angle": float(config[4]),
        }
    if isinstance(config, (list, tuple, np.ndarray)) and len(config) == 5:
        array = np.array(config, dtype=float)
        return {
            "mode": "segment",
            "segment": array[:4].reshape(2, 2),
            "angle": float(array[4]),
        }

    code = int(config)
    bnd = group[group["IDB"] != 0].copy()
    if bnd.empty:
        return None
    if code == 1:
        angle, anchor = _boundary_mrr_direction(bnd, reference_xy)
        return {"mode": "angle", "angle": angle, "anchor": anchor}
    if code == 2:
        angle, anchor = _boundary_even_direction(group, anchor_xy, reference_xy, initial)
        return {"mode": "angle", "angle": angle, "anchor": anchor}
    if code == 3:
        vector = np.array(reference_xy, dtype=float) - np.array(anchor_xy, dtype=float)
        angle = _normalize_angle(math.atan2(vector[1], vector[0]) + math.pi / 2.0)
        return {"mode": "angle", "angle": angle, "anchor": np.array(anchor_xy, dtype=float)}
    raise ValueError("Unsupported end-strip configuration. Use 1, 2, 3, or explicit line geometry.")


def _group_anchor_center(group: pd.DataFrame, idb_max: int) -> np.ndarray:
    bnd = group[group["IDB"] != 0].copy()
    sides = f_pnts_bnd_sides(bnd, idb_max)
    centers = f_pnts_bnd_sides_centerLs(sides)
    if centers is not None:
        return centers[2]
    return group[["Cx", "Cy"]].mean().to_numpy(dtype=float)


def _order_groups_for_path(
    groups: list[pd.DataFrame],
    idb_max: int,
    start_anchor: np.ndarray,
    end_anchor: np.ndarray,
) -> list[pd.DataFrame]:
    groups = [group.copy() for group in groups]
    if len(groups) <= 2:
        return _renumber_groups(groups)

    direction = _path_direction_unit(start_anchor, end_anchor, groups)
    interior: list[tuple[float, float, pd.DataFrame]] = []
    for group in groups[1:-1]:
        center = _group_anchor_center(group, idb_max)
        projection = _path_projection(center, start_anchor, direction)
        elevation = float(f_pnts_zip(group, center[0], center[1]))
        interior.append((-elevation, projection, group.copy()))

    interior.sort(key=lambda item: (item[0], item[1]))
    ordered = [groups[0]] + [item[2] for item in interior] + [groups[-1]]
    return _renumber_groups(ordered)


def _path_direction_unit(start_xy: np.ndarray, end_xy: np.ndarray, groups: list[pd.DataFrame]) -> np.ndarray:
    direction = np.array(end_xy, dtype=float) - np.array(start_xy, dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm > 0:
        return direction / norm
    if groups:
        first_center = groups[0][["Cx", "Cy"]].mean().to_numpy(dtype=float)
        last_center = groups[-1][["Cx", "Cy"]].mean().to_numpy(dtype=float)
        direction = last_center - first_center
        norm = float(np.linalg.norm(direction))
        if norm > 0:
            return direction / norm
    return np.array([1.0, 0.0], dtype=float)


def _path_projection(point_xy: np.ndarray, origin_xy: np.ndarray, direction_xy: np.ndarray) -> float:
    point_xy = np.array(point_xy, dtype=float)
    origin_xy = np.array(origin_xy, dtype=float)
    direction_xy = np.array(direction_xy, dtype=float)
    return float((point_xy - origin_xy) @ direction_xy)


def _select_progressive_candidate(
    candidate_xy: np.ndarray,
    candidate_z: np.ndarray,
    target_xy: np.ndarray,
    prev_xy: np.ndarray,
    prev_z: float,
    origin_xy: np.ndarray,
    direction_xy: np.ndarray,
) -> np.ndarray:
    if len(candidate_xy) == 0:
        return np.array(target_xy, dtype=float)

    target_xy = np.array(target_xy, dtype=float)
    prev_xy = np.array(prev_xy, dtype=float)
    origin_xy = np.array(origin_xy, dtype=float)
    direction_xy = np.array(direction_xy, dtype=float)
    prev_proj = _path_projection(prev_xy, origin_xy, direction_xy)
    target_proj = _path_projection(target_xy, origin_xy, direction_xy)

    scores: list[float] = []
    for coord, z_value in zip(candidate_xy, candidate_z):
        proj = _path_projection(coord, origin_xy, direction_xy)
        backtrack = max(prev_proj - proj, 0.0)
        uphill = max(float(z_value) - float(prev_z), 0.0)
        target_gap = abs(proj - max(target_proj, prev_proj))
        offset = float(np.linalg.norm(coord - target_xy))
        scores.append(uphill * 1e6 + backtrack * 1e6 + target_gap * 10.0 + offset)

    return np.array(candidate_xy[int(np.argmin(scores))], dtype=float)


def _select_progressive_group_anchor(
    group: pd.DataFrame,
    default_xy: np.ndarray,
    prev_xy: np.ndarray,
    prev_z: float,
    origin_xy: np.ndarray,
    direction_xy: np.ndarray,
) -> np.ndarray:
    point_xy = group[["Cx", "Cy"]].to_numpy(dtype=float)
    point_z = group["zip"].to_numpy(dtype=float)
    default_xy = np.array(default_xy, dtype=float).reshape(1, 2)
    default_z = np.array([f_pnts_zip(group, float(default_xy[0, 0]), float(default_xy[0, 1]))], dtype=float)
    candidate_xy = np.vstack([default_xy, point_xy])
    candidate_z = np.concatenate([default_z, point_z])
    return _select_progressive_candidate(candidate_xy, candidate_z, default_xy[0], prev_xy, prev_z, origin_xy, direction_xy)


def _select_progressive_terminal_anchor(
    group: pd.DataFrame,
    default_xy: np.ndarray,
    prev_xy: np.ndarray,
    prev_z: float,
    origin_xy: np.ndarray,
    direction_xy: np.ndarray,
) -> np.ndarray:
    bnd = group[group["IDB"] != 0].copy()
    source = bnd if not bnd.empty else group
    point_xy = source[["Cx", "Cy"]].to_numpy(dtype=float)
    point_z = source["zip"].to_numpy(dtype=float)
    default_xy = np.array(default_xy, dtype=float).reshape(1, 2)
    default_z = np.array([f_pnts_zip(group, float(default_xy[0, 0]), float(default_xy[0, 1]))], dtype=float)
    candidate_xy = np.vstack([default_xy, point_xy])
    candidate_z = np.concatenate([default_z, point_z])
    return _select_progressive_candidate(candidate_xy, candidate_z, default_xy[0], prev_xy, prev_z, origin_xy, direction_xy)


def _remove_duplicate_path_vertices(anchors: pd.DataFrame) -> pd.DataFrame:
    keep_rows = [0]
    coords = anchors[["Cx", "Cy"]].to_numpy(dtype=float)
    for idx in range(1, len(coords)):
        if not np.allclose(coords[idx], coords[keep_rows[-1]]):
            keep_rows.append(idx)
    return anchors.iloc[keep_rows].reset_index(drop=True)


def _is_horizontal_plane(pnts: pd.DataFrame) -> bool:
    row = pnts.iloc[0]
    return bool(np.isclose(row["Dx"], 0.0) and np.isclose(row["Dy"], 0.0) and np.isclose(row["Dz"], 0.0))


def _ensure_polygon(geometry: Any) -> Polygon | None:
    if geometry is None or geometry.is_empty:
        return None
    if geometry.geom_type == "Polygon":
        return _repair_polygon(geometry)
    if geometry.geom_type == "MultiPolygon":
        polygons = list(geometry.geoms)
        polygon = max(polygons, key=lambda item: item.area)
        return _repair_polygon(polygon)
    return None


def _repair_polygon(geometry: Polygon) -> Polygon | None:
    if geometry is None or geometry.is_empty:
        return None
    repaired = geometry
    if not repaired.is_valid:
        repaired = repaired.buffer(0)
    if repaired.is_empty:
        return None
    if repaired.geom_type == "Polygon":
        return repaired
    if repaired.geom_type == "MultiPolygon":
        polygons = [geom for geom in repaired.geoms if isinstance(geom, Polygon) and geom.area > 0]
        if not polygons:
            return None
        return max(polygons, key=lambda item: item.area)
    return None


def _path_stations(path_line: LineString, min_length: float) -> list[float]:
    length = float(path_line.length)
    if length == 0:
        return [0.0]
    if min_length <= 0:
        return [0.0, length]
    count = max(1, int(math.ceil(length / min_length)))
    stations = np.linspace(0.0, length, count + 1)
    return [float(item) for item in stations]


def _line_tangent(path_line: LineString, chainage: float) -> np.ndarray:
    delta = min(max(path_line.length / 1000.0, 1e-6), 1.0)
    s0 = max(0.0, chainage - delta)
    s1 = min(path_line.length, chainage + delta)
    p0 = path_line.interpolate(s0)
    p1 = path_line.interpolate(s1)
    tangent = np.array([p1.x - p0.x, p1.y - p0.y], dtype=float)
    norm = np.linalg.norm(tangent)
    if norm == 0.0:
        return np.array([1.0, 0.0], dtype=float)
    return tangent / norm


def _perpendicular_segment(polygon: Polygon, point: Point, tangent: np.ndarray) -> np.ndarray | None:
    normal = np.array([-tangent[1], tangent[0]], dtype=float)
    span = max(polygon.bounds[2] - polygon.bounds[0], polygon.bounds[3] - polygon.bounds[1]) * 2.0 + 1.0
    line = LineString([
        (point.x - normal[0] * span, point.y - normal[1] * span),
        (point.x + normal[0] * span, point.y + normal[1] * span),
    ])
    clipped = polygon.intersection(line)
    segment = _longest_linestring(clipped)
    if segment is None or segment.length == 0:
        return None
    coords = np.array(segment.coords, dtype=float)
    if len(coords) < 2:
        return None
    return np.vstack([coords[0], coords[-1]])


def _segment_from_end_constraint(polygon: Polygon, center: np.ndarray, constraint: dict[str, Any]) -> np.ndarray | None:
    if constraint is None:
        return None
    if constraint.get("mode") == "segment":
        segment = np.array(constraint["segment"], dtype=float)
        clipped = polygon.intersection(LineString(segment.tolist()))
        line = _longest_linestring(clipped)
        if line is not None and line.length > 0:
            coords = np.array(line.coords, dtype=float)
            if len(coords) >= 2:
                return _extend_cut_segment(np.vstack([coords[0], coords[-1]]), polygon)
        if _segment_hits_polygon_twice(polygon, segment):
            return _extend_cut_segment(segment, polygon)
        return _extend_segment_to_polygon(polygon, segment)

    angle = float(constraint.get("angle", 0.0))
    anchor = np.array(constraint.get("anchor", center), dtype=float)
    direction = np.array([math.cos(angle), math.sin(angle)], dtype=float)
    return _line_segment_through_point_onr(polygon, anchor, direction)


def _longest_linestring(geometry: Any) -> LineString | None:
    if geometry is None or geometry.is_empty:
        return None
    if isinstance(geometry, LineString):
        return geometry
    if isinstance(geometry, MultiLineString):
        return max(geometry.geoms, key=lambda item: item.length, default=None)
    if isinstance(geometry, GeometryCollection):
        lines = [geom for geom in geometry.geoms if isinstance(geom, LineString)]
        return max(lines, key=lambda item: item.length, default=None)
    return None


def _segment_is_valid_cut(polygon: Polygon, path_line: LineString, segment: np.ndarray, is_end: bool = False) -> bool:
    if segment is None:
        return False
    if not _segment_hits_polygon_twice(polygon, segment):
        return False
    if not _segment_cuts_path(path_line, segment, is_end=is_end):
        return False
    split_polygons = _split_polygon_segment(polygon, segment, extend_ratio=1e-3)
    if len(split_polygons) != 2:
        return False
    polygon_area = float(polygon.area)
    if polygon_area <= 0:
        return False
    return min(part.area for part in split_polygons) > polygon_area / 1e6


def _extend_cut_segment(segment: np.ndarray, polygon: Polygon, factor: float = 1.0) -> np.ndarray:
    segment = np.array(segment, dtype=float)
    if segment.shape != (2, 2):
        return segment
    vector = segment[1] - segment[0]
    norm = float(np.linalg.norm(vector))
    if norm == 0.0:
        return segment
    direction = vector / norm
    extension = max(1000.0 * factor, 1.0)
    return np.vstack([segment[0] - direction * extension, segment[1] + direction * extension])


def _segment_hits_polygon_twice(polygon: Polygon, segment: np.ndarray) -> bool:
    line = LineString(np.array(segment, dtype=float).tolist())
    inter = polygon.boundary.intersection(line)
    points = _extract_intersection_points(inter)
    return len(points) >= 2


def _segment_cuts_path(path_line: LineString, segment: np.ndarray, is_end: bool = False) -> bool:
    line = LineString(np.array(segment, dtype=float).tolist())
    inter = path_line.intersection(line)
    points = _extract_intersection_points(inter)
    if is_end:
        return len(points) >= 1
    return len(points) >= 1


def _extract_intersection_points(geometry: Any) -> list[Point]:
    if geometry is None or geometry.is_empty:
        return []
    if isinstance(geometry, Point):
        return [geometry]
    if geometry.geom_type == "MultiPoint":
        return list(geometry.geoms)
    if isinstance(geometry, LineString):
        coords = list(geometry.coords)
        return [Point(coords[0]), Point(coords[-1])] if coords else []
    if isinstance(geometry, MultiLineString):
        points: list[Point] = []
        for geom in geometry.geoms:
            coords = list(geom.coords)
            if coords:
                points.extend([Point(coords[0]), Point(coords[-1])])
        return points
    if isinstance(geometry, GeometryCollection):
        points: list[Point] = []
        for geom in geometry.geoms:
            points.extend(_extract_intersection_points(geom))
        return points
    return []


def _nearest_boundary_intersection_along_ray(polygon: Polygon, origin_xy: np.ndarray, direction: np.ndarray) -> np.ndarray | None:
    direction = np.array(direction, dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm == 0.0:
        return None
    direction = direction / norm
    span = max(polygon.bounds[2] - polygon.bounds[0], polygon.bounds[3] - polygon.bounds[1]) * 4.0 + 1.0
    ray = LineString([
        (float(origin_xy[0]), float(origin_xy[1])),
        (float(origin_xy[0] + direction[0] * span), float(origin_xy[1] + direction[1] * span)),
    ])
    inter = polygon.boundary.intersection(ray)
    points = _extract_intersection_points(inter)
    candidates: list[tuple[float, np.ndarray]] = []
    for point in points:
        coords = np.array([point.x, point.y], dtype=float)
        distance = float(np.dot(coords - origin_xy, direction))
        if distance > 1e-8:
            candidates.append((distance, coords))
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0])
    return candidates[0][1]


def _extended_boundary_endpoint_along_ray_onr(
    polygon: Polygon,
    origin_xy: np.ndarray,
    direction: np.ndarray,
    extnd: float = 1000.0,
) -> np.ndarray | None:
    direction = np.array(direction, dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm == 0.0:
        return None
    direction = direction / norm
    span = max(polygon.bounds[2] - polygon.bounds[0], polygon.bounds[3] - polygon.bounds[1]) * 4.0 + 1.0
    ray = LineString([
        (float(origin_xy[0]), float(origin_xy[1])),
        (float(origin_xy[0] + direction[0] * span), float(origin_xy[1] + direction[1] * span)),
    ])
    inter = polygon.boundary.intersection(ray)
    points = _extract_intersection_points(inter)
    candidates: list[tuple[float, np.ndarray]] = []
    for point in points:
        coords = np.array([point.x, point.y], dtype=float)
        distance = float(np.dot(coords - origin_xy, direction))
        if distance > 1e-8:
            candidates.append((distance, coords))
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0])
    nearest = candidates[0][1]
    if len(candidates) == 1:
        return nearest + (nearest - origin_xy) / float(extnd)
    second = candidates[1][1]
    return nearest + (second - nearest) / float(extnd)


def _line_segment_through_point_onr(polygon: Polygon, point_xy: np.ndarray, direction: np.ndarray) -> np.ndarray | None:
    direction = np.array(direction, dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm == 0.0:
        return None
    direction = direction / norm
    point_a = _extended_boundary_endpoint_along_ray_onr(polygon, point_xy, -direction)
    point_b = _extended_boundary_endpoint_along_ray_onr(polygon, point_xy, direction)
    if point_a is None or point_b is None:
        return None
    return np.vstack([point_a, point_b])


def _split_polygon_segment(polygon: Polygon, segment: np.ndarray, extend_ratio: float = 0.0) -> list[Polygon]:
    segment = np.array(segment, dtype=float)
    if extend_ratio > 0:
        vector = segment[1] - segment[0]
        segment = np.vstack([segment[0] - vector * extend_ratio, segment[1] + vector * extend_ratio])
    line = LineString(segment.tolist())
    result = shapely_split(polygon, line)
    polygons: list[Polygon] = []
    for geom in result.geoms:
        if not isinstance(geom, Polygon) or geom.area <= 0:
            continue
        repaired = _repair_polygon(geom)
        if repaired is not None and repaired.area > 0:
            polygons.append(repaired)
    return polygons


def _extend_segment_to_polygon(polygon: Polygon, segment: np.ndarray) -> np.ndarray | None:
    segment = np.array(segment, dtype=float)
    center = segment.mean(axis=0)
    vector = segment[1] - segment[0]
    norm = np.linalg.norm(vector)
    if norm == 0:
        return None
    direction = vector / norm
    return _line_segment_through_point(polygon, center, direction)


def _line_segment_through_point(polygon: Polygon, point_xy: np.ndarray, direction: np.ndarray) -> np.ndarray | None:
    direction = np.array(direction, dtype=float)
    norm = np.linalg.norm(direction)
    if norm == 0:
        return None
    direction = direction / norm
    point_a = _nearest_boundary_intersection_along_ray(polygon, point_xy, -direction)
    point_b = _nearest_boundary_intersection_along_ray(polygon, point_xy, direction)
    if point_a is None or point_b is None:
        return None
    return _extend_cut_segment(np.vstack([point_a, point_b]), polygon)


def _boundary_mrr_direction(bnd: pd.DataFrame, reference_xy: np.ndarray) -> tuple[float, np.ndarray]:
    coords = bnd[["Cx", "Cy"]].to_numpy(dtype=float)
    if len(coords) < 3:
        return _boundary_even_direction(bnd, coords.mean(axis=0), reference_xy, True)
    hull = Polygon(coords).convex_hull
    rect = hull.minimum_rotated_rectangle
    rect_coords = np.array(rect.exterior.coords[:-1], dtype=float)
    ref_vector = coords.mean(axis=0) - np.array(reference_xy, dtype=float)
    ref_angle = math.atan2(ref_vector[1], ref_vector[0])
    candidates: list[tuple[float, np.ndarray, float]] = []
    for idx in range(len(rect_coords)):
        p0 = rect_coords[idx]
        p1 = rect_coords[(idx + 1) % len(rect_coords)]
        edge = p1 - p0
        if np.linalg.norm(edge) == 0:
            continue
        angle = _normalize_angle(math.atan2(edge[1], edge[0]))
        midpoint = (p0 + p1) / 2.0
        candidates.append((_angle_difference(angle, ref_angle), midpoint, angle))
    if not candidates:
        return _boundary_even_direction(bnd, coords.mean(axis=0), reference_xy, True)
    _, midpoint, angle = min(candidates, key=lambda item: item[0])
    return angle, midpoint


def _boundary_even_direction(group: pd.DataFrame, anchor_xy: np.ndarray, reference_xy: np.ndarray, initial: bool) -> tuple[float, np.ndarray]:
    if isinstance(group, pd.DataFrame) and "IDB" in group.columns:
        bnd = group[group["IDB"] != 0].copy()
    else:
        bnd = group.copy()
    bnd = f_pnts_bnd_ascending(bnd)
    coords = bnd[["Cx", "Cy"]].to_numpy(dtype=float)
    if len(coords) == 0:
        vec = np.array(reference_xy, dtype=float) - np.array(anchor_xy, dtype=float)
        angle = _normalize_angle(math.atan2(vec[1], vec[0]) + math.pi / 2.0)
        return angle, np.array(anchor_xy, dtype=float)
    center = np.array(anchor_xy, dtype=float)
    best_idx = 0
    best_score = float("inf")
    for idx, candidate in enumerate(coords):
        vec = candidate - center
        if np.linalg.norm(vec) == 0:
            continue
        side = np.cross(np.column_stack([coords - center, np.zeros(len(coords))]), np.array([vec[0], vec[1], 0.0]))[:, 2]
        left = int(np.sum(side > 0))
        right = int(np.sum(side < 0))
        score = abs(left - right)
        if score < best_score:
            best_score = score
            best_idx = idx
    end_anchor = coords[best_idx]
    vec = end_anchor - center
    angle = _normalize_angle(math.atan2(vec[1], vec[0]))
    return angle, end_anchor


def _ordered_strip_endpoints(segment: np.ndarray, center: np.ndarray, tangent: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    p0 = np.array(segment[0], dtype=float)
    p1 = np.array(segment[1], dtype=float)
    cross0 = np.cross(np.array([tangent[0], tangent[1], 0.0]), np.array([p0[0] - center[0], p0[1] - center[1], 0.0]))[2]
    if cross0 >= 0:
        return p0, p1
    return p1, p0


def _sample_dem(src: rasterio.io.DatasetReader, x: float, y: float) -> float:
    value = list(src.sample([(x, y)]))[0][0]
    return float(value) if np.isfinite(value) else float("nan")


def _write_points_products(pnts: pd.DataFrame, crs: Any, prefix: Path) -> None:
    _write_excel(pnts.drop(columns=[col for col in ["_row", "_col", "_a", "_b", "_c"] if col in pnts.columns]), prefix.with_suffix(".xlsx"))
    gdf = gpd.GeoDataFrame(
        pnts.copy(),
        geometry=gpd.points_from_xy(pnts["Cx"], pnts["Cy"]),
        crs=crs,
    )
    _safe_write_gdf(gdf, prefix)


def _safe_write_gdf(gdf: gpd.GeoDataFrame, prefix: Path) -> None:
    if gdf is None or gdf.empty:
        return
    prefix.parent.mkdir(parents=True, exist_ok=True)
    gdf = gdf.copy()
    for column in list(gdf.columns):
        if column.startswith("_"):
            gdf = gdf.drop(columns=column)
    gdf.to_file(prefix.with_suffix(".shp"), driver="ESRI Shapefile")


def _write_excel(df: pd.DataFrame | None, file_path: Path) -> None:
    if df is None:
        return
    file_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(file_path, index=False)


def _config_to_string(value: Any) -> str:
    if isinstance(value, np.ndarray):
        return ",".join(f"{item:.6f}" for item in value.tolist())
    if isinstance(value, (list, tuple)):
        return ",".join(str(item) for item in value)
    return str(value)


def _normalize_angle(angle: float) -> float:
    while angle > math.pi:
        angle -= math.pi * 2.0
    while angle < -math.pi:
        angle += math.pi * 2.0
    return angle


def _angle_difference(angle1: float, angle2: float) -> float:
    diff = abs(_normalize_angle(angle1 - angle2))
    return min(diff, abs(math.pi - diff))


def _build_group_boundary_events(
    list_pnts: list[pd.DataFrame],
    anchors: pd.DataFrame,
    path_line: LineString,
    idb_max: int,
) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    for idx, group in enumerate(list_pnts, start=1):
        segments = _group_boundary_segments(group, anchors, idx, idb_max)
        if not segments:
            continue
        include_upper = idx != 1
        include_lower = idx != len(list_pnts)
        for key, label in (("upper", "Intersection: Group boundary upper"), ("lower", "Intersection: Group boundary lower")):
            if key == "upper" and not include_upper:
                continue
            if key == "lower" and not include_lower:
                continue
            segment = segments.get(key)
            if segment is None:
                continue
            intersection = _polyline_segment_intersection(path_line, segment)
            if intersection is None:
                continue
            chainage = float(path_line.project(Point(intersection[0], intersection[1])))
            events.append(
                {
                    "Cx": float(intersection[0]),
                    "Cy": float(intersection[1]),
                    "Czp": 0.0,
                    "type": label,
                    "grp": idx,
                    "station": 0,
                    "chainage": chainage,
                    "angle": _segment_angle(segment),
                }
            )
    events.sort(key=lambda item: item["chainage"])
    return events


def _group_boundary_segments(group: pd.DataFrame, anchors: pd.DataFrame, group_idx: int, idb_max: int) -> dict[str, np.ndarray] | None:
    bnd = group[group["IDB"] != 0].copy()
    sides = f_pnts_bnd_sides(bnd, idb_max)
    if len(sides) < 2:
        return None

    anchor_coords = anchors[["Cx", "Cy"]].to_numpy(dtype=float)
    if len(anchor_coords) < 2:
        return None

    grp_values = anchors["grp"].to_numpy(dtype=int)
    matched = np.where(grp_values == int(group_idx))[0]
    if matched.size:
        center_idx = int(matched[0])
    else:
        center_idx = min(max(group_idx, 1), len(anchor_coords) - 2)

    prev_idx = max(center_idx - 1, 0)
    next_idx = min(center_idx + 1, len(anchor_coords) - 1)
    center = anchor_coords[center_idx]
    prev_anchor = anchor_coords[prev_idx]
    next_anchor = anchor_coords[next_idx]

    if np.allclose(next_anchor, prev_anchor):
        group_center = group[["Cx", "Cy"]].mean().to_numpy(dtype=float)
        if not np.allclose(prev_anchor, group_center):
            next_anchor = group_center
        elif len(anchor_coords) >= 2:
            prev_anchor = anchor_coords[0]
            next_anchor = anchor_coords[-1]

    direction = next_anchor - prev_anchor
    norm = np.linalg.norm(direction)
    if norm == 0:
        direction = next_anchor - center
        norm = np.linalg.norm(direction)
    if norm == 0:
        direction = np.array([1.0, 0.0], dtype=float)
    else:
        direction = direction / norm
    lateral_sides = []
    for side in sides[:2]:
        coords = side[["Cx", "Cy"]].to_numpy(dtype=float)
        centroid = coords.mean(axis=0)
        side_sign = np.cross(np.array([direction[0], direction[1], 0.0]), np.array([centroid[0] - center[0], centroid[1] - center[1], 0.0]))[2]
        projections = (coords - center) @ direction
        lateral_sides.append(
            {
                "coords": coords,
                "sign": side_sign,
                "upstream": coords[int(np.argmin(projections))],
                "downstream": coords[int(np.argmax(projections))],
            }
        )
    lateral_sides.sort(key=lambda item: item["sign"])
    if len(lateral_sides) != 2:
        return None
    right_side, left_side = lateral_sides[0], lateral_sides[1]
    return {
        "upper": np.vstack([right_side["upstream"], left_side["upstream"]]),
        "lower": np.vstack([left_side["downstream"], right_side["downstream"]]),
    }


def _segment_angle(segment: np.ndarray) -> float:
    segment = np.array(segment, dtype=float)
    vector = segment[1] - segment[0]
    return _normalize_angle(math.atan2(vector[1], vector[0]))


def _polyline_segment_intersection(path_line: LineString, segment: np.ndarray) -> np.ndarray | None:
    line = LineString(np.array(segment, dtype=float).tolist())
    inter = path_line.intersection(line)
    points = _extract_intersection_points(inter)
    if not points:
        return None
    path_points = np.array([[point.x, point.y] for point in points], dtype=float)
    chainages = [path_line.project(Point(row[0], row[1])) for row in path_points]
    return path_points[int(np.argmin(chainages))]


def _path_cut_from_constraint(path_line: LineString, polygon: Polygon, constraint: dict[str, Any] | None) -> dict[str, Any] | None:
    if constraint is None or constraint.get("mode") != "segment":
        return None
    segment = _segment_from_end_constraint(polygon, np.zeros(2, dtype=float), constraint)
    if segment is None:
        return None
    point = _polyline_segment_intersection(path_line, segment)
    if point is None:
        return None
    return {
        "segment": segment,
        "angle": float(constraint.get("angle", _segment_angle(segment))),
        "point": point,
        "chainage": float(path_line.project(Point(point[0], point[1]))),
    }


def _build_station_layout(
    path_line: LineString,
    anchors: pd.DataFrame,
    boundary_events: list[dict[str, Any]],
    start_cut: dict[str, Any] | None,
    end_cut: dict[str, Any] | None,
    min_strip_horizontal_length: float,
) -> dict[str, Any] | None:
    start_chainage = float(start_cut["chainage"]) if start_cut is not None else 0.0
    end_chainage = float(end_cut["chainage"]) if end_cut is not None else float(path_line.length)
    if end_chainage <= start_chainage:
        return None

    interior_length = end_chainage - start_chainage
    if min_strip_horizontal_length <= 0:
        count_strip = 2
    else:
        count_strip = max(2, int(math.floor(interior_length / min_strip_horizontal_length - 1e-9)))
    station_chainages = np.linspace(start_chainage, end_chainage, count_strip + 1)[1:-1].tolist()
    if start_cut is not None:
        station_chainages = [start_chainage] + station_chainages
    if end_cut is not None:
        station_chainages = station_chainages + [end_chainage]
    station_chainages = sorted(set(float(item) for item in station_chainages if item > 0 and item < path_line.length))

    endpoint_events = []
    endpoint_events.append({
        "chainage": start_chainage,
        "angle": float(start_cut["angle"]) if start_cut is not None else _path_normal_angle(path_line, 0.0),
    })
    endpoint_events.extend({"chainage": event["chainage"], "angle": event["angle"]} for event in boundary_events)
    endpoint_events.append({
        "chainage": end_chainage,
        "angle": float(end_cut["angle"]) if end_cut is not None else _path_normal_angle(path_line, path_line.length),
    })
    endpoint_events = sorted(endpoint_events, key=lambda item: item["chainage"])

    station_events: list[dict[str, Any]] = []
    for idx, chainage in enumerate(station_chainages, start=1):
        point = path_line.interpolate(chainage)
        angle = _interpolate_station_angle(chainage, endpoint_events)
        station_events.append(
            {
                "Cx": float(point.x),
                "Cy": float(point.y),
                "Czp": 0.0,
                "type": "Station: Strip boundary",
                "grp": 0,
                "station": idx,
                "chainage": float(chainage),
                "angle": float(angle),
            }
        )

    anchor_events = []
    for _, row in anchors.iterrows():
        point = Point(float(row["Cx"]), float(row["Cy"]))
        anchor_events.append(
            {
                "Cx": float(row["Cx"]),
                "Cy": float(row["Cy"]),
                "Czp": float(row["Cz"]),
                "type": row["type"],
                "grp": int(row["grp"]),
                "station": 0,
                "chainage": float(path_line.project(point)),
            }
        )

    return {
        "anchor_events": anchor_events,
        "boundary_events": boundary_events,
        "station_events": station_events,
        "start_cut": start_cut,
        "end_cut": end_cut,
    }


def _path_normal_angle(path_line: LineString, chainage: float) -> float:
    tangent = _line_tangent(path_line, chainage)
    return _normalize_angle(math.atan2(tangent[1], tangent[0]) + math.pi / 2.0)


def _interpolate_station_angle(chainage: float, events: list[dict[str, float]]) -> float:
    if chainage <= events[0]["chainage"]:
        return events[0]["angle"]
    if chainage >= events[-1]["chainage"]:
        return events[-1]["angle"]
    for idx in range(1, len(events)):
        upper = events[idx - 1]
        lower = events[idx]
        if chainage <= lower["chainage"]:
            length = lower["chainage"] - upper["chainage"]
            if length == 0:
                return lower["angle"]
            ratio = (chainage - upper["chainage"]) / length
            diff = _normalize_angle(lower["angle"] - upper["angle"])
            return _normalize_angle(upper["angle"] + diff * ratio)
    return events[-1]["angle"]


def _build_station_cutlines(polygon: Polygon, path_line: LineString, station_info: dict[str, Any]) -> list[dict[str, Any]]:
    cutlines: list[dict[str, Any]] = []
    start_cut = station_info.get("start_cut")
    end_cut = station_info.get("end_cut")
    for event in station_info["station_events"]:
        segment = None
        protected = False
        if start_cut is not None and math.isclose(event["chainage"], start_cut["chainage"], abs_tol=1e-8):
            segment = np.array(start_cut["segment"], dtype=float)
            protected = True
        elif end_cut is not None and math.isclose(event["chainage"], end_cut["chainage"], abs_tol=1e-8):
            segment = np.array(end_cut["segment"], dtype=float)
            protected = True
        if segment is None:
            direction = np.array([math.cos(event["angle"]), math.sin(event["angle"])] , dtype=float)
            segment = _line_segment_through_point_onr(polygon, np.array([event["Cx"], event["Cy"]], dtype=float), direction)
        if segment is None:
            continue
        cutlines.append({"segment": segment, "protected": protected, "event": event})
    return cutlines


def _point_inside_polygon_tolerant(polygon: Polygon, point_xy: np.ndarray) -> bool:
    point = Point(float(point_xy[0]), float(point_xy[1]))
    repaired = _repair_polygon(polygon)
    if repaired is None:
        return False
    try:
        test_geom = repaired.buffer(1e-9)
        return bool(test_geom.contains(point) or test_geom.touches(point))
    except GEOSException:
        try:
            repaired2 = repaired.buffer(0)
            test_geom = repaired2.buffer(1e-9)
            return bool(test_geom.contains(point) or test_geom.touches(point))
        except GEOSException:
            return False


def _filter_cutlines_inside_polygon(polygon: Polygon, cutlines: list[dict[str, Any]]) -> list[dict[str, Any]]:
    filtered: list[dict[str, Any]] = []
    for item in cutlines:
        point_xy = np.array([item["event"]["Cx"], item["event"]["Cy"]], dtype=float)
        if _point_inside_polygon_tolerant(polygon, point_xy):
            filtered.append(item)
    return filtered


def _remove_intersecting_cutlines(cutlines: list[dict[str, Any]]) -> list[dict[str, Any]]:
    cutlines = sorted(cutlines, key=lambda item: item["event"]["chainage"])
    changed = True
    while changed and len(cutlines) >= 2:
        changed = False
        counts = [0 for _ in cutlines]
        for i in range(len(cutlines) - 1):
            line_i = LineString(cutlines[i]["segment"].tolist())
            for j in range(i + 1, len(cutlines)):
                line_j = LineString(cutlines[j]["segment"].tolist())
                inter = line_i.intersection(line_j)
                if inter.is_empty:
                    continue
                counts[i] += 1
                counts[j] += 1
        protected = [item["protected"] for item in cutlines]
        for idx, flag in enumerate(protected):
            if flag:
                counts[idx] = 0
        if max(counts, default=0) > 0:
            idxs = [idx for idx, value in enumerate(counts) if value == max(counts)]
            scores = []
            for idx in idxs:
                upper = 0.0 if idx == 0 else abs(_normalize_angle(cutlines[idx]["event"]["angle"] - cutlines[idx - 1]["event"]["angle"]))
                lower = 0.0 if idx == len(cutlines) - 1 else abs(_normalize_angle(cutlines[idx]["event"]["angle"] - cutlines[idx + 1]["event"]["angle"]))
                scores.append(upper + lower)
            drop_idx = idxs[int(np.argmax(scores))]
            cutlines.pop(drop_idx)
            changed = True
    return cutlines


def _split_polygon_into_strips(
    polygon: Polygon,
    path_line: LineString,
    anchors: pd.DataFrame,
    cutlines: list[dict[str, Any]],
) -> tuple[list[Polygon], list[dict[str, Any]]] | None:
    cutlines = sorted(cutlines, key=lambda item: item["event"]["chainage"])
    remaining = polygon
    strips: list[Polygon] = []
    anchor_events = [
        {
            "Cx": float(row["Cx"]),
            "Cy": float(row["Cy"]),
            "Czp": float(row["Cz"]),
            "type": row["type"],
            "grp": int(row["grp"]),
            "station": 0,
            "chainage": float(path_line.project(Point(float(row["Cx"]), float(row["Cy"])))),
        }
        for _, row in anchors.iterrows()
    ]
    upstream_point = np.array([anchor_events[0]["Cx"], anchor_events[0]["Cy"]], dtype=float)
    station_events = [item["event"].copy() for item in cutlines]
    kept_station_events: list[dict[str, Any]] = []

    for idx, cutline in enumerate(cutlines):
        parts = _split_polygon_segment(remaining, cutline["segment"], extend_ratio=1e-3)
        if len(parts) != 2:
            return None
        station_xy = np.array([cutline["event"]["Cx"], cutline["event"]["Cy"]], dtype=float)
        strip_idx = _select_upstream_polygon(parts, upstream_point, station_xy)
        strip = _repair_polygon(parts[strip_idx])
        remaining = _repair_polygon(parts[1 - strip_idx])
        if strip is None or remaining is None:
            continue
        strips.append(strip)

        intersections = _path_polygon_boundary_intersections(path_line, strip)
        if not intersections:
            return None
        if len(intersections) >= 3:
            return None
        nearest_idx = int(np.argmin([np.linalg.norm(point - station_xy) for point in intersections]))
        station_event = station_events[idx]
        station_event["Cx"] = float(intersections[nearest_idx][0])
        station_event["Cy"] = float(intersections[nearest_idx][1])
        station_event["chainage"] = float(path_line.project(Point(*intersections[nearest_idx])))
        if not kept_station_events and len(intersections) >= 2:
            other_idx = 1 - nearest_idx if len(intersections) == 2 else int(np.argmax([np.linalg.norm(point - intersections[nearest_idx]) for point in intersections]))
            anchor_events[0]["Cx"] = float(intersections[other_idx][0])
            anchor_events[0]["Cy"] = float(intersections[other_idx][1])
            anchor_events[0]["chainage"] = float(path_line.project(Point(*intersections[other_idx])))
        kept_station_events.append(station_event)
        upstream_point = np.array([station_event["Cx"], station_event["Cy"]], dtype=float)

    strips.append(remaining)
    distal_intersections = _path_polygon_boundary_intersections(path_line, remaining)
    if distal_intersections:
        if len(distal_intersections) >= 3:
            return None
        current = np.array([anchor_events[-1]["Cx"], anchor_events[-1]["Cy"]], dtype=float)
        best_idx = int(np.argmin([np.linalg.norm(point - current) for point in distal_intersections]))
        anchor_events[-1]["Cx"] = float(distal_intersections[best_idx][0])
        anchor_events[-1]["Cy"] = float(distal_intersections[best_idx][1])
        anchor_events[-1]["chainage"] = float(path_line.project(Point(*distal_intersections[best_idx])))

    if not strips:
        return None

    node_events = anchor_events + kept_station_events
    node_events.sort(key=lambda item: item["chainage"])
    return strips, node_events


def _select_upstream_polygon(parts: list[Polygon], upstream_point: np.ndarray, station_xy: np.ndarray) -> int:
    probe = upstream_point + (station_xy - upstream_point) * 1e-6
    probe_point = Point(float(probe[0]), float(probe[1]))
    contains = []
    for part in parts:
        repaired = _repair_polygon(part)
        if repaired is None:
            contains.append(False)
            continue
        try:
            test_geom = repaired.buffer(1e-9)
            contains.append(test_geom.contains(probe_point) or test_geom.touches(probe_point))
        except GEOSException:
            try:
                repaired2 = repaired.buffer(0)
                test_geom = repaired2.buffer(1e-9)
                contains.append(test_geom.contains(probe_point) or test_geom.touches(probe_point))
            except GEOSException:
                contains.append(False)
    if contains.count(True) == 1:
        return contains.index(True)
    distances = [part.distance(probe_point) for part in parts]
    return int(np.argmin(distances))


def _path_polygon_boundary_intersections(path_line: LineString, polygon: Polygon) -> list[np.ndarray]:
    inter = path_line.intersection(polygon.boundary)
    points = _extract_intersection_points(inter)
    unique_points: list[np.ndarray] = []
    for point in points:
        coords = np.array([point.x, point.y], dtype=float)
        if not any(np.allclose(coords, existing, atol=1e-8) for existing in unique_points):
            unique_points.append(coords)
    unique_points.sort(key=lambda item: path_line.project(Point(item[0], item[1])))
    return unique_points


def _build_triangle_surface_package(
    strip_polygons: list[Polygon],
    path_line: LineString,
    node_events: list[dict[str, Any]],
    anchors: pd.DataFrame,
    src: rasterio.io.DatasetReader,
) -> dict[str, Any]:
    triangle_geometries: list[Polygon] = []
    triangle_records: list[dict[str, Any]] = []
    triangle_planes: list[dict[str, Any]] = []
    strip_surface_area: list[float] = []
    strip_drop: list[float] = []

    for strip_id, strip in enumerate(strip_polygons, start=1):
        triangles = _triangulate_strip_polygon(strip)
        area3d_sum = 0.0
        z_values: list[float] = []
        for triangle in triangles:
            vertices = np.array(triangle.exterior.coords[:-1], dtype=float)[:3, :2]
            z = np.array([_sample_dem(src, x, y) for x, y in vertices], dtype=float)
            area2d = float(triangle.area)
            area3d = _triangle_area_3d(vertices, z)
            triangle_geometries.append(triangle)
            triangle_records.append({"StripID": strip_id, "Area2D": area2d, "Area3D": area3d})
            triangle_planes.append({"geometry": triangle, "xy": vertices, "z": z})
            area3d_sum += area3d
            z_values.extend(z.tolist())
        strip_surface_area.append(area3d_sum)
        strip_drop.append(max(z_values) - min(z_values) if z_values else 0.0)

    triangle_events = _triangle_intersection_events(path_line, triangle_planes)
    merged_events = _merge_path_events(node_events, triangle_events, path_line)
    strip_nodes = _finalize_strip_nodes(merged_events, anchors, triangle_planes, src)
    strip_df = _build_strip_dataframe(strip_polygons, strip_nodes, strip_surface_area, strip_drop)

    paras_strip = None
    if not strip_df.empty:
        paras_strip = (
            float(strip_df["Lhrz"].sum()),
            float(strip_df["Ahrz"].sum()),
            float(strip_df["Lall"].sum()),
            float(strip_df["Aall"].sum()),
        )

    return {
        "strip_nodes": strip_nodes,
        "strip_df": strip_df,
        "triangle_df": pd.DataFrame(triangle_records),
        "triangle_geometries": triangle_geometries,
        "paras_strip": paras_strip,
    }


def _triangulate_strip_polygon(strip: Polygon) -> list[Polygon]:
    strip = _repair_polygon(strip)
    if strip is None:
        return []

    container = strip.buffer(1e-9)
    triangles: list[Polygon] = []
    for triangle in shapely_triangulate(strip):
        if triangle.is_empty:
            continue
        if triangle.area <= 0:
            continue
        try:
            inside = container.contains(triangle.representative_point())
        except GEOSException:
            try:
                container = strip.buffer(0).buffer(1e-9)
                inside = container.contains(triangle.representative_point())
            except GEOSException:
                inside = True
        if not inside:
            continue
        triangles.append(triangle)
    return triangles


def _triangle_area_3d(vertices_xy: np.ndarray, vertices_z: np.ndarray) -> float:
    p1 = np.array([vertices_xy[0, 0], vertices_xy[0, 1], vertices_z[0]], dtype=float)
    p2 = np.array([vertices_xy[1, 0], vertices_xy[1, 1], vertices_z[1]], dtype=float)
    p3 = np.array([vertices_xy[2, 0], vertices_xy[2, 1], vertices_z[2]], dtype=float)
    return float(np.linalg.norm(np.cross(p2 - p1, p3 - p1)) / 2.0)


def _triangle_intersection_events(path_line: LineString, triangle_planes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    for plane in triangle_planes:
        xy = plane["xy"]
        for idx in range(3):
            segment = np.vstack([xy[idx], xy[(idx + 1) % 3]])
            inter = _polyline_segment_intersection(path_line, segment)
            if inter is None:
                continue
            events.append(
                {
                    "Cx": float(inter[0]),
                    "Cy": float(inter[1]),
                    "Czp": 0.0,
                    "type": "Intersection: Triangle boundary",
                    "grp": 0,
                    "station": 0,
                    "chainage": float(path_line.project(Point(inter[0], inter[1]))),
                }
            )
    return events


def _merge_path_events(base_events: list[dict[str, Any]], extra_events: list[dict[str, Any]], path_line: LineString) -> list[dict[str, Any]]:
    merged = [event.copy() for event in base_events]
    for event in extra_events:
        coords = np.array([event["Cx"], event["Cy"]], dtype=float)
        duplicated = False
        for existing in merged:
            existing_coords = np.array([existing["Cx"], existing["Cy"]], dtype=float)
            if np.allclose(coords, existing_coords, atol=1e-8):
                duplicated = True
                break
        if not duplicated:
            merged.append(event.copy())
    merged.sort(key=lambda item: item["chainage"])
    for event in merged:
        event["chainage"] = float(path_line.project(Point(event["Cx"], event["Cy"])))
    merged.sort(key=lambda item: item["chainage"])
    return merged


def _finalize_strip_nodes(
    events: list[dict[str, Any]],
    anchors: pd.DataFrame,
    triangle_planes: list[dict[str, Any]],
    src: rasterio.io.DatasetReader,
) -> pd.DataFrame:
    df = pd.DataFrame(events).copy()
    df = df.sort_values("chainage").reset_index(drop=True)
    df["Czd"] = df.apply(lambda row: _sample_dem(src, float(row["Cx"]), float(row["Cy"])), axis=1)
    df["Czp"] = _interpolate_anchor_elevations(df, anchors)
    df["Czs"] = df.apply(lambda row: _surface_z_at_node(row, df, triangle_planes), axis=1)
    df["Lhrz"] = _cumulative_length(df[["Cx", "Cy"]].to_numpy(dtype=float))
    df["Lalld"] = _cumulative_length(df[["Cx", "Cy", "Czd"]].to_numpy(dtype=float))
    df["Lallp"] = _cumulative_length(df[["Cx", "Cy", "Czp"]].to_numpy(dtype=float))
    df["Lalls"] = _cumulative_length(df[["Cx", "Cy", "Czs"]].to_numpy(dtype=float))
    df["Cz"] = df["Czd"]
    df["Lall"] = df["Lalld"]
    return df[["Cx", "Cy", "Cz", "Czd", "Czp", "Czs", "type", "grp", "station", "Lhrz", "Lall", "Lalld", "Lallp", "Lalls"]]


def _interpolate_anchor_elevations(df: pd.DataFrame, anchors: pd.DataFrame) -> pd.Series:
    anchor_mask = df["type"].astype(str).str.startswith("Anchor:")
    anchor_chainages = df.loc[anchor_mask, "chainage"].to_numpy(dtype=float)
    anchor_z = df.loc[anchor_mask, "Czp"].to_numpy(dtype=float)
    if len(anchor_chainages) == 0:
        return pd.Series(np.zeros(len(df)), index=df.index)
    values = np.interp(df["chainage"].to_numpy(dtype=float), anchor_chainages, anchor_z)
    return pd.Series(values, index=df.index)


def _surface_z_at_node(row: pd.Series, df: pd.DataFrame, triangle_planes: list[dict[str, Any]]) -> float:
    point = Point(float(row["Cx"]), float(row["Cy"]))
    idx = df.index[df.index == row.name][0]
    probe_points = [point]
    if idx > 0:
        probe_points.append(Point((float(row["Cx"]) + float(df.iloc[idx - 1]["Cx"])) / 2.0, (float(row["Cy"]) + float(df.iloc[idx - 1]["Cy"])) / 2.0))
    if idx < len(df) - 1:
        probe_points.append(Point((float(row["Cx"]) + float(df.iloc[idx + 1]["Cx"])) / 2.0, (float(row["Cy"]) + float(df.iloc[idx + 1]["Cy"])) / 2.0))
    for probe in probe_points:
        for plane in triangle_planes:
            if plane["geometry"].buffer(1e-9).contains(probe) or plane["geometry"].buffer(1e-9).touches(probe):
                return _plane_z(plane["xy"], plane["z"], float(row["Cx"]), float(row["Cy"]))
    return float(row["Czd"])


def _plane_z(xy: np.ndarray, z: np.ndarray, x: float, y: float) -> float:
    p1 = np.array([xy[0, 0], xy[0, 1], z[0]], dtype=float)
    p2 = np.array([xy[1, 0], xy[1, 1], z[1]], dtype=float)
    p3 = np.array([xy[2, 0], xy[2, 1], z[2]], dtype=float)
    normal = np.cross(p2 - p1, p3 - p1)
    if np.isclose(normal[2], 0.0):
        return float(z.mean())
    return float(-(normal[0] * (x - p1[0]) + normal[1] * (y - p1[1])) / normal[2] + p1[2])


def _cumulative_length(coords: np.ndarray) -> np.ndarray:
    if len(coords) == 0:
        return np.array([], dtype=float)
    values = [0.0]
    for idx in range(1, len(coords)):
        values.append(values[-1] + float(np.linalg.norm(coords[idx] - coords[idx - 1])))
    return np.array(values, dtype=float)


def _build_strip_dataframe(
    strip_polygons: list[Polygon],
    strip_nodes: pd.DataFrame,
    strip_surface_area: list[float],
    strip_drop: list[float],
) -> pd.DataFrame:
    station_nodes = strip_nodes[strip_nodes["type"] == "Station: Strip boundary"].copy()
    boundary_lhrz = [0.0] + station_nodes["Lhrz"].tolist() + [float(strip_nodes.iloc[-1]["Lhrz"])]
    boundary_lall = [0.0] + station_nodes["Lall"].tolist() + [float(strip_nodes.iloc[-1]["Lall"])]
    boundary_lalls = [0.0] + station_nodes["Lalls"].tolist() + [float(strip_nodes.iloc[-1]["Lalls"])]
    records = []
    for idx, strip in enumerate(strip_polygons):
        l_hrz = boundary_lhrz[idx + 1] - boundary_lhrz[idx]
        l_all = boundary_lall[idx + 1] - boundary_lall[idx]
        l_alls = boundary_lalls[idx + 1] - boundary_lalls[idx]
        a_hrz = float(strip.area)
        a_srf = float(strip_surface_area[idx]) if idx < len(strip_surface_area) else a_hrz
        w_hrz = a_hrz / l_hrz if l_hrz else 0.0
        w_all = a_srf / l_alls if l_alls else 0.0
        a_all = l_all * w_all
        records.append(
            {
                "Lhrz": l_hrz,
                "Lall": l_all,
                "Ahrz": a_hrz,
                "Aall": a_all,
                "Whrz": w_hrz,
                "Wall": w_all,
                "Drop": float(strip_drop[idx]) if idx < len(strip_drop) else 0.0,
                "PATH": 1,
            }
        )
    return pd.DataFrame(records)


# ==================== 批处理配置 ====================
_SHP_FOLDER   = Path("/Users/chenrc/Nutstore Files/我的坚果云/近期工作/[01]近期科研工作/论文-滑坡主滑线提取方法_data/wenchuan/Wenchuan_polygon_WGS84")
_DEM_FOLDER   = Path("/Users/chenrc/Nutstore Files/我的坚果云/近期工作/[01]近期科研工作/论文-滑坡主滑线提取方法_data/wenchuan/dem_folder")
_OUTPUT_ROOT  = Path("/Users/chenrc/Desktop/滑坡主轴线提取_data/ALPA_py版本测试")


def _build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="ALPA Python version: supports single-case mode aligned with the R runner and batch mode.",
    )
    parser.add_argument("--dem", type=Path, help="Input DEM raster path for single-case mode.")
    parser.add_argument("--landslides", type=Path, help="Input landslide shapefile path for single-case mode.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_OUTPUT_ROOT,
        help="Output directory. Defaults to the requested desktop test folder.",
    )
    parser.add_argument(
        "--shp-folder",
        type=Path,
        default=_SHP_FOLDER,
        help="Folder containing polygon shapefiles for batch mode.",
    )
    parser.add_argument(
        "--dem-folder",
        type=Path,
        default=_DEM_FOLDER,
        help="Folder containing DEM rasters for batch mode.",
    )
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Run all polygon shapefiles under --shp-folder. If omitted, single-case mode is used when --dem and --landslides are provided.",
    )
    return parser


def _run_single_case(dem_path: Path, landslide_path: Path, output_root: Path) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    result = ALPA(
        InDEM=dem_path,
        InLandslides=landslide_path,
        output_dir=output_root,
    )

    png_path = output_root / f"{landslide_path.stem}_paths2d.png"
    _generate_png(landslide_path, result, png_path)

    print("\n========== 单案例运行完成 ==========")
    print(f"DEM: {dem_path}")
    print(f"Landslides: {landslide_path}")
    print(f"Output dir: {output_root}")
    print(f"PNG: {png_path}")
    print("===================================")


def _run_batch_cases(shp_folder: Path, dem_folder: Path, output_root: Path) -> None:
    # ---- 创建输出目录结构 ----
    for sub in ["png", "shp", "csv", "cases"]:
        (output_root / sub).mkdir(parents=True, exist_ok=True)

    # ---- 扫描所有面要素 shp 文件 ----
    shp_files = []
    for f in shp_folder.iterdir():
        if f.suffix.lower() != '.shp':
            continue
        try:
            gdf = gpd.read_file(f)
            if gdf.geometry.dtype.name != 'geometry':
                continue
            if gdf.geometry.geom_type.iloc[0] in ('Polygon', 'MultiPolygon'):
                shp_files.append(f)
        except Exception:
            pass
    shp_files.sort(key=lambda p: p.name)
    print(f"共找到 {len(shp_files)} 个面要素 shp 文件")

    csv_records: list[dict] = []
    failed_cases: list[dict] = []
    total_start = time.time()
    temp_dem_path = output_root / "temp_dem_proj.tif"

    for shp_path in shp_files:
        name = shp_path.stem
        case_start = time.time()
        print(f"\n处理: {name}")

        dem_path = _find_matching_dem(shp_path, dem_folder)
        if dem_path is None:
            elapsed = time.time() - case_start
            msg = "No matching DEM found"
            print(f"  [{name}] 失败：{msg}（{elapsed:.1f}s）")
            failed_cases.append({"File Name": shp_path.name, "Error Message": msg})
            continue

        try:
            reproj_dem = _reproject_dem_to_match_shp(shp_path, dem_path, temp_dem_path)
        except Exception as e:
            elapsed = time.time() - case_start
            print(f"  [{name}] DEM 重投影失败：{e}（{elapsed:.1f}s）")
            failed_cases.append({"File Name": shp_path.name, "Error Message": f"DEM reproject: {e}"})
            continue

        case_output_dir = output_root / "cases" / name
        case_output_dir.mkdir(parents=True, exist_ok=True)
        try:
            result = ALPA(
                InDEM=reproj_dem,
                InLandslides=shp_path,
                output_dir=case_output_dir,
            )
        except Exception as e:
            elapsed = time.time() - case_start
            print(f"  [{name}] ALPA 失败：{e}（{elapsed:.1f}s）")
            failed_cases.append({"File Name": shp_path.name, "Error Message": f"ALPA: {e}"})
            if temp_dem_path.exists():
                temp_dem_path.unlink()
            continue

        png_path = output_root / "png" / f"{name}.png"
        try:
            _generate_png(shp_path, result, png_path)
            print(f"  [{name}] PNG 已保存")
        except Exception as e:
            print(f"  [{name}] PNG 生成失败：{e}")

        path_gdf: gpd.GeoDataFrame | None = result.get("paths")
        failed_df = result.get("failed_cases")
        total_cases = int(result.get("total_cases", 0))
        success_cases = int(result.get("success_cases", 0))
        failed_count = len(failed_df) if isinstance(failed_df, pd.DataFrame) else 0
        has_final_path = _has_valid_path_gdf(path_gdf)
        batch_case_ok = has_final_path and failed_count == 0 and success_cases == total_cases

        if not batch_case_ok:
            if not has_final_path:
                msg = "No final path geometry was generated"
            else:
                msg = f"{failed_count}/{total_cases} subcases failed to generate final paths"
            failed_cases.append({"File Name": shp_path.name, "Error Message": msg})
            if temp_dem_path.exists():
                temp_dem_path.unlink()
            elapsed = time.time() - case_start
            print(f"  [{name}] 失败：{msg}（{elapsed:.1f}s）")
            continue

        if has_final_path:
            try:
                _safe_write_gdf(path_gdf, output_root / "shp" / f"{name}_path")
                print(f"  [{name}] 路径线 shp 已保存")
            except Exception as e:
                print(f"  [{name}] 路径线 shp 保存失败：{e}")

        rec: dict = {"Name": name, "success": 1}
        if has_final_path:
            non_geom_cols = [c for c in path_gdf.columns if c != 'geometry']
            for col in non_geom_cols:
                val = path_gdf[col].iloc[0]
                if isinstance(val, float):
                    val = round(val, 4)
                rec[col] = val
        csv_records.append(rec)

        if temp_dem_path.exists():
            temp_dem_path.unlink()

        elapsed = time.time() - case_start
        print(f"  [{name}] 完成，耗时 {elapsed:.1f} 秒")

    summary_csv = output_root / "csv" / "ALPA_summary.csv"
    pd.DataFrame(csv_records).to_csv(summary_csv, index=False, encoding='utf-8-sig')
    print(f"\n汇总 CSV 已保存至: {summary_csv}")

    if failed_cases:
        fail_csv = output_root / "failed_files_record.csv"
        pd.DataFrame(failed_cases).to_csv(fail_csv, index=False, encoding='utf-8')
        print(f"失败记录已保存至: {fail_csv}")

    total_elapsed = time.time() - total_start
    print(f"\n========== 运行汇总 ==========")
    print(f"总案例数：{len(shp_files)} 个")
    print(f"成功完成：{len(csv_records)} 个")
    print(f"失败案例：{len(failed_cases)} 个")
    print(f"总耗时：{int(total_elapsed // 60)} 分 {total_elapsed % 60:.1f} 秒")
    print("==============================")


def _find_matching_dem(shp_path: Path, dem_folder: Path) -> Path | None:
    """在 dem_folder 中查找空间范围覆盖该 shp 的 DEM 文件"""
    gdf = gpd.read_file(shp_path)  # WGS84 用于范围匹配
    b = gdf.total_bounds  # [minx, miny, maxx, maxy]
    for dem_file in dem_folder.iterdir():
        if dem_file.suffix.lower() != '.tif':
            continue
        with rasterio.open(dem_file) as src:
            db = src.bounds
            if db.left <= b[0] and db.bottom <= b[1] and db.right >= b[2] and db.top >= b[3]:
                return dem_file
    return None


def _reproject_dem_to_match_shp(shp_path: Path, dem_path: Path, temp_dem_path: Path) -> Path:
    """将 DEM 重投影到与 shp 质心对应的 UTM 坐标系，保存为 temp_dem_path"""
    gdf = gpd.read_file(shp_path)
    centroid = gdf.geometry.unary_union.centroid
    zone = int((centroid.x + 180) / 6) + 1
    target_crs = f'EPSG:{"326" if centroid.y >= 0 else "327"}{zone}'

    if temp_dem_path.exists():
        temp_dem_path.unlink()
    with rasterio.open(dem_path) as src:
        t, w, h = calculate_default_transform(src.crs, target_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({'crs': target_crs, 'transform': t, 'width': w, 'height': h})
        with rasterio.open(temp_dem_path, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform, src_crs=src.crs,
                    dst_transform=t, dst_crs=target_crs,
                    resampling=Resampling.nearest
                )
    return temp_dem_path


def _generate_png(shp_path: Path, result: dict, output_png_path: Path) -> None:
    """生成平面图 PNG：滑坡轮廓 + ALPA 路径线"""
    path_gdf = result.get("paths")  # ALPA 返回字典中合并路径的 key 是 'paths'
    # 读取滑坡多边形，轮廓用原始 WGS84
    gdf_poly = gpd.read_file(shp_path)
    target_crs = path_gdf.crs if path_gdf is not None and not path_gdf.empty else gdf_poly.crs
    if gdf_poly.crs != target_crs:
        gdf_poly = gdf_poly.to_crs(target_crs)

    fig, ax = plt.subplots(figsize=(8, 8))
    gdf_poly.plot(ax=ax, facecolor='lightyellow', edgecolor='black', linewidth=1.5, label='Landslide')

    if path_gdf is not None and not path_gdf.empty:
        path_gdf.plot(ax=ax, color='blue', linewidth=2.0, label='ALPA Path')

    ax.set_aspect('equal', adjustable='datalim')
    ax.legend(fontsize=8, loc='best')
    ax.set_title(shp_path.stem, fontsize=10)
    plt.tight_layout()
    plt.savefig(output_png_path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    parser = _build_cli_parser()
    cli_args = parser.parse_args()

    if cli_args.batch:
        _run_batch_cases(cli_args.shp_folder, cli_args.dem_folder, cli_args.output_dir)
    else:
        if cli_args.dem is None or cli_args.landslides is None:
            parser.error("single-case mode requires both --dem and --landslides, or use --batch")
        _run_single_case(cli_args.dem, cli_args.landslides, cli_args.output_dir)
