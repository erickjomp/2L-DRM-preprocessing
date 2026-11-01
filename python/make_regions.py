"""Source regions and target region of 2L-DRM from shapefiles, on the grid of the WRF files.

    python make_regions.py --source_regions basins.shp [--field ID] --target_region 1 3 \
        --wrf wrfout_d01_2018-05-01_00:00:00 --prefix inputs/amazon [--coarse_factor 2] [--keep_ocean]
    python make_regions.py --source_regions countries.shp --field ID --target_region la_plata.shp ...

Writes, as 2L-DRM reads them (float32, Fortran order (nx, ny), i.e. numpy (ny, nx)):
- <prefix>_regions.dat: map of the source regions. The region of a cell is the --field value
  of the polygon of --source_regions that contains its center (1 for every polygon without
  --field), and nregions = (largest region) + 1 for the cells outside the polygons and,
  without --keep_ocean, over the ocean (LANDMASK and LAKEMASK both 0).
- <prefix>_il.dat, <prefix>_jl.dat: x and y indices of the target-region cells, 1-based as the
  Fortran indices of 2L-DRM. --target_region is either region numbers of --source_regions,
  or a shapefile, independent of the source regions: the cells whose center is inside any
  of its polygons (without the ocean, unless --keep_ocean).
It prints nx, ny, nregions and the number of target cells (domsize), which 2L-DRM needs.

The grid is that of prepro2LDRM-slabvars with the same --coarse_factor: XLAT and XLONG
averaged over blocks of coarse_factor x coarse_factor cells (the last rows and columns that
do not fill a block are dropped); a coarse cell is land if at least half of its cells are.

Needs numpy, netCDF4, geopandas and regionmask.
"""
import argparse
import sys
from pathlib import Path

import geopandas as gpd
import netCDF4
import numpy as np
import regionmask


def coarsen(field, cf):
    """block means of cf x cf cells of a (ny, nx) field, as prepro2LDRM-slabvars"""
    ny, nx = (s // cf for s in field.shape)
    return field[:ny * cf, :nx * cf].reshape(ny, cf, nx, cf).mean(axis=(1, 3), dtype=np.float64)


def read_grid(wrf_file, cf, keep_ocean):
    """lat, lon and land mask (or None with keep_ocean) of the coarse grid"""
    with netCDF4.Dataset(wrf_file) as nc:
        lat = coarsen(np.ma.filled(nc["XLAT"][0], np.nan), cf)
        lon = coarsen(np.ma.filled(nc["XLONG"][0], np.nan), cf)
        if keep_ocean:
            return lat, lon, None
        try:
            land = (np.ma.filled(nc["LANDMASK"][0], 0) + np.ma.filled(nc["LAKEMASK"][0], 0)) > 0
        except IndexError:
            sys.exit(f"{wrf_file}: LANDMASK or LAKEMASK not found (use --keep_ocean to keep the ocean cells)")
    return lat, lon, coarsen(land.astype(np.float64), cf) >= 0.5


def read_shapes(shapefile, field=None):
    """polygons of a shapefile in longitude-latitude, one (multi)polygon per region: the
    field values (integers >= 1), or 1 for all without field"""
    shapes = gpd.read_file(shapefile)
    if field is None:
        shapes["region"] = 1
    else:
        if field not in shapes.columns:
            sys.exit(f"{shapefile}: no field {field!r} (fields: {', '.join(c for c in shapes.columns if c != 'geometry')})")
        values = shapes[field].to_numpy()
        if not (np.issubdtype(values.dtype, np.number) and np.all(values == np.round(values)) and values.min() >= 1):
            sys.exit(f"{shapefile}: the {field!r} values must be integers >= 1 (2L-DRM numbers the regions from 1)")
        shapes["region"] = values.astype(int)
    if shapes.crs is not None and not shapes.crs.is_geographic:
        shapes = shapes.to_crs("EPSG:4326")
    return shapes[["region", "geometry"]].dissolve(by="region").reset_index()


def mask_cells(shapes, lat, lon):
    """region of each cell (by its center), NaN outside the polygons"""
    return np.asarray(regionmask.mask_geopandas(shapes, lon, lat, numbers="region"))


def region_map(shapefile, field, lat, lon, land):
    """regions of the cells (integer array) and nregions"""
    shapes = read_shapes(shapefile, field)
    mask = mask_cells(shapes, lat, lon)
    nregions = int(shapes["region"].max()) + 1
    regions = np.where(np.isfinite(mask), mask, nregions)
    if land is not None:
        regions[~land] = nregions
    return regions.astype(int), nregions


def target_cells(target_region, regions, nregions, lat, lon, land):
    """target cells (boolean array) from region numbers or a shapefile, and a description"""
    is_file = len(target_region) == 1 and Path(target_region[0]).exists()
    try:
        ids = None if is_file else [int(v) for v in target_region]
    except ValueError:
        sys.exit(f"--target_region: {' '.join(target_region)} is neither region numbers nor an existing shapefile")
    if ids is not None:
        missing = [i for i in ids if not np.any(regions == i) or i >= nregions]
        if missing:
            sys.exit(f"--target_region: regions {missing} have no cells in the source regions (1 to {nregions - 1})")
        return np.isin(regions, ids), f"source regions {' '.join(map(str, ids))}"
    # all its polygons as one region
    target = np.isfinite(mask_cells(read_shapes(target_region[0]), lat, lon))
    if land is not None:
        target &= land
    return target, f"cells inside {target_region[0]}"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--source_regions", "-sr", required=True, type=Path,
                        help="shapefile with the polygons of the source regions")
    parser.add_argument("--field", "-f", help="field of --source_regions with the region of each polygon "
                                              "(integers >= 1); without it every polygon is region 1")
    parser.add_argument("--target_region", "-tr", required=True, nargs="+", metavar="IDS_OR_SHAPEFILE",
                        help="target region: region numbers of --source_regions (e.g. 1 3), or a shapefile")
    parser.add_argument("--wrf", required=True, type=Path,
                        help="a wrfout file of the run (for XLAT, XLONG, LANDMASK, LAKEMASK)")
    parser.add_argument("--prefix", required=True, help="prefix of the output files, e.g. inputs/amazon")
    parser.add_argument("--coarse_factor", "-cf", type=int, default=1, help="as in prepro2LDRM-slabvars (default 1)")
    parser.add_argument("--keep_ocean", action="store_true", help="keep the ocean cells in the regions and the target")
    args = parser.parse_args(argv)

    lat, lon, land = read_grid(args.wrf, args.coarse_factor, args.keep_ocean)
    regions, nregions = region_map(args.source_regions, args.field, lat, lon, land)
    target, target_desc = target_cells(args.target_region, regions, nregions, lat, lon, land)
    if not target.any():
        sys.exit(f"No target cells ({target_desc}): check --target_region and the WRF grid")
    jl, il = np.nonzero(target)

    prefix = Path(args.prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    regions.astype(np.float32).tofile(f"{prefix}_regions.dat")
    (il + 1).astype(np.float32).tofile(f"{prefix}_il.dat")
    (jl + 1).astype(np.float32).tofile(f"{prefix}_jl.dat")

    ny, nx = regions.shape
    print(f"Written {prefix}_regions.dat, {prefix}_il.dat, {prefix}_jl.dat")
    print(f"Target region: {target_desc}")
    print(f"nx = {nx}, ny = {ny}, nregions = {nregions}, domsize (target cells) = {target.sum()}")


if __name__ == "__main__":
    main()
