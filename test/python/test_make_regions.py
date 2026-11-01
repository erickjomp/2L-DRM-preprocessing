"""Tests of python/make_regions.py on the grid of the regression test case.

    pytest test/python

Needs numpy, netCDF4, geopandas and regionmask (skipped otherwise). The regions are checked
against shapely's point-in-polygon test of the cell centers.
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

gpd = pytest.importorskip("geopandas")
pytest.importorskip("regionmask")
import netCDF4  # noqa: E402
import shapely  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
import make_regions  # noqa: E402

WRF_FILE = ROOT / "test/regression/data/wrfout_d01_2018-05-31_00_00_00.nc"
LAND_FROM_X = 6  # land for x >= 6 (west_east), ocean to the west


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    """a WRF file with LANDMASK and LAKEMASK, a shapefile with source regions 1 and 3, and a
    target polygon across both regions and the coast"""
    d = tmp_path_factory.mktemp("regions")
    wrf = d / "wrfout.nc"
    shutil.copy(WRF_FILE, wrf)
    with netCDF4.Dataset(wrf, "a") as nc:
        lat, lon = nc["XLAT"][0], nc["XLONG"][0]
        land = np.zeros((1,) + lat.shape, np.float32)
        land[..., LAND_FROM_X:] = 1
        lake = np.zeros_like(land)
        lake[0, 0, LAND_FROM_X - 1] = 1  # one lake cell on the coast: kept
        for name, data in [("LANDMASK", land), ("LAKEMASK", lake)]:
            nc.createVariable(name, "f4", ("Time", "south_north", "west_east"))[:] = data
    # region 1: west half of the domain, region 3: north-east quarter, between the cell centers
    lon_mid = float(np.mean(lon[:, 7:9]))
    lat_mid = float(np.mean(lat[6:8, :]))
    box1 = shapely.box(float(lon.min()) - 1, float(lat.min()) - 1, lon_mid, float(lat.max()) + 1)
    box3 = shapely.box(lon_mid, lat_mid, float(lon.max()) + 1, float(lat.max()) + 1)
    shp = d / "regions.shp"
    gpd.GeoDataFrame({"ID": [1, 3]}, geometry=[box1, box3], crs="EPSG:4326").to_file(shp)
    # target: two polygons (one cell of overlap would not matter) in the center of the domain
    target_polys = [shapely.box(float(np.mean(lon[:, 2:4])), float(np.mean(lat[3:5, :])),
                                float(np.mean(lon[:, 11:13])), float(np.mean(lat[9:11, :]))),
                    shapely.box(float(np.mean(lon[:, 13:15])), float(lat.min()) - 1,
                                float(lon.max()) + 1, float(np.mean(lat[1:3, :])))]
    target_shp = d / "target.shp"
    gpd.GeoDataFrame({"name": ["a", "b"]}, geometry=target_polys, crs="EPSG:4326").to_file(target_shp)
    return {"dir": d, "wrf": wrf, "shp": shp, "boxes": {1: box1, 3: box3},
            "target_shp": target_shp, "target_polys": target_polys}


def grid(case, cf):
    with netCDF4.Dataset(case["wrf"]) as nc:
        lat = make_regions.coarsen(np.asarray(nc["XLAT"][0]), cf)
        lon = make_regions.coarsen(np.asarray(nc["XLONG"][0]), cf)
        land = (np.asarray(nc["LANDMASK"][0]) + np.asarray(nc["LAKEMASK"][0])) > 0
    return lat, lon, make_regions.coarsen(land.astype(float), cf) >= 0.5


def expected_regions(case, cf, keep_ocean):
    lat, lon, land = grid(case, cf)
    regions = np.full(lat.shape, 4)
    for value, box in case["boxes"].items():
        regions[shapely.contains_xy(box, lon, lat)] = value
    if not keep_ocean:
        regions[~land] = 4
    return regions


def run(case, name, target, *args, field=("--field", "ID")):
    prefix = case["dir"] / name
    make_regions.main(["--source_regions", str(case["shp"]), *field, "--target_region", *target,
                       "--wrf", str(case["wrf"]), "--prefix", str(prefix), "--format", "bin", *args])
    regions = np.fromfile(f"{prefix}_regions.dat", np.float32)
    il = np.fromfile(f"{prefix}_il.dat", np.float32).astype(int)
    jl = np.fromfile(f"{prefix}_jl.dat", np.float32).astype(int)
    # the target grid has the cells of il, jl
    cf = int(args[args.index("-cf") + 1]) if "-cf" in args else 1
    nx = 16 // cf
    target = np.fromfile(f"{prefix}_target.dat", np.float32)
    assert target.size == regions.size and set(np.unique(target)) <= {0, 1}
    assert np.array_equal(np.flatnonzero(target), (jl - 1) * nx + il - 1)
    return regions, il, jl


def target_mask(shape, il, jl):
    """target cells from 1-based il, jl"""
    assert il.min() >= 1 and il.max() <= shape[1] and jl.min() >= 1 and jl.max() <= shape[0]
    mask = np.zeros(shape, bool)
    mask[jl - 1, il - 1] = True
    assert mask.sum() == len(il)  # no repeated cells
    return mask


@pytest.mark.parametrize("cf", [1, 2])
@pytest.mark.parametrize("keep_ocean", [False, True])
def test_regions(case, cf, keep_ocean):
    """region map as the point-in-polygon test; target: source regions 1 and 3, 1-based"""
    regions, il, jl = run(case, f"cf{cf}_{keep_ocean}", ["1", "3"], "-cf", str(cf),
                          *(["--keep_ocean"] if keep_ocean else []))
    expected = expected_regions(case, cf, keep_ocean)
    assert set(np.unique(expected)) == {1, 3, 4}  # the case has cells of every region
    assert np.array_equal(regions.reshape(expected.shape), expected)
    assert np.array_equal(target_mask(expected.shape, il, jl), expected < 4)


def test_target_ids(case):
    """--target_region with some of the source regions"""
    _, il, jl = run(case, "target3", ["3"])
    expected = expected_regions(case, 1, False)
    assert np.array_equal(target_mask(expected.shape, il, jl), expected == 3)


@pytest.mark.parametrize("keep_ocean", [False, True])
def test_target_shapefile(case, keep_ocean):
    """--target_region with a shapefile independent of the source regions: the cells inside
    any of its polygons (land only without --keep_ocean)"""
    regions, il, jl = run(case, f"target_shp_{keep_ocean}", [str(case["target_shp"])],
                          *(["--keep_ocean"] if keep_ocean else []))
    lat, lon, land = grid(case, 1)
    expected = np.zeros(lat.shape, bool)
    for poly in case["target_polys"]:
        expected |= shapely.contains_xy(poly, lon, lat)
    if not keep_ocean:
        expected &= land
    target = target_mask(lat.shape, il, jl)
    assert np.array_equal(target, expected)
    # it crosses both source regions and cells outside them (land, and ocean with --keep_ocean)
    assert set(np.unique(regions.reshape(lat.shape)[target])) == {1, 3, 4}
    assert np.any(target & ~land) == keep_ocean
    assert np.array_equal(regions.reshape(lat.shape), expected_regions(case, 1, keep_ocean))


def test_without_field(case):
    """without --field every polygon is region 1, and nregions is 2"""
    regions, _, _ = run(case, "nofield", ["1"], field=())
    expected = expected_regions(case, 1, False).ravel()
    assert np.array_equal(regions, np.where(expected == 4, 2, 1))


@pytest.mark.parametrize("target", [["2"], ["5"], ["not_a_file.shp"], ["1", "x"]])
def test_bad_target(case, target):
    """region numbers without cells (2) or not source regions (5), no such file, or both"""
    with pytest.raises(SystemExit) as error:
        run(case, "bad", target)
    assert "--target_region" in str(error.value)


@pytest.mark.parametrize("target", [["1", "3"], "shapefile"])
def test_netcdf(case, target):
    """--format nc (default): same regions and target as the binary files, with XLAT, XLONG"""
    target = [str(case["target_shp"])] if target == "shapefile" else target
    regions, il, jl = run(case, "bin_for_nc", target, "-cf", "2")
    prefix = case["dir"] / "nc"
    make_regions.main(["-sr", str(case["shp"]), "-f", "ID", "-tr", *target, "--wrf", str(case["wrf"]),
                       "--prefix", str(prefix), "-cf", "2"])
    lat, lon, _ = grid(case, 2)
    with netCDF4.Dataset(f"{prefix}_regions.nc") as nc:
        assert np.array_equal(nc["regions"][:], regions.reshape(lat.shape))
        assert np.array_equal(nc["target"][:], target_mask(lat.shape, il, jl))
        assert "target_i" not in nc.variables  # the mask only
        assert np.allclose(nc["XLAT"][:], lat) and np.allclose(nc["XLONG"][:], lon)
        assert nc.nregions == 4 and nc.domsize == len(il) and nc.coarse_factor == 2
