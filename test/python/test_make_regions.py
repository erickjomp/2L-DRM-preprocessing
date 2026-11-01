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
    """a WRF file with LANDMASK and LAKEMASK, and a shapefile with regions 1 and 3"""
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
    return {"dir": d, "wrf": wrf, "shp": shp, "boxes": {1: box1, 3: box3}}


def expected_regions(case, cf, keep_ocean):
    with netCDF4.Dataset(case["wrf"]) as nc:
        lat = make_regions.coarsen(np.asarray(nc["XLAT"][0]), cf)
        lon = make_regions.coarsen(np.asarray(nc["XLONG"][0]), cf)
        land = (np.asarray(nc["LANDMASK"][0]) + np.asarray(nc["LAKEMASK"][0])) > 0
    regions = np.full(lat.shape, 4)
    for value, box in case["boxes"].items():
        regions[shapely.contains_xy(box, lon, lat)] = value
    if not keep_ocean:
        regions[make_regions.coarsen(land.astype(float), cf) < 0.5] = 4
    return regions


def run(case, name, *args):
    prefix = case["dir"] / name
    make_regions.main(["--shapefile", str(case["shp"]), "--wrf", str(case["wrf"]), "--prefix", str(prefix),
                       "--field", "ID", *args])
    regions = np.fromfile(f"{prefix}_regions.dat", np.float32)
    il = np.fromfile(f"{prefix}_il.dat", np.float32).astype(int)
    jl = np.fromfile(f"{prefix}_jl.dat", np.float32).astype(int)
    return regions, il, jl


@pytest.mark.parametrize("cf", [1, 2])
@pytest.mark.parametrize("keep_ocean", [False, True])
def test_regions(case, cf, keep_ocean):
    """region map as the point-in-polygon test; target cells: all the regions, 1-based"""
    regions, il, jl = run(case, f"cf{cf}_{keep_ocean}", "-cf", str(cf), *(["--keep_ocean"] if keep_ocean else []))
    expected = expected_regions(case, cf, keep_ocean)
    assert set(np.unique(expected)) == {1, 3, 4}  # the case has cells of every region
    assert np.array_equal(regions.reshape(expected.shape), expected)
    assert il.min() >= 1 and il.max() <= expected.shape[1] and jl.min() >= 1 and jl.max() <= expected.shape[0]
    assert np.array_equal(np.sort(expected[jl - 1, il - 1]), np.sort(expected[expected < 4]))


def test_target(case):
    """--target selects the cells of the given regions"""
    _, il, jl = run(case, "target3", "--target", "3")
    expected = expected_regions(case, 1, False)
    assert len(il) == np.sum(expected == 3) and np.all(expected[jl - 1, il - 1] == 3)


def test_without_field(case):
    """without --field every polygon is region 1, and nregions is 2"""
    prefix = case["dir"] / "nofield"
    make_regions.main(["--shapefile", str(case["shp"]), "--wrf", str(case["wrf"]), "--prefix", str(prefix)])
    regions = np.fromfile(f"{prefix}_regions.dat", np.float32)
    expected = expected_regions(case, 1, False).ravel()
    assert np.array_equal(regions, np.where(expected == 4, 2, 1))
