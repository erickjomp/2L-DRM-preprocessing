"""Build the small regression-test case from the AMZ2 CPM 4 km run.

The original wrfout files (1055 x 863 x 60, hourly, ~2.5 GB each) are cropped
to a NXC x NYC window over the Ecuadorian Andes and to 2 days that cross a
month boundary (2018-05-31 00:00 .. 2018-06-02 00:00, 49 files). Only the
variables read by prepro2LDRM-slabvars are kept. The same window and time
steps are cut from the outputs of the original run
(scratch/2LDRM_AMZ2/scripts/01_process_slabsvars_CPM-4km__2levels500.sh),
which become the reference of the "pl50000" case. The pressure of the lowest
model level is also saved (data/p_lowest_level.dat.gz), so the tests can tell
where the ground is above an interslab level without reading NetCDF.

Every calculation in the program works on one grid column (destaggering only
uses the next point, kept in the staggered dimensions), so the cropped run
must reproduce the cropped reference.

Only needed to regenerate the case (requires access to keeling, numpy and netCDF4):

    python make_test_data.py [--out-data DIR] [--out-ref DIR]

The reference files are raw float32 in Fortran order: (nx, ny, nslabs, nt).
"""
import argparse
import gzip
from pathlib import Path

import netCDF4
import numpy as np

HERE = Path(__file__).resolve().parent

LIST_FILES = Path("/data/keeling/a/erickkc2/scratch/2LDRM_AMZ2/00_raw_data/list_files/AMZ2_CPM01.txt")
REF_SRC = Path("/data/keeling/a/erickkc2/scratch/2LDRM_AMZ2/01_inputs/CPM_2levels500")
REF_BASENAME = "AMZ2_CPM01"
NX, NY, NSLABS = 1055, 863, 2

# window (0-based, unstaggered) and size of the cropped case
I0, J0 = 116, 582
NXC, NYC = 16, 14

FIRST_FILE = "2018-05-31_00:00:00"
NFILES = 49           # 48 hourly files used (-rf 2 drops the last one)
TIME_FACTOR = 3       # -tf of the original run

WRF_VARS = ["P", "PB", "U", "V", "W", "T", "QVAPOR", "XLAT", "XLONG"]

# name -> number of levels (slabs or interslabs) in the reference outputs
REF_VARS = {"PW": NSLABS, "U": NSLABS, "V": NSLABS,
            "W": NSLABS - 1, "Q": NSLABS - 1, "PWflux": NSLABS - 1, "moistdensity": NSLABS - 1}


def crop_slices(dims):
    """Index tuple cropping the WRF dimensions (Time, [z], y, x) to the window."""
    index = []
    for dim in dims:
        if dim == "Time":
            index.append(0)
        elif dim.startswith("west_east"):
            index.append(slice(I0, I0 + NXC + dim.endswith("_stag")))
        elif dim.startswith("south_north"):
            index.append(slice(J0, J0 + NYC + dim.endswith("_stag")))
        else:  # vertical, kept whole
            index.append(slice(None))
    return tuple(index)


def crop_wrfout(src, dst):
    with netCDF4.Dataset(src) as fin, netCDF4.Dataset(dst, "w", format="NETCDF4") as fout:
        fin.set_auto_mask(False)
        sizes = {"Time": None, "DateStrLen": 19, "west_east": NXC, "west_east_stag": NXC + 1,
                 "south_north": NYC, "south_north_stag": NYC + 1,
                 "bottom_top": len(fin.dimensions["bottom_top"]),
                 "bottom_top_stag": len(fin.dimensions["bottom_top_stag"])}
        for name, size in sizes.items():
            fout.createDimension(name, size)
        fout.title = f"Cropped from {src.name}: west_east {I0}:{I0 + NXC}, south_north {J0}:{J0 + NYC}"

        times = fout.createVariable("Times", "S1", ("Time", "DateStrLen"))
        times[0, :] = fin["Times"][0, :]
        for name in WRF_VARS:
            var = fin[name]
            out = fout.createVariable(name, "f4", var.dimensions, zlib=True, complevel=9, shuffle=True)
            out.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            out[0, ...] = var[crop_slices(var.dimensions)]
        # pressure of the lowest model level, summed in float32 like the program does
        return fout["P"][0, 0] + fout["PB"][0, 0]


def write_gz(path, array):
    # fixed mtime and level keep the .gz bytes identical when regenerated with the same data
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0) as f:
        f.write(np.ascontiguousarray(array, dtype="<f4").tobytes())


def crop_reference(name, nlev, start, count, out_dir):
    data = np.memmap(REF_SRC / f"{REF_BASENAME}_{name}.dat", dtype="<f4", mode="r").reshape(-1, nlev, NY, NX)
    # C order (nt, nlev, ny, nx) is the Fortran order (nx, ny, nlev, nt) of the program outputs
    write_gz(out_dir / f"{name}.dat.gz", data[start:start + count, :, J0:J0 + NYC, I0:I0 + NXC])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out-data", type=Path, default=HERE / "data")
    parser.add_argument("--out-ref", type=Path, default=HERE / "reference" / "pl50000")
    args = parser.parse_args()
    args.out_data.mkdir(parents=True, exist_ok=True)
    args.out_ref.mkdir(parents=True, exist_ok=True)

    files = [Path(f) for f in LIST_FILES.read_text().split()]
    first = next(i for i, f in enumerate(files) if f.name.endswith(FIRST_FILE))
    p_lowest = []
    for src in files[first:first + NFILES]:
        # ':' is not allowed in Windows file names
        dst = args.out_data / (src.name.replace(":", "_") + ".nc")
        print(f"{src.name} -> {dst.name}")
        p_lowest.append(crop_wrfout(src, dst))

    # the reference run also used -rf 2, so output step n starts at file n * time factor
    if first % TIME_FACTOR:
        raise SystemExit(f"{FIRST_FILE} must start an output step of the reference outputs")
    nused = NFILES - 1
    # used by the tests to know where the ground is above an interslab level, (nused, ny, nx)
    write_gz(args.out_data / "p_lowest_level.dat.gz", np.array(p_lowest[:nused]))
    for name, nlev in REF_VARS.items():
        crop_reference(name, nlev, first // TIME_FACTOR, nused // TIME_FACTOR, args.out_ref)
        print(f"reference {name}: steps {first // TIME_FACTOR}..{(first + nused) // TIME_FACTOR - 1}")


if __name__ == "__main__":
    main()
