"""Case definition and helpers shared by conftest.py and test_regression.py."""
import gzip
import os
import shutil
import subprocess
from pathlib import Path

import netCDF4
import numpy as np

HERE = Path(__file__).resolve().parent
CASE_DIR = HERE / "data"
REF_DIR = HERE / "reference"

# Cropped grid written by make_test_data.py
NX, NY, NZ, NSLABS = 16, 14, 60, 2
TIME_FACTOR = 3
NFILES_USED = 48  # 49 files in data/, -rf 2 drops the last one

LIST_NAME = "wrf_files.txt"  # its stem is the prefix of the output files
BASENAME = Path(LIST_NAME).stem
# the same files from 01:00, used with -rf 0: the output step 2018-05-31 22:00 .. 2018-06-01 00:00
# then crosses the month boundary
LIST_SHIFTED = "wrf_files_shifted.txt"

COMMON_ARGS = [
    "-tf", str(TIME_FACTOR),
    "-nx", str(NX), "-ny", str(NY), "-nz", str(NZ), "-ns", str(NSLABS),
]

# Each case adds the options that differ. pl50000 is the original AMZ2 run
# (reference cut from its outputs); at 950 hPa the ground is above the
# interslab level in about half of the window, so pl95000 covers the NaN path.
# NetCDF is the default format; the reference is compared with the binary outputs.
CASES = {
    "pl50000": ["-pl", "50000", "-fmt", "bin"],
    "pl95000": ["-pl", "95000", "-fmt", "bin"],
    "pl50000_cf2": ["-pl", "50000", "-cf", "2", "-fmt", "bin"],
    "pl50000_etf": ["-pl", "50000", "-etf", "24", "-fmt", "bin"],  # deprecated, must be ignored
    "pl50000_nc": ["-pl", "50000"],
    "pl95000_nc": ["-pl", "95000"],
    "pl50000_cf2_nc": ["-pl", "50000", "-cf", "2"],
    "pl50000_month": ["-pl", "50000", "-fmt", "bin", "-sp", "month"],
    "pl50000_month_nc": ["-pl", "50000", "-sp", "month"],
    "pl50000_year_nc": ["-pl", "50000", "-sp", "year"],
}
PRESSURE_LEVEL = {"pl50000": 50000.0, "pl95000": 95000.0}
REFERENCE_CASES = ["pl50000", "pl95000"]

# output name -> number of levels (slabs, or interslabs for variables at the slab boundaries)
SLAB_OUTPUTS = {
    "PW": NSLABS, "U": NSLABS, "V": NSLABS,
    "W": NSLABS - 1, "Q": NSLABS - 1, "PWflux": NSLABS - 1, "moistdensity": NSLABS - 1,
}
# PP and ET have no levels: read as one level, (nt, 1, ny, nx), which is the same binary layout
SURFACE_OUTPUTS = {"PP": 1, "ET": 1}
OUTPUTS = {**SLAB_OUTPUTS, **SURFACE_OUTPUTS}
NO_PP_ET = "-noppet"
SECONDS_PER_FILE = 3600

# Same compiler and flags reproduce the reference bit for bit; the tolerance
# leaves room for other compilers, optimization levels and platforms.
RTOL = 1e-5
ATOL_FRACTION = 1e-6  # of the largest absolute value of each output


def prepare_inputs(workdir):
    """Copy the case into workdir and write the list of WRF files.

    The list has relative names: the program stores paths in fixed-length strings."""
    files = sorted(CASE_DIR.glob("wrfout_*.nc"))
    for f in files:
        shutil.copy(f, workdir / f.name)
    (workdir / LIST_NAME).write_text("".join(f"{f.name}\n" for f in files))
    (workdir / LIST_SHIFTED).write_text("".join(f"{f.name}\n" for f in files[1:]))


def run_model(exe, workdir, outname, extra_args, nprocs=1, mpiexec=None, list_name=LIST_NAME, rem_file=2):
    """Run the program in workdir, writing to workdir/outname. Returns (outdir, error).

    outname must not exist yet: the program creates it (and its parents)."""
    outdir = workdir / outname
    cmd = [str(exe), "-po", outname, "-lf", list_name, "-rf", str(rem_file)] + COMMON_ARGS + extra_args
    if nprocs > 1:
        cmd = [str(mpiexec), "-n", str(nprocs)] + cmd
    # Open MPI: the UCX transport fails on hosts without a network fabric (e.g. login
    # nodes); ob1 uses shared memory, which is all a single-node run needs
    env = dict(os.environ)
    env.setdefault("OMPI_MCA_pml", "ob1")
    # allow more processes than cores (e.g. GitHub runners): Open MPI 4 and Open MPI 5 (PRRTE)
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    env.setdefault("PRTE_MCA_rmaps_default_mapping_policy", ":oversubscribe")
    result = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True)
    outdir.mkdir(parents=True, exist_ok=True)  # for the log, also if the run failed early
    (outdir / "stdout.txt").write_text(result.stdout + result.stderr)
    expected = [] if "-sp" in extra_args else expected_files(outdir, extra_args, Path(list_name).stem)
    missing = [f.name for f in expected if not f.is_file()]
    if result.returncode != 0 or missing:
        log = (result.stdout + result.stderr)[-3000:]
        return outdir, f"program failed (exit {result.returncode}, missing: {missing})\n{log}"
    return outdir, None


def output_path(outdir, name, period=None, basename=BASENAME):
    prefix = basename if period is None else f"{basename}_{period}"
    return outdir / f"{prefix}_{name}.dat"


def nc_path(outdir, period=None, basename=BASENAME):
    return outdir / (f"{basename}.nc" if period is None else f"{basename}_{period}.nc")


def expected_files(outdir, args, basename):
    """Output files of a run without --split."""
    if "bin" in args:
        names = SLAB_OUTPUTS if NO_PP_ET in args else OUTPUTS
        return [output_path(outdir, name, basename=basename) for name in names]
    return [nc_path(outdir, basename=basename)]


def reference_path(case, name):
    """PP and ET do not depend on the pressure levels: one reference for all cases."""
    return REF_DIR / ("surface" if name in SURFACE_OUTPUTS else case) / f"{name}.dat.gz"


def to_grid(raw, name, nx=NX, ny=NY):
    """Raw float32 bytes in Fortran order (nx, ny, nlev, nt) -> array (nt, nlev, ny, nx)."""
    return np.frombuffer(raw, dtype="<f4").reshape(-1, OUTPUTS[name], ny, nx)


def read_output(outdir, name, nx=NX, ny=NY, period=None):
    return to_grid(output_path(outdir, name, period).read_bytes(), name, nx, ny)


def read_nc(path):
    """Variables of a NetCDF output as numpy arrays (NaN kept), plus the global attributes."""
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_mask(False)
        data = {name: var[...] for name, var in nc.variables.items()}
        for name in SURFACE_OUTPUTS:  # (nt, ny, nx) -> (nt, 1, ny, nx), as read_output
            if name in data:
                data[name] = data[name][:, None]
        attrs = {name: {a: var.getncattr(a) for a in var.ncattrs()} for name, var in nc.variables.items()}
        attrs["global"] = {a: nc.getncattr(a) for a in nc.ncattrs()}
        dims = {name: len(dim) for name, dim in nc.dimensions.items()}
    return data, attrs, dims


def read_wrf(var):
    """A variable of the first used input file, without the Time dimension."""
    with netCDF4.Dataset(sorted(CASE_DIR.glob("wrfout_*.nc"))[0]) as nc:
        nc.set_auto_mask(False)
        return nc[var][0]


def read_wrf_all(var, case_dir=CASE_DIR):
    """A variable of all the input files (49, with the one -rf 2 drops), (nfiles, ny, nx)."""
    out = []
    for f in sorted(case_dir.glob("wrfout_*.nc")):
        with netCDF4.Dataset(f) as nc:
            nc.set_auto_mask(False)
            out.append(nc[var][0])
    return np.array(out)


def coarsen(data, cf):
    """Block mean over cf x cf cells of (..., ny, nx)."""
    if cf == 1:
        return data
    *lead, ny, nx = data.shape
    return data.reshape(*lead, ny // cf, cf, nx // cf, cf).mean(axis=(-3, -1))


def pp_from_wrf(cf=1, case_dir=CASE_DIR):
    """Expected PP, (nt, 1, ny, nx): accumulated rain (with buckets) at the end of each step
    minus at its start, the steps going from file n * TIME_FACTOR to (n + 1) * TIME_FACTOR."""
    rain = sum(read_wrf_all(v, case_dir).astype("f8") + 100.0 * read_wrf_all(f"I_{v}", case_dir)
               for v in ["RAINNC", "RAINC"])
    bounds = coarsen(rain[::TIME_FACTOR], cf)
    return np.diff(bounds, axis=0)[:, None]


def et_from_wrf(cf=1, first_file=0):
    """Expected ET, (nt, 1, ny, nx): sum over the files of each step of LH / Lv(T2B) times
    the seconds per file, negative values set to 0 (after coarsening, as the python program)."""
    lh, t2b = read_wrf_all("LH"), read_wrf_all("T2B")
    rate = lh / (np.float32(4185.5) * (np.float32(751.78) - np.float32(0.5655) * t2b))
    rate = coarsen(rate[first_file:first_file + NFILES_USED].astype("f8"), cf)
    steps = rate.reshape(-1, TIME_FACTOR, *rate.shape[1:]).sum(axis=1) * SECONDS_PER_FILE
    return np.maximum(steps, 0)[:, None]


def read_reference(case, name):
    return to_grid(gzip.decompress(reference_path(case, name).read_bytes()), name)


def read_p_lowest_level():
    """Pressure (Pa) of the lowest model level of every used file, (nfiles, ny, nx)."""
    raw = gzip.decompress((CASE_DIR / "p_lowest_level.dat.gz").read_bytes())
    return np.frombuffer(raw, dtype="<f4").reshape(NFILES_USED, NY, NX)


def save_reference(output, reference):
    """Write output as the new reference. Returns False if the content was unchanged."""
    data = output.read_bytes()
    if reference.is_file() and gzip.decompress(reference.read_bytes()) == data:
        return False
    reference.parent.mkdir(parents=True, exist_ok=True)
    # fixed mtime and level keep the .gz bytes identical when regenerated with the same data
    with open(reference, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0) as f:
        f.write(data)
    return True


def compare_grid(new, ref):
    """Compare two output arrays. Returns a list of error messages (empty if equal)."""
    if new.shape != ref.shape:
        return [f"shape differs: {new.shape} vs {ref.shape} in reference"]
    nan_new, nan_ref = np.isnan(new), np.isnan(ref)
    errors = []
    if (nan_new != nan_ref).any():
        errors.append(f"NaN in {(nan_new != nan_ref).sum()} different places "
                      f"({nan_new.sum()} NaN vs {nan_ref.sum()} in reference)")
    valid = ~(nan_new | nan_ref)
    a, b = new[valid].astype("f8"), ref[valid].astype("f8")
    atol = ATOL_FRACTION * (np.abs(b).max() if b.size else 0.0)
    bad = np.abs(a - b) > atol + RTOL * np.abs(b)
    if bad.any():
        errors.append(f"{bad.sum()} of {a.size} values differ "
                      f"(max abs diff {np.abs(a - b).max():.3g}, tolerance {atol:.3g} + {RTOL:g} * |ref|)")
    return errors
