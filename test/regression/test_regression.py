"""Regression tests: run prepro2LDRM-slabvars on the cropped AMZ2 CPM 4 km case.

    pytest test/regression --exe PATH [--mpiexec PATH] [--update-reference]

The pl50000 reference was cut from the outputs of the original AMZ2 run, and
the PP and ET_LH (ET with --et_from LH) references from the outputs of the former
python programs; ET (from QFX) is checked against numpy and was produced by this program;
the pl95000 reference was produced by an earlier version of this program.
A failure means results changed, not that they are wrong. After an intended
change in results, rerun with --update-reference and commit the new references.
"""
import netCDF4
import numpy as np
import pytest

from regression import (BASENAME, LIST_NAME, LIST_SHIFTED, NFILES_USED, NO_PP_ET, NX, NY, OUTPUTS, PRESSURE_LEVEL,
                        REFERENCE_CASES, SLAB_OUTPUTS, SURFACE_OUTPUTS, TIME_FACTOR, compare_grid,
                        et_from_wrf, nc_path, output_path, pp_from_wrf, prepare_inputs, read_nc, read_output,
                        read_p_lowest_level, read_reference, read_wrf, reference_path, run_model,
                        save_reference)

NPROCS = 3  # does not divide the 16 output steps, so processes get different numbers of chunks


@pytest.mark.parametrize("output", OUTPUTS)
@pytest.mark.parametrize("case", REFERENCE_CASES)
def test_matches_reference(run_case, update_reference, case, output):
    outdir = run_case(case)
    if update_reference:
        changed = save_reference(output_path(outdir, output), reference_path(case, output))
        pytest.skip(f"reference {'updated' if changed else 'unchanged'}: {case}/{output}")

    errors = compare_grid(read_output(outdir, output), read_reference(case, output))
    assert not errors, "\n".join(errors)


@pytest.mark.parametrize("case", REFERENCE_CASES)
def test_same_results_any_process_count(run_case, case):
    """Each chunk of files is processed independently, so 1 and NPROCS MPI processes
    must give bit-identical outputs; a difference points to an error when the
    per-process results are merged."""
    serial, parallel = run_case(case, 1), run_case(case, NPROCS)
    differing = [name for name in OUTPUTS
                 if output_path(serial, name).read_bytes() != output_path(parallel, name).read_bytes()]
    assert not differing, f"outputs differ between 1 and {NPROCS} processes: {differing}"


@pytest.mark.parametrize("nprocs", [1, NPROCS])
def test_temporary_files_deleted(run_case, nprocs):
    leftover = sorted(f.name for f in run_case("pl50000", nprocs).glob("tempdata_*"))
    assert not leftover, f"temporary files not deleted: {leftover}"


@pytest.mark.parametrize("case", REFERENCE_CASES)
def test_nan_where_ground_above_level(run_case, case):
    """Slab 1 and the interslab variables are NaN exactly where the lowest model
    level is above the interslab pressure in any hour of the output step."""
    outdir = run_case(case)
    below = read_p_lowest_level() < PRESSURE_LEVEL[case]  # (hours, ny, nx)
    expected = below.reshape(-1, TIME_FACTOR, NY, NX).any(axis=1)
    for name in SLAB_OUTPUTS:
        data = read_output(outdir, name)
        assert (np.isnan(data[:, 0]) == expected).all(), \
            f"{name}: NaN in {np.isnan(data[:, 0]).sum()} cells, expected {expected.sum()}"
        # the upper slabs always have data
        assert not np.isnan(data[:, 1:]).any(), f"{name}: NaN above the first level"


@pytest.mark.parametrize("output", SLAB_OUTPUTS)
def test_coarse_factor_2(run_case, output):
    """-cf 2 averages 2 x 2 cells after the slab calculations, so it must equal the
    2 x 2 block mean of the full-resolution outputs."""
    full = read_output(run_case("pl50000"), output).astype("f8")
    blocks = full.reshape(*full.shape[:2], NY // 2, 2, NX // 2, 2).mean(axis=(3, 5))
    errors = compare_grid(read_output(run_case("pl50000_cf2"), output, NX // 2, NY // 2), blocks)
    assert not errors, "\n".join(errors)


def test_extra_time_factor_ignored(run_case):
    """-etf (deprecated) is still accepted, with a warning, and does not change the outputs."""
    plain, with_etf = run_case("pl50000"), run_case("pl50000_etf")
    differing = [name for name in OUTPUTS
                 if output_path(plain, name).read_bytes() != output_path(with_etf, name).read_bytes()]
    assert not differing, f"outputs differ with -etf: {differing}"
    assert "deprecated" in (with_etf / "stdout.txt").read_text()
    assert not list(with_etf.glob("*24.dat")), "daily outputs written"


# NetCDF output (the default) and --split

NSTEPS = NFILES_USED // TIME_FACTOR
STEPS_PER_MONTH = {"2018-05": 8, "2018-06": 8}  # 2018-05-31 and 2018-06-01


@pytest.mark.parametrize("case", REFERENCE_CASES)
def test_netcdf_same_as_binary(run_case, case):
    """Every variable of the NetCDF output equals the binary output, bit for bit."""
    data, _, _ = read_nc(nc_path(run_case(f"{case}_nc")))
    binary = run_case(case)
    for name in OUTPUTS:
        assert np.array_equal(data[name], read_output(binary, name), equal_nan=True), name


def test_netcdf_metadata(run_case):
    data, attrs, dims = read_nc(nc_path(run_case("pl95000_nc")))
    assert dims == {"west_east": NX, "south_north": NY, "slab": 2, "interslab": 1, "time": NSTEPS}
    assert attrs["time"]["units"] == "hours since 2018-05-31 00:00:00"
    assert np.array_equal(data["time"], np.arange(NSTEPS) * TIME_FACTOR)
    assert np.array_equal(data["slab"], [1, 2])
    assert np.array_equal(data["interslab"], [95000.0])
    assert attrs["slab"]["bounds_description"] == "slab 1: surface to 95000 Pa; slab 2: 95000 Pa to model top"
    assert np.array_equal(data["XLAT"], read_wrf("XLAT"))
    assert np.array_equal(data["XLONG"], read_wrf("XLONG"))
    for name in OUTPUTS:
        assert attrs[name]["units"] and attrs[name]["long_name"] and attrs[name]["description"], name
        assert attrs[name]["coordinates"] == "XLONG XLAT", name
    with netCDF4.Dataset(nc_path(run_case("pl95000_nc"))) as nc:
        for name in SURFACE_OUTPUTS:
            assert nc[name].dimensions == ("time", "south_north", "west_east"), name
    assert attrs["global"]["time_factor"] == TIME_FACTOR
    assert attrs["global"]["split"] == "none"
    assert attrs["global"]["et_from"] == "QFX"
    assert "QFX" in attrs["ET"]["description"]
    assert np.array_equal(attrs["global"]["pressure_levels_Pa"], 95000.0)


def test_netcdf_coarse_lat_lon(run_case):
    """With -cf 2, XLAT and XLONG are the 2 x 2 block means, like the variables."""
    data, _, dims = read_nc(nc_path(run_case("pl50000_cf2_nc")))
    assert (dims["west_east"], dims["south_north"]) == (NX // 2, NY // 2)
    for var in ["XLAT", "XLONG"]:
        blocks = read_wrf(var).astype("f8").reshape(NY // 2, 2, NX // 2, 2).mean(axis=(1, 3))
        assert not compare_grid(data[var][None, None], blocks[None, None]), var
    binary = run_case("pl50000_cf2")
    for name in OUTPUTS:
        assert np.array_equal(data[name], read_output(binary, name, NX // 2, NY // 2), equal_nan=True), name


def test_split_month_binary(run_case):
    """--split month: one set of files per month, which together are the unsplit outputs."""
    split, whole = run_case("pl50000_month"), run_case("pl50000")
    assert sorted(f.name for f in split.glob("*.dat")) == sorted(
        output_path(split, name, period).name for name in OUTPUTS for period in STEPS_PER_MONTH)
    for name in OUTPUTS:
        parts = [read_output(split, name, period=period) for period in STEPS_PER_MONTH]
        assert [len(p) for p in parts] == list(STEPS_PER_MONTH.values()), name
        assert np.array_equal(np.concatenate(parts), read_output(whole, name), equal_nan=True), name


def test_split_month_netcdf(run_case):
    split = run_case("pl50000_month_nc")
    whole, _, _ = read_nc(nc_path(run_case("pl50000_nc")))
    assert sorted(f.name for f in split.glob("*.nc")) == [nc_path(split, p).name for p in STEPS_PER_MONTH]
    parts = [read_nc(nc_path(split, period)) for period in STEPS_PER_MONTH]
    for (data, attrs, _), nsteps in zip(parts, STEPS_PER_MONTH.values()):
        assert len(data["time"]) == nsteps
        # same time origin in every file, so the times continue
        assert attrs["time"]["units"] == "hours since 2018-05-31 00:00:00"
        assert attrs["global"]["split"] == "month"
    for name in ["time"] + list(OUTPUTS):
        joined = np.concatenate([data[name] for data, _, _ in parts])
        assert np.array_equal(joined, whole[name], equal_nan=True), name


def test_split_year_netcdf(run_case):
    """All the steps are in 2018, so --split year gives one file equal to the unsplit one."""
    split = run_case("pl50000_year_nc")
    assert [f.name for f in split.glob("*.nc")] == [f"{BASENAME}_2018.nc"]
    data, _, _ = read_nc(nc_path(split, "2018"))
    whole, _, _ = read_nc(nc_path(run_case("pl50000_nc")))
    for name in ["time"] + list(OUTPUTS):
        assert np.array_equal(data[name], whole[name], equal_nan=True), name


def test_split_step_across_boundary_fails(exe, inputs):
    """Starting at 01:00, the output step 2018-05-31 22:00 .. 2018-06-01 00:00 is in two
    months: --split month must stop with an error, without --split it runs."""
    _, error = run_model(exe, inputs, "shifted_month", ["-pl", "50000", "-sp", "month", NO_PP_ET],
                         list_name=LIST_SHIFTED, rem_file=0)
    assert error is not None and "is not within one month" in error
    _, error = run_model(exe, inputs, "shifted_none", ["-pl", "50000", NO_PP_ET],
                         list_name=LIST_SHIFTED, rem_file=0)
    assert error is None, error


def test_invalid_option_value_fails(exe, inputs):
    _, error = run_model(exe, inputs, "bad_format", ["-pl", "50000", "-fmt", "grib"])
    assert error is not None and "exit 0" not in error


# PP and ET

@pytest.mark.parametrize("cf", [1, 2])
def test_pp_et_from_wrf_fields(run_case, cf):
    """PP and ET recalculated with numpy from the WRF fields of the case (independent of
    the reference): PP from the accumulated rain with buckets, ET from LH and T2B."""
    outdir = run_case("pl50000" if cf == 1 else "pl50000_cf2")
    nx, ny = NX // cf, NY // cf
    for name, expected in [("PP", pp_from_wrf(cf)), ("ET", et_from_wrf(cf))]:
        errors = compare_grid(read_output(outdir, name, nx, ny), expected)
        assert not errors, f"{name}: " + "\n".join(errors)


def test_et_from_lh(run_case, update_reference):
    """--et_from LH: ET as the former python program (reference cut from its outputs), and
    2-3 % more than from QFX, since WRF's LH is QFX times a constant 2.5e6 J kg-1."""
    outdir = run_case("pl50000_etLH")
    if update_reference:
        changed = save_reference(output_path(outdir, "ET"), reference_path("pl50000", "ET_LH"))
        pytest.skip(f"reference {'updated' if changed else 'unchanged'}: ET_LH")
    et_lh = read_output(outdir, "ET")
    errors = compare_grid(et_lh, read_reference("pl50000", "ET_LH"))
    assert not errors, "\n".join(errors)
    errors = compare_grid(et_lh, et_from_wrf(et_from="LH"))
    assert not errors, "\n".join(errors)
    et_qfx = read_output(run_case("pl50000"), "ET")
    wet = et_qfx > 0.05 * et_qfx.max()
    ratio = et_lh[wet] / et_qfx[wet]
    assert 1.0 < np.median(ratio) < 1.04, np.median(ratio)


def test_remfile_1(exe, inputs, run_case):
    """-rf 1 drops the first file: the steps are the same instants for PP (the dropped file
    is the start of the first step), and ET sums the files from the second one."""
    outdir, error = run_model(exe, inputs, "rf1", ["-pl", "50000", "-fmt", "bin"], rem_file=1)
    assert error is None, error
    assert np.array_equal(read_output(outdir, "PP"), read_output(run_case("pl50000"), "PP"))
    errors = compare_grid(read_output(outdir, "ET"), et_from_wrf(first_file=1))
    assert not errors, "\n".join(errors)


def test_list_comments_blank_lines(exe, inputs, run_case):
    """The list of files can have # comments, blank lines, leading blanks and Windows line
    endings: same outputs as the plain list."""
    files = (inputs / LIST_NAME).read_text().split()
    lines = ["# WRF files of the test case", ""] + ["  " + f for f in files[:5]] + ["", "# rest"] + files[5:] + [""]
    (inputs / "wrf_files_messy.txt").write_bytes("\r\n".join(lines).encode())
    outdir, error = run_model(exe, inputs, "list_messy", ["-pl", "50000", "-fmt", "bin"],
                              list_name="wrf_files_messy.txt")
    assert error is None, error
    full = run_case("pl50000")
    for name in OUTPUTS:
        assert (output_path(outdir, name, basename="wrf_files_messy").read_bytes()
                == output_path(full, name).read_bytes()), name


def test_no_pp_et(exe, inputs, run_case):
    """--no_pp_et: no PP and ET, the other outputs unchanged; also with -rf 0, which
    PP cannot use (it needs the file after the last step)."""
    outdir, error = run_model(exe, inputs, "no_pp_et", ["-pl", "50000", "-fmt", "bin", NO_PP_ET])
    assert error is None, error
    assert not output_path(outdir, "PP").exists() and not output_path(outdir, "ET").exists()
    full = run_case("pl50000")
    for name in SLAB_OUTPUTS:
        assert output_path(outdir, name).read_bytes() == output_path(full, name).read_bytes(), name
    outdir, error = run_model(exe, inputs, "no_pp_et_nc", ["-pl", "50000", NO_PP_ET])
    assert error is None, error
    data, _, _ = read_nc(nc_path(outdir))
    assert "PP" not in data and "ET" not in data

    _, error = run_model(exe, inputs, "rf0_pp", ["-pl", "50000"], list_name=LIST_SHIFTED, rem_file=0)
    assert error is not None and "--remfile must be 1 or 2" in error


# netCDF4 itself triggers a numpy 2.5 DeprecationWarning when writing
@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_pp_buckets(exe, tmp_path):
    """Files without the bucket option (all the rain in RAINNC and RAINC, I_RAINNC = 0, no
    BUCKET_MM attribute): PP from the rain variables alone. Stored in float32, the rain
    is rounded (to ~6e-5 at 600 mm), so the expected PP is calculated from these files."""
    prepare_inputs(tmp_path)
    for f in tmp_path.glob("wrfout_*.nc"):
        with netCDF4.Dataset(f, "a") as nc:
            nc.set_auto_mask(False)
            for var in ["RAINNC", "RAINC"]:
                nc[var][0, :, :] = nc[var][0].astype("f8") + nc.BUCKET_MM * nc[f"I_{var}"][0]
                nc[f"I_{var}"][0, :, :] = np.zeros((NY, NX), "i4")
            nc.delncattr("BUCKET_MM")
    outdir, error = run_model(exe, tmp_path, "no_buckets", ["-pl", "50000", "-fmt", "bin"])
    assert error is None, error
    errors = compare_grid(read_output(outdir, "PP"), pp_from_wrf(case_dir=tmp_path))
    assert not errors, "\n".join(errors)


def test_nseconds(exe, inputs, run_case):
    """--nseconds is optional and must agree with the WRF Times (1 h between files)."""
    outdir, error = run_model(exe, inputs, "nsec_ok", ["-pl", "50000", "-fmt", "bin", "-nsec", "3600"])
    assert error is None, error
    for name in ["PWflux", "ET"]:
        assert output_path(outdir, name).read_bytes() == output_path(run_case("pl50000"), name).read_bytes()
    _, error = run_model(exe, inputs, "nsec_bad", ["-pl", "50000", "-nsec", "1800"])
    assert error is not None and "does not agree with the WRF Times" in error
