"""Regression tests: run prepro2LDRM-slabvars on the cropped AMZ2 CPM 4 km case.

    pytest test/regression --exe PATH [--mpiexec PATH] [--update-reference]

The pl50000 reference was cut from the outputs of the original AMZ2 run;
the pl95000 reference was produced by an earlier version of this program.
A failure means results changed, not that they are wrong. After an intended
change in results, rerun with --update-reference and commit the new references.
"""
import numpy as np
import pytest

from regression import (NX, NY, OUTPUTS, PRESSURE_LEVEL, REFERENCE_CASES, TIME_FACTOR, compare_grid,
                        output_path, read_output, read_p_lowest_level, read_reference, reference_path,
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
    for name in OUTPUTS:
        data = read_output(outdir, name)
        assert (np.isnan(data[:, 0]) == expected).all(), \
            f"{name}: NaN in {np.isnan(data[:, 0]).sum()} cells, expected {expected.sum()}"
        # the upper slabs always have data
        assert not np.isnan(data[:, 1:]).any(), f"{name}: NaN above the first level"


@pytest.mark.parametrize("output", OUTPUTS)
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
