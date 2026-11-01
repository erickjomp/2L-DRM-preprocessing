import shutil
from pathlib import Path

import pytest

import regression


def pytest_addoption(parser):
    parser.addoption("--exe", type=Path, help="path to the prepro2LDRM-metdata executable to test")
    parser.addoption("--mpiexec", type=Path, default=shutil.which("mpiexec"),
                     help="MPI launcher for the multi-process runs (default: mpiexec in PATH); "
                          "those tests are skipped without it")
    parser.addoption("--serial", action="store_true",
                     help="the executable was built without MPI (no USE_MPI): skip the multi-process runs")
    parser.addoption("--update-reference", action="store_true",
                     help="overwrite the reference outputs instead of comparing against them")


@pytest.fixture(scope="session")
def exe(request):
    path = request.config.getoption("--exe")
    if path is None:
        pytest.fail("pass the executable with --exe, e.g. "
                    "--exe $(fpm run prepro2LDRM-metdata --runner echo | tail -1)", pytrace=False)
    if not path.is_file():
        pytest.fail(f"--exe {path} does not exist", pytrace=False)
    return path.resolve()


@pytest.fixture(scope="session")
def mpiexec(request):
    if request.config.getoption("--serial"):
        return None
    return request.config.getoption("--mpiexec")


@pytest.fixture(scope="session")
def update_reference(request):
    return request.config.getoption("--update-reference")


@pytest.fixture(scope="session")
def inputs(tmp_path_factory):
    """Directory with the copied case. The program is run from here."""
    workdir = tmp_path_factory.mktemp("amz2_cpm_crop")
    regression.prepare_inputs(workdir)
    return workdir


@pytest.fixture(scope="session")
def run_case(exe, inputs, mpiexec):
    """run_case(case, nprocs) -> output directory. Each run happens only once."""
    done = {}

    def run(case, nprocs=1):
        if nprocs > 1 and mpiexec is None:
            pytest.skip("serial build (--serial) or no MPI launcher (mpiexec not in PATH, see --mpiexec)")
        if (case, nprocs) not in done:
            done[case, nprocs] = regression.run_model(
                exe, inputs, f"{case}_np{nprocs}", regression.CASES[case], nprocs, mpiexec)
        outdir, error = done[case, nprocs]
        if error:
            pytest.fail(f"[{case}, {nprocs} processes] {error}", pytrace=False)
        return outdir

    return run
