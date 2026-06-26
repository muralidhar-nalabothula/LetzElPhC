import pytest
import os

def pytest_addoption(parser):
    parser.addoption("--lelphc", action="store", default="./lelphc", help="Path to lelphc executable. Default is ./lelphc")
    parser.addoption("--mpirun", action="store", default="mpirun", help="Mpi run command. Default is mpirun")
    parser.addoption("--ncpus", action="store", default=1, type=int, help="Number of CPUs to use. Default is 1")
    parser.addoption("--download", action="store_true", help="Download test data and exit")
    parser.addoption("--unittest", "--unitest", action="store_true", dest="unittest", help="Run C unit tests")
    parser.addoption("--testsuite", action="store_true", help="Run integration test suite")

def pytest_sessionstart(session):
    """
    Called after the Session object has been created and
    before performing collection and entering the run test loop.
    Here we handle the download logic.
    """
    config = session.config
    if config.getoption("--download"):
        print("\n=> Downloading test data...")
        cwd = os.getcwd()
        dft_qe_dir = os.path.join(cwd, 'dft_qe')
        if not os.path.exists(dft_qe_dir):
            os.makedirs(dft_qe_dir)
        os.chdir(dft_qe_dir)
        if not os.path.exists(".git"):
            os.system("git init")
            os.system("git remote add origin https://github.com/muralidhar-nalabothula/Lelphc-Test-Data")
        os.system("git pull")
        os.system("git checkout main -f")
        os.system("git branch --set-upstream-to origin/main")
        os.chdir(cwd)
        pytest.exit("Download complete.")

def pytest_collection_modifyitems(config, items):
    """Skip items based on the provided flags."""
    run_unittest = config.getoption("unittest")
    run_testsuite = config.getoption("--testsuite")
    
    # If neither is provided, we skip all tests to avoid running them accidentally
    if not run_unittest and not run_testsuite and not config.getoption("--download"):
        skip_msg = pytest.mark.skip(reason="Pass --unittest or --testsuite to run tests.")
        for item in items:
            item.add_marker(skip_msg)
        return

    for item in items:
        if "test_units.py" in str(item.fspath):
            if not run_unittest:
                item.add_marker(pytest.mark.skip(reason="Not running unit tests. Use --unittest"))
        elif "test_integration.py" in str(item.fspath):
            if not run_testsuite:
                item.add_marker(pytest.mark.skip(reason="Not running test suite. Use --testsuite"))

def pytest_report_header(config):
    """Add custom header to pytest output."""
    import subprocess
    ncpus = config.getoption("--ncpus")
    mpirun = config.getoption("--mpirun")
    lelphc = config.getoption("--lelphc")
    
    lelphc_cmd = lelphc
    if lelphc_cmd.startswith(".") or "/" in lelphc_cmd:
        lelphc_cmd = os.path.abspath(lelphc_cmd)
        
    precision = "Unknown"
    try:
        ver_cmd = f"{mpirun} -n 1 {lelphc_cmd} --version"
        ver_res = subprocess.run(ver_cmd, shell=True, capture_output=True, text=True)
        if "Double" in ver_res.stdout:
            precision = "Double"
        elif "LetzElPhC" in ver_res.stdout or ver_res.stdout.strip() != "":
            precision = "Single"
    except Exception:
        pass
    
    header = []
    header.append("*" * 50)
    header.append("*" * 14 + " LetzElPhC Test Suite " + "*" * 14)
    header.append("=" * 50)
    header.append(f"NCPUS     : {ncpus}")
    header.append(f"MPIRUN    : {mpirun}")
    header.append(f"LELPHC    : {lelphc}")
    header.append(f"PRECISION : {precision}")
    header.append("=" * 50)
    return header
