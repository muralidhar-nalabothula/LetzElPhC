# LetzElPhC Test Suite Manual

Welcome to the LetzElPhC `pytest`-based test suite. This document serves as the comprehensive manual for understanding, running, and extending the testing framework.

## Table of Contents
1. [Overview & Architecture](#overview--architecture)
2. [Running the Test Suite](#running-the-test-suite)
3. [How the Test Suite Works (Under the Hood)](#how-the-test-suite-works-under-the-hood)
4. [Adding a New Test Case](#adding-a-new-test-case)
5. [Adding a New Run Mode](#adding-a-new-run-mode)

---

## Overview & Architecture

The LetzElPhC test suite has been completely modernized using the `pytest` framework. It replaces the old, monolithic `driver.py` and `run_tests.py` scripts with a dynamic, highly parallelized, and modular architecture. 

The suite automatically:
- Discovers and parses configurations from `.ini` files.
- Tests every valid parallelization strategy (MPI) across `k`-point and `q`-point pools dynamically.
- Automatically handles Single vs. Double precision data conversion.
- Cleans up temporary test files to ensure reproducible, isolated test environments.

---

## Running the Test Suite

Because we leverage `pytest`, executing your tests is highly customizable. The framework requires you to explicitly declare what you want to run.

### 1. Download Test Data
The large `.SAVE` databases and reference files are hosted externally. Before running tests for the first time, you must pull them:
```bash
pytest --download
```
*Note: This will securely clone the `Lelphc-Test-Data` repository directly into `tests/dft_qe` and immediately exit.*

### 2. Execution Commands
To run the tests, you must specify the suite type. You can combine these flags.

**Run the C Unit Tests:**
```bash
pytest --unittest
```

**Run the Integration Test Suite:**
```bash
pytest --testsuite --ncpus=4 --mpirun=mpirun --lelphc=../src/lelphc
```

**Run Both:**
```bash
pytest --unittest --testsuite --ncpus=4 --mpirun=mpirun --lelphc=../src/lelphc
```

### 3. Custom Arguments
- `--lelphc=PATH` : Absolute or relative path to the `lelphc` executable (Default: `./lelphc`).
- `--mpirun=CMD`  : MPI run command (Default: `mpirun`).
- `--ncpus=N`     : Max number of processors to test parallelization against (Default: `1`).

### Debugging Failures
If a test fails, `pytest` will print the exact traceback. For deeper inspection, the suite logs the complete `stdout` and `stderr` of every single run directly to a `TEST_RUN_OUTPUT` file located in that specific test's directory (e.g., `tests/dft_qe/Si_bulk/TEST_RUN_OUTPUT`).

---

## How the Test Suite Works (Under the Hood)

When you type `pytest --testsuite`, the following lifecycle occurs:

1. **Discovery (`test_list.py`)**: Pytest reads `test_list.py` to identify which directories contain active tests.
2. **Configuration Parsing (`test.ini`)**: It steps into each directory and parses the `test.ini` file, merging the `[DEFAULT]` block with the specific test block (e.g., `[test1 qe]`).
3. **Data Preparation**: It auto-detects the compiled precision of your executable (`lelphc --version`). It then dynamically builds the binary `ph_save` files and converts the `SAVE` databases to match your compiled precision exactly.
4. **Parallelization Math**: It determines the maximum allowed `k`-pools and `q`-pools from the reference database. It then calculates every mathematically valid factorization (using `get_triplet()`) for the requested `--ncpus`.
5. **Execution**: For every single pool combination, it writes a custom `test_input.in` file and executes the code via MPI.
6. **Validation**: It passes the generated files (e.g., `ndb.Dmats`) into `check_data.py` to strictly verify numerical accuracy against the reference files.

---

## Adding a New Test Case

Adding a standard test requires zero python scripting. You simply add your data and register it.

### Step 1: Create the Directory
Create a folder inside `tests/dft_qe/` (e.g., `tests/dft_qe/My_New_Material`). Place your `SAVE` directory, `ph_save` directory, and reference files (e.g., `ndb.elph_ref`) inside.

### Step 2: Create the `test.ini` Configuration
In your new folder, create `test.ini`.

```ini
[DEFAULT]
run_mode = standard
elph_db_ref = ndb.elph_ref
dmat_db_ref = ndb.Dmats_ref
# Any parameter defined here is injected into test_input.in automatically

[test1 qe]
# Specific parameters for test1
parameter_A = 10

[test2 qe]
# Specific parameters for test2
parameter_A = 20
```

### Step 3: Register the Test
Open `tests/test_list.py` and append your folder path to the array:
```python
tests = (
[
    ['dft_qe/Si_bulk'     , 'test.ini'],
    ['dft_qe/My_New_Test' , 'test.ini']  # <--- Your new test
])
```
Done! Pytest will automatically discover it.

---

## Adding a New Run Mode

The `run_mode` variable in `test.ini` tells the test suite how to execute the binary and how to validate the output. If you want to add a new physical workflow (e.g., `interpolation`, `wannier`), follow these steps.

### Step 1: Update the `.ini` File
In your new test directory's `test.ini`, define your new `run_mode`. Also, define any custom reference files you need.

*Crucial Note: Any parameter ending in `_ref` (e.g., `my_custom_db_ref`) is automatically safely ignored when writing the `test_input.in` file, but remains available to the python validation script!*

```ini
[DEFAULT]
run_mode = interpolation
interp_db_ref = ndb.interp_custom_ref
```

### Step 2: Add Validation Logic
Open `tests/check_data.py`. Write a new Python function that reads the output database and compares it to the reference database.

```python
def check_interp_files(file_out, file_ref):
    # Your logic using netCDF4 to compare numerical arrays
    return True # if arrays match perfectly
```

### Step 3: Hook it into `test_integration.py`
Open `tests/test_integration.py`. Locate the execution block (around line 118). Add any specific command-line flags required by your new mode:

```python
# Setup flags based on run mode
flags = f"-F test_input.in"
if run_mode == "interpolation":
    flags = f"-i -F test_input.in" # Example: injecting the -i flag
```

Next, scroll down to the validation block (around line 140) and hook in your new checking function from Step 2:

```python
# Verification based on run_mode
if run_mode == "standard":
    assert check_dmat_files('ndb.Dmats', dmat_db_ref), "Dmats mismatch"
    assert check_elph_files('ndb.elph', elph_db_ref), "elph DB mismatch"

elif run_mode == "interpolation":
    # Grab the reference filename from the config
    interp_ref = cfg.get('interp_db_ref', 'ndb.interp_ref')
    # Assert against your custom logic
    assert check_interp_files('ndb.interp', interp_ref), "Interpolation mismatch!"

else:
    pytest.fail(f"Unknown run_mode: {run_mode}")
```
You have successfully implemented a brand new physics run mode!
