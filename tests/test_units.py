import pytest
import subprocess
import os
import glob

def test_c_unit_tests(pytestconfig):
    if not pytestconfig.getoption("unittest"):
        pytest.skip("Skipping C unit tests. Pass --unittest to execute them.")

    base_dir = os.path.dirname(__file__)
    unit_test_dir = os.path.join(base_dir, "unit_tests")
    
    # 1. Compile the unit tests
    print("\n=> Compiling C unit tests...")
    compile_cmd = "make clean && make"
    res_compile = subprocess.run(compile_cmd, cwd=unit_test_dir, shell=True, capture_output=True, text=True)
    
    if res_compile.returncode != 0:
        pytest.fail(f"Failed to compile unit tests.\nStdout:\n{res_compile.stdout}\nStderr:\n{res_compile.stderr}")
        
    # 2. Find the generated executables
    # The Makefile generates binaries with no extension from the .c files
    c_files = glob.glob(os.path.join(unit_test_dir, "*.c"))
    binaries = [f[:-2] for f in c_files]
    
    # 3. Run each executable
    for binary in binaries:
        if not os.path.exists(binary):
            pytest.fail(f"Expected compiled binary {binary} was not found.")
            
        bin_name = os.path.basename(binary)
        print(f"\n=> Running unit test: {bin_name}")
        
        res_run = subprocess.run(f"./{bin_name}", cwd=unit_test_dir, shell=True, capture_output=True, text=True)
        
        # Log to the same TEST_RUN_OUTPUT for consistency
        with open(os.path.join(base_dir, "TEST_RUN_OUTPUT"), "a") as out_f:
            out_f.write(f"\n--- UNIT TEST: {bin_name} ---\n")
            out_f.write(res_run.stdout)
            out_f.write(res_run.stderr)
            
        if res_run.returncode != 0:
            pytest.fail(f"Unit test {bin_name} failed.\nStdout:\n{res_run.stdout}\nStderr:\n{res_run.stderr}")
