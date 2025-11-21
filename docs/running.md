# Running the Code

## Step 0: Running DFT and DFPT

Before using LetzElPhC, ensure you have:

- Kohn-Sham wavefunctions from an SCF calculation
- Phonon eigenvectors and perturbed potentials from DFPT

### Steps with Quantum Espresso

1. Run **SCF** calculation.  
2. Run **DFPT** calculation using `ph.x` to obtain dynamical matrices and potential changes.  
3. Run **NSCF** calculation on a uniform k-point grid.  

**Note:** The q-grid must be commensurate with the k-grid. Use `dvscf` flag to save changes in potentials.  

### Prepare SAVE folder

```bash
$ p2y
$ yambo
```

This generates the **SAVE** folder needed by LetzElPhC.

## Step 1: Preprocessor

Create the `ph_save` folder:

```bash
$ cd /path/to/phonon_calc_dir
$ lelphc -pp --code=qe -F PH.X_input_file
```

Optional environment variables:

```bash
$ export ELPH_PH_SAVE_DIR=ph_save_name_you_want
# Use symlinks instead of copying files
$ export ELPH_COPY_CMD="ln -sr"
```

## Final Step: Perform ELPH Calculation

Run the calculation in any folder with the LetzElPhC input file:

```bash
$ mpirun -n 4 lelphc -F LetzElPhC_input_file
```

### Example input variables:

```make
nkpool      = 1   # k point parallelization
nqpool      = 1   # q point parallelization
start_bnd   = 1   # starting band
end_bnd     = 40  # last band
save_dir    = SAVE
ph_save_dir = ph_save
kernel      = dfpt
convention  = standard
```
## El-ph Calculation for Yambo via LetzElPhC using Yambopy

This section shows how **Yambopy** can be used to run LetzElPhC and generate NetCDF databases for Yambo.

### Requirements

- LetzElPhC installed
- Yambopy installed
- `pw.x`, `ph.x`, `p2y`, `yambo` available

### SCF Calculation

Run a standard SCF calculation with symmorphic symmetries:

```text
force_symmorphic = .true.  # in system card
```

### NSCF Calculation

Copy the SCF `save` directory and run NSCF for desired empty states.

#### DVSCF Calculation

```text
prefix_dvscf
&inputph
  tr2_ph = 1.0d-14,
  verbosity = 'high',
  prefix = 'prefix',
  fildvscf = 'prefix-dvscf',
  electron_phonon = 'dvscf',
  fildyn = 'prefix.dyn',
  epsil = .false.,
  ldisp = .true.,
  recover = .true.,
  nq1 = Nx,
  nq2 = Ny,
  nq3 = Nz
/
```

### Yambo SAVE Directory

Run:

```bash
$ p2y
$ yambo
```

Move the **SAVE** directory to a convenient location.

### Obtain El-ph Databases

```bash
yambopy l2y
```

For serial run (bands `n_i` to `n_f`):

```bash
yambopy l2y -ph path/of/ph_input.in -b n_i n_f
```

For parallel calculation with 4 qpools and 2 kpools:

```bash
yambopy l2y -ph path/of/ph_input.in -b n_i n_f -par 4 2 -lelphc path/to/lelphc_exe -D
```

Check the `SAVE` folder:

```bash
ls SAVE/ndb.elph*
```

This contains the Yambo-compatible databases.

