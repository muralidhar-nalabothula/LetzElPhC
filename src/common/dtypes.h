/**
 * @file data_structures.h
 * @brief Core data structures for electron-phonon coupling calculations
 * 
 * Defines all major data structures used throughout the ELPH code including:
 * lattice information, phonon data, pseudopotentials, wavefunctions, MPI
 * communication topology, and user input parameters.
 * 
 * @note Convention: Always use struct/enum keywords explicitly when declaring
 *       variables (e.g., "struct lattice a;" not "lattice a;")
 */

#pragma once
#include <mpi.h>
#include <stdbool.h>

#include "elphC.h"
#include "error.h"

/**
 * @enum ELPH_dft_code
 * @brief Supported DFT code backends
 */
enum ELPH_dft_code
{
    DFT_CODE_QE /**< Quantum ESPRESSO */
};

/**
 * @enum ELPH_screening
 * @brief Electronic screening types for electron-phonon coupling
 */
enum ELPH_screening
{
    ELPH_NO_SCREENING,   /**< Bare (unscreened) interaction */
    ELPH_DFPT_SCREENING  /**< DFPT-calculated screening */
};

/**
 * @struct kernel_info
 * @brief Information about the electron-phonon coupling kernel type
 */
struct kernel_info
{
    char name_str[32];              /**< Human-readable kernel name */
    bool bare_loc;                  /**< Include bare local potential */
    bool non_loc;                   /**< Include non-local pseudopotential */
    enum ELPH_screening screening;  /**< Type of electronic screening */
};

/**
 * @enum calc_type
 * @brief Type of calculation to perform
 */
enum calc_type
{
    CALC_ELPH,            /**< Electron-phonon coupling calculation */
    CALC_PH_SAVE_CREATE,  /**< Preprocessing: create phonon save directory */
    CALC_HELP,            /**< Display help information */
    CALC_VERSION,         /**< Print version information */
    CALC_INTERPOLATION    /**< Interpolation calculation */
};

/**
 * @struct calc_details
 * @brief Calculation configuration and input file information
 */
struct calc_details
{
    enum calc_type calc;        /**< Type of calculation to perform */
    enum ELPH_dft_code code;    /**< DFT code backend */
    char input_file[512];       /**< Path to input file (ELPH or DFT-Phonon) */
};

/**
 * @struct symmetry
 * @brief Crystallographic symmetry operation
 * 
 * Represents a symmetry operation as: \f$ \mathbf{r}' = R\mathbf{r} + \boldsymbol{\tau} \f$.
 * In case of time reversal symmetry, R -> -R and tau -> -tau and time_rev is true 
 */
struct symmetry
{
    ELPH_float Rmat[9];  /**< 3×3 rotation matrix (row-major, Cartesian coords) */
    ELPH_float tau[3];   /**< Fractional translation vector (Cartesian coords) */
    bool time_rev;       /**< True if operation is time reversal */
};

/**
 * @struct Lattice
 * @brief Crystal structure and electronic structure information
 * 
 * Contains all information about the crystal lattice, atomic positions,
 * symmetries, k-point sampling, and band structure parameters.
 */
struct Lattice
{
    int natom;                   /**< Number of atoms in unit cell */
    int nsym;                    /**< Number of symmetries in NSCF calculation */
    int timerev;                 /**< System has time-reversal symmetry */
    int nspin;                   /**< Number of spin components (1 or 2) */
    int nspinor;                 /**< Number of spinor components (1 or 2) */
    int total_bands;             /**< Total bands in NSCF calculation */
    int start_band;              /**< Starting band for ELPH calculation */
    int end_band;                /**< Ending band for ELPH calculation */
    int nbnds;                   /**< Number of bands used in ELPH (end-start+1) */
    char dimension;              /**< '1','2','3' for 1D/2D/3D (for Coulomb cutoff) */
    ND_int nmag;                 /**< Magnetic components: 1 (none or non-magnetic non-collinear), 2 (LSDA), 4 (magnetic non-collinear) */
    ND_int nmodes;               /**< Number of phonon modes (3×natom) */
    ND_int nkpts_iBZ;            /**< Number of k-points in irreducible BZ */
    ND_int nkpts_BZ;             /**< Number of k-points in full BZ */
    ND_int npw_max;              /**< Maximum plane waves in spherical grid for all wavefunctions */
    ND_int fft_dims[3];          /**< FFT grid dimensions [nx, ny, nz] */
    ND_int nfftz_loc;            /**< Number of FFT z-vectors on this MPI rank */
    ND_int nfftz_loc_shift;      /**< Global index of first FFT z-vector on this rank */
    ELPH_float alat_vec[9];      /**< Lattice vectors in Cartesian (row-major: a₁,a₂,a₃) */
    ELPH_float blat_vec[9];      /**< Reciprocal lattice vectors (includes 2π) */
    ELPH_float volume;           /**< Unit cell volume: \f$ V = \det(\mathbf{a}) \f$ */
    ELPH_float* atomic_pos;      /**< Atomic positions in Cartesian coords (natom×3) */
    int* atom_type;              /**< Atom type indices (natom) */
    ELPH_float* kpt_iredBZ;      /**< k-points in irreducible BZ (Cartesian, nkpts_iBZ×3) */
    ELPH_float* kpt_fullBZ;      /**< k-points in full BZ (Cartesian, nkpts_BZ×3) */
    ELPH_float* kpt_fullBZ_crys; /**< k-points in full BZ (crystal coords, nkpts_BZ×3) */
    int* kmap;                   /**< Map full BZ to iBZ: (nkpts_BZ×2) → [iBZ_idx, sym_idx] */
    struct symmetry* syms;       /**< Array of symmetry operations (nsym) */
    bool is_soc_present;         /**< True if spin-orbit coupling is present */
};

/**
 * @struct Phonon
 * @brief Phonon structure and properties
 * 
 * Stores phonon q-point sampling, symmetries, Born effective charges,
 * and dielectric properties for polar materials.
 */
struct Phonon
{
    ND_int nq_iBZ;               /**< Number of q-points in irreducible BZ */
    ND_int nq_BZ;                /**< Number of q-points in full BZ */
    ND_int nq_iBZ_loc;           /**< Number of q-points on this q-pool */
    ND_int nq_shift;             /**< Global q-index offset: [nq_shift, nq_shift+nq_iBZ_loc) */
    ND_int nph_sym;              /**< Number of phonon symmetries */
    ELPH_float* qpts_iBZ;        /**< q-points in iBZ (crystal units, nq_iBZ×3) */
    ELPH_float* qpts_BZ;         /**< q-points in full BZ (crystal units, nq_BZ×3) */
    struct symmetry* ph_syms;    /**< Phonon symmetry operations (nph_sym) */
    int* qmap;                   /**< Map full BZ to iBZ: (nq_BZ×2) → [iBZ_idx, sym_idx] */
    ND_int* nqstar;              /**< Number of q-points in each star (nq_iBZ) */
    ELPH_float* epsilon;         /**< Dielectric tensor ε∞ (3×3), NULL if not available */
    ELPH_float* Zborn;           /**< Born effective charges Z* (natom×3×3), NULL if not available */
    ELPH_float* Qpole;           /**< Quadrupole tensor (natom×3×3×3), currently unused (NULL) */
};

/**
 * @struct Vloc_table
 * @brief Interpolation table for local pseudopotential in reciprocal space
 * 
 * Provides tabulated values of local pseudopotential V(G) on a coarse grid
 * for efficient interpolation during calculations.
 */
struct Vloc_table
{
    ELPH_float qmax_abs;   /**< Maximum absolute q-coordinate (crystal units): max(|qₓ|,|qᵧ|,|qᵤ|) */
    ND_int npts_co;        /**< Number of points in interpolation table */
    ELPH_float* g_co;      /**< |G| grid points (coarse grid, npts_co) */
    ELPH_float dg;         /**< Spacing between grid points: Δ|G| */
    ELPH_float* vlocg;     /**< V(G) on coarse grid (ntype×npts_co) */
    ELPH_float* vploc_co;  /**< dV/d|G| on coarse grid (ntype×npts_co) */
};

/**
 * @struct local_pseudo
 * @brief Local part of pseudopotential in real space
 */
struct local_pseudo
{
    ELPH_float* Vloc_atomic; /**< Local pseudopotential (ntype) */
    ELPH_float* r_grid;      /**< Radial grid points r (ntype×ngrid) */
    ELPH_float* rab_grid;    /**< Radial grid derivative dr/dξ (ntype×ngrid) */
    ND_int ngrid;            /**< Number of radial grid points */
    ELPH_float Zval;         /**< Number of valence electrons */
};

/**
 * @struct Pseudo
 * @brief Complete pseudopotential information
 * 
 * Contains both local and non-local (Kleinman-Bylander) parts of the
 * pseudopotential for all atom types.
 */
struct Pseudo
{
    struct local_pseudo* loc_pseudo; /**< Local pseudopotential data (ntype) */
    ELPH_float* PP_table;            /**< PP table from Yambo format (nltimesj×ntype×3) */
    ELPH_float* Fsign;               /**< Sign of KB coefficients (nlj_max×ntype) */
    ELPH_cmplx** fCoeff;             /**< KB form factors per Eq.9 of PRB 71, 115106 (2005)
                                          (ntype)[array of (spin×spin×(2l+1)×(2l+1))] */
    struct Vloc_table vloc_table[1]; /**< Interpolation table for local potential */
    ND_int nltimesj;                 /**< Maximum number of projectors (n×j) */
    ND_int ngrid_max;                /**< Maximum radial grid size: max(len(r_grid)) */
    ND_int ntype;                    /**< Number of atom types */
    int lmax;                        /**< Maximum angular momentum l in PP table */
};

/**
 * @struct WFC
 * @brief Wavefunction data in reciprocal space
 * 
 * Stores wavefunctions and associated plane-wave coefficients for k-points
 * in the irreducible Brillouin zone.
 */
struct WFC
{
    ELPH_cmplx* wfc;    /**< Wavefunction coefficients (nspin×nbnd×nspinor×npw_loc) */
    ELPH_float* gvec;   /**< G-vectors (Cartesian, no 2π, npw_loc×3) */
    ELPH_float* Fk;     /**< Kleinman-Bylander factors in k-space (nltimesj×ntype×npw_loc)
                             Fₖ in Eq.6 of ABINIT pseudopotential theory docs */
    ND_int npw_total;   /**< Total number of G-vectors for this wavefunction */
    ND_int npw_loc;     /**< Number of G-vectors on this MPI rank */
};

/**
 * @struct ELPH_MPI_Comms
 * @brief MPI communication topology for parallel calculations
 * 
 * Defines a hierarchical parallelization scheme:
 * - Level 1: q-point parallelization (q-pools)
 * - Level 2: k-point parallelization (k-pools within each q-pool)
 * 
 * Communication structure:
 * @code
 *                    Global MPI_COMM_WORLD {0,1,2,3,4,5,6,7}
 *                              |
 *            __________________|__________________
 *            |                                    |
 *     {0,1,2,3} (commQ)                    {4,5,6,7} (commQ)
 *    _____|_____                          _____|_____
 *    |         |                          |         |
 * {0,1}     {2,3}                      {4,5}     {6,7}
 * (commK)   (commK)                    (commK)   (commK)
 * 
 * commR:  {0,2,4,6}, {1,3,5,7}  (same rank across all k-pools)
 * commRq: {0,2}, {1,3} (in 1st q-pool), {4,6}, {5,7} (in 2nd q-pool)
 * @endcode
 */
struct ELPH_MPI_Comms
{
    MPI_Comm commW;   /**< MPI_COMM_WORLD (do not free!) */
    MPI_Comm commQ;   /**< Communicator for q-pool */
    MPI_Comm commK;   /**< Communicator for k-pool */
    MPI_Comm commR;   /**< Equal-rank CPUs across all k-pools in COMM_WORLD */
    MPI_Comm commRq;  /**< Equal-rank CPUs across all k-pools within commQ */
    
    int nqpools;      /**< Total number of q-pools */
    int nkpools;      /**< Total number of k-pools (per q-pool) */
    
    /* Ranks in each communicator */
    int commW_rank;   /**< Rank in COMM_WORLD */
    int commQ_rank;   /**< Rank in q-pool */
    int commK_rank;   /**< Rank in k-pool */
    int commR_rank;   /**< Rank in equal-rank communicator */
    int commRq_rank;  /**< Rank in q-pool equal-rank communicator */
    
    /* Sizes of each communicator */
    int commW_size;   /**< Size of COMM_WORLD */
    int commQ_size;   /**< Size of q-pool */
    int commK_size;   /**< Size of k-pool */
    int commR_size;   /**< Size of equal-rank communicator */
    int commRq_size;  /**< Size of q-pool equal-rank communicator */
};

/**
 * @struct usr_input
 * @brief User input parameters from input file
 */
struct usr_input
{
    int nkpool;          /**< Number of k-point pools */
    int nqpool;          /**< Number of q-point pools */
    int start_bnd;       /**< Starting band index */
    int end_bnd;         /**< Ending band index */
    char* save_dir;      /**< Path to DFT save directory */
    char* ph_save_dir;   /**< Path to phonon save directory */
    char* kernel_str;    /**< Kernel screening level specification */
    bool kminusq;        /**< True for Yambo convention (k-q), false otherwise */
};
