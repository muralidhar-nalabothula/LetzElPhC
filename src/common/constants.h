/**
 * @file constants.h
 * @brief Physical and mathematical constants with precision selection
 *
 * Defines fundamental constants used throughout the ELPH library.
 * Constants are defined with appropriate precision (float/double) based on
 * the COMPILE_ELPH_DOUBLE compilation flag.
 */

#pragma once
#include "elphC.h"

#if defined(COMPILE_ELPH_DOUBLE)
/**
 * @def ELPH_EPS
 * @brief Numerical tolerance/epsilon for comparisons (double precision)
 *
 * Used for floating-point comparisons and convergence criteria.
 * Value: \f$ 10^{-6} \f$
 */
#define ELPH_EPS 1e-6

/**
 * @def ELPH_PI
 * @brief Mathematical constant pi (double precision)
 *
 * Value: \f$ \pi \approx 3.1415927 \f$
 */
#define ELPH_PI 3.1415927

/**
 * @def ELPH_SQRT2
 * @brief Square root of 2 (double precision)
 *
 * Value: \f$ \sqrt{2} \approx 1.4142136 \f$
 */
#define ELPH_SQRT2 1.4142136

/**
 * @def ELPH_e2
 * @brief Elementary charge squared in Rydberg atomic units (double precision)
 *
 * In Rydberg units: \f$ e^2 = 2 \f$ (Ry)
 *
 * @note This is the electronic charge squared used in electron-phonon
 *       coupling calculations in atomic Rydberg units
 */
#define ELPH_e2 2.0

#else
/**
 * @def ELPH_EPS
 * @brief Numerical tolerance/epsilon for comparisons (single precision)
 *
 * Used for floating-point comparisons and convergence criteria.
 * Value: \f$ 10^{-6} \f$
 */
#define ELPH_EPS 1e-6

/**
 * @def ELPH_PI
 * @brief Mathematical constant pi (single precision)
 *
 * Value: \f$ \pi \approx 3.1415927 \f$
 */
#define ELPH_PI 3.1415927f

/**
 * @def ELPH_SQRT2
 * @brief Square root of 2 (single precision)
 *
 * Value: \f$ \sqrt{2} \approx 1.4142136 \f$
 */
#define ELPH_SQRT2 1.4142136f

/**
 * @def ELPH_e2
 * @brief Elementary charge squared in Rydberg atomic units (single precision)
 *
 * In Rydberg units: \f$ e^2 = 2 \f$ (Ry)
 *
 * @note This is the electronic charge squared used in electron-phonon
 *       coupling calculations in atomic Rydberg units
 */
#define ELPH_e2 2.0f

#endif
