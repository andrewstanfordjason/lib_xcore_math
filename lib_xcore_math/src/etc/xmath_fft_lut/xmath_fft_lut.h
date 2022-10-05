// Copyright 2021 XMOS LIMITED. This Software is subject to the terms of the
// XMOS Public License: Version 1
#pragma once
#include "xmath/xmath.h"

extern const complex_s16_t xmath_dit_fft_lut_16[1024]; // 4096 bytes

/** @brief Maximum FFT length (log2) that can be performed using decimation-in-time. */
#define MAX_DIT_FFT_LOG2_16 10

/** @brief Convenience macro to index into the decimation-in-time FFT look-up table. 

	This will return the address at which the coefficients for the final pass of the real DIT
	FFT algorithm begin. 

	@param N	The FFT length.
*/
#define XMATH_DIT_REAL_FFT_LUT_16(N) &xmath_dit_fft_lut_16[(N)-8]

extern const complex_s16_t xmath_dif_fft_lut_16[1024]; // 4096 bytes

/** @brief Maximum FFT length (log2) that can be performed using decimation-in-frequency. */
#define MAX_DIF_FFT_LOG2_16 10

/** @brief Convenience macro to index into the decimation-in-frequency FFT look-up table. 

	Use this to retrieve the correct address for the DIF FFT look-up table when performing
	an FFT (or IFFT) using the DIF algorithm. (@see fft_dif_forward).

	@param N	The FFT length.
*/
#define XMATH_DIF_FFT_LUT_16(N) &xmath_dif_fft_lut_16[(1<<(MAX_DIF_FFT_LOG2)) - (N)]

extern const complex_s32_t xmath_dit_fft_lut_32[1020]; // 8160 bytes

/** @brief Maximum FFT length (log2) that can be performed using decimation-in-time. */
#define MAX_DIT_FFT_LOG2 10

/** @brief Convenience macro to index into the decimation-in-time FFT look-up table. 

	This will return the address at which the coefficients for the final pass of the real DIT
	FFT algorithm begin. 

	@param N	The FFT length.
*/
#define XMATH_DIT_REAL_FFT_LUT(N) &xmath_dit_fft_lut_32[(N)-8]

extern const complex_s32_t xmath_dif_fft_lut_32[1020]; // 8160 bytes

/** @brief Maximum FFT length (log2) that can be performed using decimation-in-frequency. */
#define MAX_DIF_FFT_LOG2 10

/** @brief Convenience macro to index into the decimation-in-frequency FFT look-up table. 

	Use this to retrieve the correct address for the DIF FFT look-up table when performing
	an FFT (or IFFT) using the DIF algorithm. (@see fft_dif_forward).

	@param N	The FFT length.
*/
#define XMATH_DIF_FFT_LUT(N) &xmath_dif_fft_lut_32[(1<<(MAX_DIF_FFT_LOG2)) - (N)]

