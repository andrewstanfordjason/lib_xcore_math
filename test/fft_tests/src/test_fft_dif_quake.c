// Copyright 2020-2022 XMOS LIMITED.
// This Software is subject to the terms of the XMOS Public Licence: Version 1.

#include "xmath/xmath.h"
#include "testing.h"
#include "floating_fft.h"
#include "tst_common.h"
#include "fft.h"
#include "unity_fixture.h"
#include "xmath_fft_lut.h"


TEST_GROUP_RUNNER(fft_dif_quake) {
  RUN_TEST_CASE(fft_dif_quake, fft_dif_forward_quake);
//   RUN_TEST_CASE(fft_dif_quake, fft_dif_inverse_quake);
//   RUN_TEST_CASE(fft_dif_quake, fft_dif_forward_complete_quake);
//   RUN_TEST_CASE(fft_dif_quake, fft_dif_inverse_complete_quake);
}

TEST_GROUP(fft_dif_quake);
TEST_SETUP(fft_dif_quake) { fflush(stdout); }
TEST_TEAR_DOWN(fft_dif_quake) {}


#define MAX_PROC_FRAME_LENGTH_LOG2 (MAX_DIF_FFT_LOG2)
#define MAX_PROC_FRAME_LENGTH (1<<MAX_PROC_FRAME_LENGTH_LOG2)


#define EXPONENT_SIZE 5
#define BASIC_HEADROOM 2
#define EXTRA_HEADROOM_MAX 3
#define WIGGLE 20

#define MIN_FFT_N_LOG2  (3)

#define LOOPS_LOG2 8


TEST(fft_dif_quake, fft_dif_forward_complete_quake)
{
#define FUNC_NAME "fft_dif_forward_complete_quake"

#if PRINT_FUNC_NAMES
    printf("\n%s..\n", FUNC_NAME);
#endif

    unsigned r = 0x6999B20C;
    conv_error_e error = 0;

    for(unsigned k = MIN_FFT_N_LOG2; k <= MAX_PROC_FRAME_LENGTH_LOG2; k++){
        const unsigned FFT_N = (1<<k);
        unsigned worst_case = 0;

        double sine_table[(MAX_PROC_FRAME_LENGTH/2) + 1];

        flt_make_sine_table_double(sine_table, FFT_N);
        
        for(unsigned t = 0; t < (1 << LOOPS_LOG2); t++){
            const unsigned seed = r;

            complex_s16_t DWORD_ALIGNED a[MAX_PROC_FRAME_LENGTH];
            complex_double_t DWORD_ALIGNED A[MAX_PROC_FRAME_LENGTH];
            double real[MAX_PROC_FRAME_LENGTH], imag[MAX_PROC_FRAME_LENGTH];

            exponent_t exponent = sext(pseudo_rand_int32(&r), EXPONENT_SIZE);
            right_shift_t extra_hr = BASIC_HEADROOM + (pseudo_rand_uint32(&r) % (EXTRA_HEADROOM_MAX+1));
            
            rand_vect_complex_s16(a, FFT_N, extra_hr, &r);
            conv_vect_complex_s16_to_complex_double(A, a, FFT_N, exponent, &error);
            conv_vect_complex_s16_to_complex_double_v2(real, imag, a, FFT_N, exponent, &error);
            TEST_ASSERT_CONVERSION(error);

            headroom_t headroom = vect_complex_s16_headroom_quake(a, FFT_N);

            Fft_transform(real, imag, FFT_N);

            flt_bit_reverse_indexes_double(A, FFT_N);
            flt_fft_forward_double(A, FFT_N, sine_table);
            fft_dif_forward_quake_s16(a, FFT_N, &headroom, &exponent);
            fft_index_bit_reversal_quake(a, FFT_N);

            unsigned diff = 0;
            print_vect_complex_double(A, FFT_N, &error);
            diff = abs_diff_vect_complex_s16(a, exponent, A, FFT_N, &error);
            TEST_ASSERT_CONVERSION(error);

            if(diff > worst_case) { worst_case = diff;  }
            TEST_ASSERT_LESS_OR_EQUAL_UINT32_MESSAGE(k+WIGGLE, diff, "Output delta is too large");

            if(diff > worst_case) { worst_case = diff;  }
            diff = abs_diff_vect_complex_s16(a, exponent, A, FFT_N, &error);
            TEST_ASSERT_CONVERSION(error);

            TEST_ASSERT_LESS_OR_EQUAL_UINT32_MESSAGE(k+WIGGLE, diff, "Output delta is too large");

            TEST_ASSERT_EQUAL_MESSAGE(vect_complex_s16_headroom_quake(a, FFT_N), headroom, "Reported headroom was incorrect.");
        }

#if PRINT_ERRORS
        printf("    %s worst error (%u-point): %u\n", FUNC_NAME, FFT_N, worst_case);
#endif

#if WRITE_PERFORMANCE_INFO
        fprintf(perf_file, "%s, %u, %u, -,\n", "fft_dif_forward", FFT_N, worst_case);
#endif

#undef FUNC_NAME
    }
}


TEST(fft_dif_quake, fft_dif_inverse_complete_quake)
{
#define FUNC_NAME "fft_dif_inverse_complete_quake"

#if PRINT_FUNC_NAMES
    printf("\n%s..\n", FUNC_NAME);
#endif

    unsigned r = 1;
    conv_error_e error = 0;

    for(unsigned k = MIN_FFT_N_LOG2; k <= MAX_PROC_FRAME_LENGTH_LOG2; k++){
        const unsigned FFT_N = (1<<k);
        unsigned worst_case = 0;

        double sine_table[(MAX_PROC_FRAME_LENGTH/2) + 1];

        flt_make_sine_table_double(sine_table, FFT_N);

        for(unsigned t = 0; t < (1<<LOOPS_LOG2); t++){
            const unsigned seed = r;

            complex_s16_t DWORD_ALIGNED a[MAX_PROC_FRAME_LENGTH];
            complex_double_t DWORD_ALIGNED A[MAX_PROC_FRAME_LENGTH];
            double real[MAX_PROC_FRAME_LENGTH], imag[MAX_PROC_FRAME_LENGTH];

            exponent_t exponent = sext(pseudo_rand_int32(&r), EXPONENT_SIZE);
            right_shift_t extra_hr = BASIC_HEADROOM + (pseudo_rand_uint32(&r) % (EXTRA_HEADROOM_MAX+1));
            
            rand_vect_complex_s16(a, FFT_N, extra_hr, &r);
            conv_vect_complex_s16_to_complex_double(A, a, FFT_N, exponent, &error);
            conv_vect_complex_s16_to_complex_double_v2(real, imag, a, FFT_N, exponent, &error);
            TEST_ASSERT_CONVERSION(error);
            
            headroom_t headroom = vect_complex_s16_headroom_quake(a, FFT_N);

            flt_bit_reverse_indexes_double(A, FFT_N);
            flt_fft_inverse_double(A, FFT_N, sine_table);

            fft_dif_inverse_quake_s16(a, FFT_N, &headroom, &exponent);
            fft_index_bit_reversal_quake(a, FFT_N);

            Fft_inverseTransform(real, imag, FFT_N);

            unsigned diff = 0;

            diff = abs_diff_vect_complex_s16(a, exponent, A, FFT_N, &error);
            TEST_ASSERT_CONVERSION(error);
            
            if(diff > worst_case) { worst_case = diff;  }
            TEST_ASSERT_LESS_OR_EQUAL_UINT32_MESSAGE(k+WIGGLE, diff, "Output delta is too large");
            TEST_ASSERT_CONVERSION(error);

            for(unsigned i = 0; i < FFT_N; i++){
                A[i].re = real[i];
                A[i].im = imag[i];
            }

            diff = abs_diff_vect_complex_s16(a, exponent, A, FFT_N, &error);
            TEST_ASSERT_CONVERSION(error);

            if(diff > worst_case) { worst_case = diff;  }
            TEST_ASSERT_LESS_OR_EQUAL_UINT32_MESSAGE(k+WIGGLE, diff, "Output delta is too large");

            TEST_ASSERT_EQUAL_MESSAGE(vect_complex_s16_headroom_quake(a, FFT_N), headroom, "Reported headroom was incorrect.");
        }


#if PRINT_ERRORS
        printf("    %s worst error (%u-point): %u\n", FUNC_NAME, FFT_N, worst_case);
#endif

#if WRITE_PERFORMANCE_INFO
        fprintf(perf_file, "%s, %u, %u, -,\n", "fft_dif_inverse", FFT_N, worst_case);
#endif

#undef FUNC_NAME
    }
}


TEST(fft_dif_quake, fft_dif_forward_quake)
{
#define FUNC_NAME "fft_dif_forward_quake"

#if PRINT_FUNC_NAMES
    printf("\n%s..\n", FUNC_NAME);
#endif

    unsigned r = 1;
    conv_error_e error = 0;

    for(unsigned k = MIN_FFT_N_LOG2; k <= MAX_PROC_FRAME_LENGTH_LOG2; k++){

        unsigned FFT_N = (1<<k);
        float worst_timing = 0.0f;
        double sine_table[(MAX_PROC_FRAME_LENGTH/2) + 1];

        flt_make_sine_table_double(sine_table, FFT_N);
        
        for(unsigned t = 0; t < (1<<LOOPS_LOG2); t++){

            complex_s16_t DWORD_ALIGNED a[MAX_PROC_FRAME_LENGTH];
            complex_double_t DWORD_ALIGNED A[MAX_PROC_FRAME_LENGTH];

            exponent_t exponent = sext(pseudo_rand_int32(&r), EXPONENT_SIZE);
            right_shift_t extra_hr = BASIC_HEADROOM + (pseudo_rand_uint32(&r) % (EXTRA_HEADROOM_MAX+1));
            
            rand_vect_complex_s16(a, FFT_N, extra_hr, &r);

            // for (int i=0;i<FFT_N;i++){
            //     a[i].re = INT16_MAX/4;
            //     a[i].im = 0;
            // }

            conv_vect_complex_s16_to_complex_double(A, a, FFT_N, exponent, &error);
            TEST_ASSERT_CONVERSION(error);
            
            headroom_t headroom = vect_complex_s16_headroom_quake(a, FFT_N);

            flt_bit_reverse_indexes_double(A, FFT_N);
            flt_fft_forward_double(A, FFT_N, sine_table);

            unsigned ts1 = getTimestamp();
            fft_dif_forward_quake_s16(a, FFT_N, &headroom, &exponent);
            unsigned ts2 = getTimestamp();
            fft_index_bit_reversal_quake(a, FFT_N);

            float timing = (ts2-ts1)/100.0;
            if(timing > worst_timing) worst_timing = timing;

            // print_vect_complex_double(A, FFT_N, &error);
            // print_vect_complex_s16(a, exponent, FFT_N, &error);

            // for (int i=0;i<8;i++)
            //     printf("(%d+%dj) ", a[i].re , a[i].im);
            // printf("\n");

            unsigned diff = abs_diff_vect_complex_s16(a, exponent, A, FFT_N, &error);
            TEST_ASSERT_CONVERSION(error);
            TEST_ASSERT_LESS_OR_EQUAL_UINT32_MESSAGE(k+WIGGLE, diff, "Output delta is too large");

            TEST_ASSERT_EQUAL_MESSAGE(vect_complex_s16_headroom_quake(a, FFT_N), headroom, "Reported headroom was incorrect.");
        }

#if TIME_FUNCS
        printf("    %s (%u-point): %f us\n", FUNC_NAME, FFT_N, worst_timing);
#endif

#if WRITE_PERFORMANCE_INFO
        fprintf(perf_file, "%s, %u,, %0.02f,\n", &(__func__[5]), FFT_N, worst_timing);
#endif

#undef FUNC_NAME
    }
}


TEST(fft_dif_quake, fft_dif_inverse_quake)
{
#define FUNC_NAME "fft_dif_inverse_quake"

#if PRINT_FUNC_NAMES
    printf("\n%s..\n", FUNC_NAME);
#endif

    unsigned r = 1;
    conv_error_e error = 0;

    for(unsigned k = MIN_FFT_N_LOG2; k <= MAX_PROC_FRAME_LENGTH_LOG2; k++){

        unsigned FFT_N = (1<<k);
        float worst_timing = 0.0f;
        double sine_table[(MAX_PROC_FRAME_LENGTH/2) + 1];

        flt_make_sine_table_double(sine_table, FFT_N);
        
        for(unsigned t = 0; t < (1<<LOOPS_LOG2); t++){
            
            complex_s16_t DWORD_ALIGNED a[MAX_PROC_FRAME_LENGTH];
            complex_double_t DWORD_ALIGNED A[MAX_PROC_FRAME_LENGTH];

            exponent_t exponent = sext(pseudo_rand_int32(&r), EXPONENT_SIZE);
            right_shift_t extra_hr = BASIC_HEADROOM + (pseudo_rand_uint32(&r) % (EXTRA_HEADROOM_MAX+1));
            
            rand_vect_complex_s16(a, FFT_N, extra_hr, &r);
            conv_vect_complex_s16_to_complex_double(A, a, FFT_N, exponent, &error);
            TEST_ASSERT_CONVERSION(error);
            
            headroom_t headroom = vect_complex_s16_headroom_quake(a, FFT_N);
            
            flt_bit_reverse_indexes_double(A, FFT_N);
            flt_fft_inverse_double(A, FFT_N, sine_table);

            unsigned ts1 = getTimestamp();
            fft_dif_inverse_quake_s16(a, FFT_N, &headroom, &exponent);
            unsigned ts2 = getTimestamp();
            fft_index_bit_reversal_quake(a, FFT_N);

            float timing = (ts2-ts1)/100.0;
            if(timing > worst_timing) worst_timing = timing;

            unsigned diff = abs_diff_vect_complex_s16(a, exponent, A, FFT_N, &error);
            TEST_ASSERT_CONVERSION(error);
            TEST_ASSERT_LESS_OR_EQUAL_UINT32_MESSAGE(k+WIGGLE, diff, "Output delta is too large");

            TEST_ASSERT_EQUAL_MESSAGE(vect_complex_s16_headroom_quake(a, FFT_N), headroom, "Reported headroom was incorrect.");
        }
#if TIME_FUNCS
        printf("    %s (%u-point): %f us\n", FUNC_NAME, FFT_N, worst_timing);
#endif

#if WRITE_PERFORMANCE_INFO
        fprintf(perf_file, "%s, %u,, %0.02f,\n", &(__func__[5]), FFT_N, worst_timing);
#endif
    }

#undef FUNC_NAME
}
