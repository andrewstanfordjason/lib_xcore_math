// Copyright 2020-2022 XMOS LIMITED.
// This Software is subject to the terms of the XMOS Public Licence: Version 1.

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>

#include "xmath/xmath.h"

#include "../tst_common.h"

#include "unity_fixture.h"

TEST_GROUP_RUNNER(float_sXX_add) {
  RUN_TEST_CASE(float_sXX_add, float_s32_add);
  RUN_TEST_CASE(float_sXX_add, float_complex_s16_add);
  RUN_TEST_CASE(float_sXX_add, float_complex_s32_add);
}

TEST_GROUP(float_sXX_add);
TEST_SETUP(float_sXX_add) { fflush(stdout); }
TEST_TEAR_DOWN(float_sXX_add) {}


#if SMOKE_TEST
#  define REPS       (1000)
#else
#  define REPS       (10000)
#endif



TEST(float_sXX_add, float_s32_add)
{
  unsigned seed = SEED_FROM_FUNC_NAME();

  float_s32_t x, y;
  conv_error_e error = 0;
  
  for(int v = 0; v < REPS; v++){

    setExtraInfo_RS(v, seed);

    x.exp = pseudo_rand_int(&seed, -40, 40);
    y.exp = x.exp + pseudo_rand_int(&seed, -5, 6);

    right_shift_t x_shr = pseudo_rand_uint(&seed, 0, 14);
    right_shift_t y_shr = pseudo_rand_uint(&seed, 0, 14);

    x.mant = pseudo_rand_int32(&seed) >> x_shr;
    y.mant = pseudo_rand_int32(&seed) >> y_shr;
    
    double expected_f = ldexp(x.mant, x.exp) + ldexp(y.mant, y.exp);

    float_s32_t actual = float_s32_add(x, y);
    
    int32_t expected = conv_double_to_s32(expected_f, actual.exp, &error);

    TEST_ASSERT_EQUAL_MESSAGE(0, error, "Conversion error.");

    TEST_ASSERT_INT32_WITHIN(2, expected, actual.mant);

  }
}


TEST(float_sXX_add, float_complex_s16_add)
{
  unsigned seed = SEED_FROM_FUNC_NAME();

  float_complex_s16_t x, y;
  conv_error_e error = 0;
  
  for(int v = 0; v < REPS; v++){

    setExtraInfo_RS(v, seed);

    x.exp = pseudo_rand_int(&seed, -20, 20);
    y.exp = x.exp + pseudo_rand_int(&seed, -6, 5);

    x.mant.re = pseudo_rand_int16(&seed);
    x.mant.im = pseudo_rand_int16(&seed);
    y.mant.re = pseudo_rand_int16(&seed);
    y.mant.im = pseudo_rand_int16(&seed);

    double expected_re = ldexp(x.mant.re, x.exp) + ldexp(y.mant.re, y.exp);
    double expected_im = ldexp(x.mant.im, x.exp) + ldexp(y.mant.im, y.exp);

    float_complex_s16_t actual = float_complex_s16_add(x, y);
    
    float_complex_s16_t expected = { 
      { conv_double_to_s16(expected_re, actual.exp, &error),
        conv_double_to_s16(expected_im, actual.exp, &error) },
      actual.exp
    };

    TEST_ASSERT_EQUAL_MESSAGE(0, error, "Conversion error.");

    TEST_ASSERT_INT16_WITHIN(2, expected.mant.re, actual.mant.re);
    TEST_ASSERT_INT16_WITHIN(2, expected.mant.im, actual.mant.im);

  }
}


TEST(float_sXX_add, float_complex_s32_add)
{
  unsigned seed = SEED_FROM_FUNC_NAME();

  float_complex_s32_t x, y;
  conv_error_e error = 0;
  
  for(int v = 0; v < REPS; v++){

    setExtraInfo_RS(v, seed);

    x.exp = pseudo_rand_int(&seed, -20, 20);
    y.exp = x.exp + pseudo_rand_int(&seed, -6, 5);

    x.mant.re = pseudo_rand_int32(&seed);
    x.mant.im = pseudo_rand_int32(&seed);
    y.mant.re = pseudo_rand_int32(&seed);
    y.mant.im = pseudo_rand_int32(&seed);

    double expected_re = ldexp(x.mant.re, x.exp) + ldexp(y.mant.re, y.exp);
    double expected_im = ldexp(x.mant.im, x.exp) + ldexp(y.mant.im, y.exp);

    float_complex_s32_t actual = float_complex_s32_add(x, y);
    
    float_complex_s32_t expected = { 
      { conv_double_to_s32(expected_re, actual.exp, &error),
        conv_double_to_s32(expected_im, actual.exp, &error) },
      actual.exp
    };

    TEST_ASSERT_EQUAL_MESSAGE(0, error, "Conversion error.");

    TEST_ASSERT_INT32_WITHIN(2, expected.mant.re, actual.mant.re);
    TEST_ASSERT_INT32_WITHIN(2, expected.mant.im, actual.mant.im);

  }
}

