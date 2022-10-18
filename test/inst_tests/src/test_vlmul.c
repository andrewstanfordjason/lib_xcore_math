// Copyright 2020-2022 XMOS LIMITED.
// This Software is subject to the terms of the XMOS Public Licence: Version 1.

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>

#include "xmath/xmath.h"

#include "tst_common.h"

#include "unity_fixture.h"

TEST_GROUP_RUNNER(inst_vlmul) {
  RUN_TEST_CASE(inst_vlmul, vlmul_s8x8_basic);
  RUN_TEST_CASE(inst_vlmul, vlmul_s16x16_basic);
//   RUN_TEST_CASE(inst_vlmul, vlmul_s16x16_sweep);
  RUN_TEST_CASE(inst_vlmul, vlmul_s32x32_basic);
}

TEST_GROUP(inst_vlmul);
TEST_SETUP(inst_vlmul) { fflush(stdout); }
TEST_TEAR_DOWN(inst_vlmul) {}

static char msg_buff[200];

#define TEST_ASSERT_EQUAL_MSG(EXPECTED, ACTUAL, LINE_NUM)   do{       \
    if((EXPECTED)!=(ACTUAL)) {                                        \
      sprintf(msg_buff, "(test vector @ line %u)", (LINE_NUM));       \
      TEST_ASSERT_EQUAL_MESSAGE((EXPECTED), (ACTUAL), msg_buff);      \
    }} while(0)


#if SMOKE_TEST
#  define REPS       (100)
#else
#  define REPS       (100000)
#endif


void ref_vlmul_s8(int8_t * vR, int8_t * mem){
    for(int i=0;i<32;i++){
        int16_t v = (int16_t)(vR[i]) * (int16_t)(mem[i]);
        vR[i] = ( v + (1<<6)) >> 7;
    }
}

void ref_vlmul_s16(int16_t * vR, int16_t * mem){
    for(int i=0;i<16;i++){
        int32_t v = (int32_t)(vR[i]) * (int32_t)(mem[i]);
        vR[i] = ( v + (1<<14)) >> 15;
    }
}

// void ref_vlmul_s16xu16(int16_t * vR, uint16_t * mem){
//     for(int i=0;i<16;i++){
//         int32_t v = (int32_t)(vR[i]) * (uint32_t)mem[i];
//         vR[i] = ( v + (1<<15)) >> 16;
//     }
// }

void ref_vlmul_s32(int32_t * vR, int32_t * mem){
    for(int i=0;i<8;i++){
        int64_t v = (int64_t)(vR[i]) * (int64_t)(mem[i]);
        vR[i] = ( v + (1<<30)) >> 31;
    }
}

void VLMUL_0(int8_t * vR_out, int8_t * vD_out, const int8_t * vR_in, const int8_t * Mem_in, int mode);
void VLMUL_1(int8_t * vR_out, int8_t * vD_out, const int8_t * vD_in, const int8_t * vR_in, const int8_t * Mem_in, int mode);
void VLMUL_2(int8_t * vR_out, int8_t * vD_out, const int8_t * vD_in, const int8_t * vR_in, const int8_t * Mem_in, int mode);
void VLMUL_3(int8_t * vR_out, int8_t * vD_out, const int8_t * vD_in, const int8_t * vR_in, const int8_t * Mem_in, int mode);

#define S32_LEN 8
#define S16_LEN 16
#define S8_LEN 32
TEST(inst_vlmul, vlmul_s8x8_basic)
{
    
    unsigned seed = SEED_FROM_FUNC_NAME();

    int8_t A_ref[S8_LEN];
    int8_t A[S8_LEN];
    int8_t B[S8_LEN];

    for(int v = 0; v < REPS; v++){

        setExtraInfo_R(v);

        const unsigned len = S8_LEN;

        for(int i = 0; i < len; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            A[i] = pseudo_rand_int8(&seed) >> shr;
        }

        memcpy(A_ref, A, sizeof(A));

        for(int i = 0; i < len; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            B[i] = pseudo_rand_int8(&seed) >> shr;
        }

        VLMUL_0(A, 0, A, B, 8);
        ref_vlmul_s8(A_ref, B);

        for(int i = 0; i < len; i++)
            TEST_ASSERT_EQUAL_INT8(A[i], A_ref[i]);

    }
}


TEST(inst_vlmul, vlmul_s16x16_basic)
{
    
    unsigned seed = SEED_FROM_FUNC_NAME();

    int16_t A_ref[S16_LEN];
    int16_t A[S16_LEN];
    int16_t B[S16_LEN];
    int16_t vD[S16_LEN];

    for(int v = 0; v < REPS; v++){

        setExtraInfo_R(v);

        const unsigned len = S16_LEN;

        for(int i = 0; i < len; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            A[i] = pseudo_rand_int16(&seed) >> shr;
        }

        memcpy(A_ref, A, sizeof(A));

        for(int i = 0; i < len; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            B[i] = pseudo_rand_int16(&seed) >> shr;
        }

        VLMUL_0(A, vD, A, B, 16);
        VLMUL_1(A, vD, vD, A, B, 16);
        ref_vlmul_s16(A_ref, B);

        for(int i = 0; i < len; i++)
            TEST_ASSERT_EQUAL_INT16(A[i], A_ref[i]);

    }
}

TEST(inst_vlmul, vlmul_s16x16_sweep)
{

    int16_t A_ref[S16_LEN];
    int16_t A[S16_LEN];
    int16_t B[S16_LEN];
    int16_t vD[S16_LEN];

    for(int32_t a = INT16_MIN; a < INT16_MAX; a++){

        for(int32_t b = INT16_MIN; b < INT16_MAX;){

            for(int i = 0; i < S16_LEN; i++)
                A[i] = a;
            memcpy(A_ref, A, sizeof(A));

            for(int i = 0; i < S16_LEN; i++)
                B[i] = b++;

            VLMUL_0(A, vD, A, B, 16);
            VLMUL_1(A, vD, vD, A, B, 16);

            ref_vlmul_s16(A_ref, B);

            for(int i = 0; i < S16_LEN; i++)
                TEST_ASSERT_EQUAL_INT16(A[i], A_ref[i]);
        }
    }
}


// TEST(inst_vlmul, vlmul_s16xu16_basic)
// {
    
//     unsigned seed = SEED_FROM_FUNC_NAME();

//     int16_t A_ref[S16_LEN];
//     int16_t A[S16_LEN];
//     uint16_t B[S16_LEN];
//     int16_t vD[S16_LEN];

//     for(int v = 0; v < REPS; v++){

//         setExtraInfo_R(v);

//         const unsigned len = S16_LEN;

//         for(int i = 0; i < len; i++){
//             unsigned shr = pseudo_rand_uint32(&seed) % 5;
//             A[i] = pseudo_rand_int16(&seed) >> shr;
//         }

//         memcpy(A_ref, A, sizeof(A));

//         for(int i = 0; i < len; i++){
//             unsigned shr = pseudo_rand_uint32(&seed) % 5;
//             B[i] = pseudo_rand_uint16(&seed) >> shr;
//         }

//         VLMUL_0(A, vD, A, B, 16);
//         VLMUL_3(A, vD, vD, A, B, 16);
//         ref_vlmul_s16xu16(A_ref, B);

//         for(int i = 0; i < len; i++)
//             TEST_ASSERT_INT16_WITHIN(1, A[i], A_ref[i]);

//     }
// }

TEST(inst_vlmul, vlmul_s32x32_basic)
{
    
    unsigned seed = SEED_FROM_FUNC_NAME();

    int32_t A_ref[S32_LEN];
    int32_t A[S32_LEN];
    int32_t B[S32_LEN];
    int32_t vD[S32_LEN];

    for(int v = 0; v < REPS; v++){

        setExtraInfo_R(v);

        const unsigned len = S32_LEN;

        for(int i = 0; i < len; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            A[i] = pseudo_rand_int32(&seed) >> shr;
        }

        memcpy(A_ref, A, sizeof(A));

        for(int i = 0; i < len; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            B[i] = pseudo_rand_int32(&seed) >> shr;
        }

        VLMUL_0(A, vD, A, B, 32);
        VLMUL_1(A, vD, vD, A, B, 32);
        VLMUL_2(A, vD, vD, A, B, 32);
        VLMUL_3(A, vD, vD, A, B, 32);
        ref_vlmul_s32(A_ref, B);

        for(int i = 0; i < len; i++)
            TEST_ASSERT_EQUAL_INT32( A[i], A_ref[i]);

    }
}

