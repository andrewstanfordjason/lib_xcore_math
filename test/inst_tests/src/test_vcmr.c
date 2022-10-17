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

TEST_GROUP_RUNNER(inst_vcmr) {
  RUN_TEST_CASE(inst_vcmr, vcmr_s16_basic);
  RUN_TEST_CASE(inst_vcmr, vcmr_s32_basic);
}

TEST_GROUP(inst_vcmr);
TEST_SETUP(inst_vcmr) { fflush(stdout); }
TEST_TEAR_DOWN(inst_vcmr) {}

static char msg_buff[200];

#define TEST_ASSERT_EQUAL_MSG(EXPECTED, ACTUAL, LINE_NUM)   do{       \
    if((EXPECTED)!=(ACTUAL)) {                                        \
      sprintf(msg_buff, "(test vector @ line %u)", (LINE_NUM));       \
      TEST_ASSERT_EQUAL_MESSAGE((EXPECTED), (ACTUAL), msg_buff);      \
    }} while(0)


#if SMOKE_TEST
#  define REPS       (100)
#else
#  define REPS       (1000)
#endif

void VCMR_0(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode);
void VCMR_1(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in, 
    const int8_t * vR_in, int mode);
void VCMR_2(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in, 
    const int8_t * vR_in, int mode);
void VCMR_3(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in, 
    const int8_t * vR_in, int mode);

void ref_vcmr_s16(int16_t * accu, int16_t * vC, int16_t * mem){
    for(int i=0;i<8;i++){

        // int32_t v = (int32_t)(vC[i]) * (int32_t)(mem[i]);
        
        // int64_t a = (int64_t)accu[i] + v;

        // if (a > INT32_MAX) accu[i] = INT32_MAX;
        // else if (a < INT32_MIN) accu[i] = INT32_MIN;
        // else accu[i] = a;

    }
}

void ref_vcmr_s32(int32_t * accu, int32_t * vC, int32_t * mem){
    for(int i=0;i<4;i++){

    }
}

#define S8_LEN 32
#define S16_LEN 16
#define S32_LEN 8

TEST(inst_vcmr, vcmr_s16_basic)
{
    
    unsigned seed = SEED_FROM_FUNC_NAME();

    int32_t accu_ref[S16_LEN]; 
    int16_t vR[S16_LEN];    
    int16_t vD[S16_LEN];

    int16_t vC[S16_LEN];
    int16_t Mem[S16_LEN];

    memset(accu_ref, 0, sizeof(accu_ref));

    int32_t error = 0;
    int32_t abs_error = 0;

    for(int v = 0; v < REPS; v++){

        setExtraInfo_R(v);
        for(int i = 0; i < S16_LEN; i++){
            int32_t accu = accu_ref[i];
            ((int16_t*)vR)[i] = ((int16_t*)&accu)[0];
            ((int16_t*)vD)[i] = ((int16_t*)&accu)[1];
        }

        for(int i = 0; i < S16_LEN; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 1;
            vC[i] = pseudo_rand_int16(&seed) >> shr;
        }

        for(int i = 0; i < S16_LEN; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 1;
            Mem[i] = pseudo_rand_int16(&seed) >> shr;
        }

        // ref_vcmr_s16(accu_ref, vC, Mem);

        // vcmr_0(vD, vR, vD, vR, Mem, vC, 16);
        // vcmr_1(vD, vR, vD, vR, Mem, vC, 16);

        for(int i = 0; i < S16_LEN; i++){
            int32_t accu = 0;
            ((int16_t*)&accu)[0] = ((int16_t*)vR)[i];
            ((int16_t*)&accu)[1] = ((int16_t*)vD)[i];
            TEST_ASSERT_INT32_WITHIN(1, accu, accu_ref[i]);
            int e = accu_ref[i] - accu;
            error += e;
            if (e<0) e=-e;
            abs_error += e;
        }
    }


    // printf("error:%f abs_error:%f\n", (double)error / (REPS*16), (double)abs_error / (REPS*16));
}


TEST(inst_vcmr, vcmr_s32_basic)
{
    
    unsigned seed = SEED_FROM_FUNC_NAME();

    int64_t accu_ref[S32_LEN]; 
    int32_t vR[S32_LEN];    
    int32_t vD[S32_LEN];

    int32_t vC[S32_LEN];
    int32_t Mem[S32_LEN];

    memset(accu_ref, 0, sizeof(accu_ref));

    for(int v = 0; v < REPS; v++){

        setExtraInfo_R(v);


        for(int i = 0; i < S32_LEN; i++){
            int64_t accu = accu_ref[i];
            ((int32_t*)vR)[i] = ((int32_t*)&accu)[0];
            ((int32_t*)vD)[i] = ((int32_t*)&accu)[1];
        }

        for(int i = 0; i < S32_LEN; i++){
            unsigned shr = 20 + pseudo_rand_uint32(&seed) % 5;
            vC[i] = pseudo_rand_int32(&seed) >> shr;
        }

        for(int i = 0; i < S32_LEN; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            Mem[i] = pseudo_rand_int32(&seed) >> shr;
        }

        // ref_vcmr_s32(accu_ref, vC, Mem);

        // vcmr_0(vD, vR, vD, vR, Mem, vC, 32);
        // vcmr_1(vD, vR, vD, vR, Mem, vC, 32);
        // vcmr_2(vD, vR, vD, vR, Mem, vC, 32);
        // vcmr_3(vD, vR, vD, vR, Mem, vC, 32);

        for(int i = 0; i < S32_LEN; i++){
            int64_t accu = 0;
            ((int32_t*)&accu)[0] = ((int32_t*)vR)[i];
            ((int32_t*)&accu)[1] = ((int32_t*)vD)[i];
            TEST_ASSERT_INT64_WITHIN(2, accu, accu_ref[i]);
        }
    }
}

