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

TEST_GROUP_RUNNER(inst_vlmaccr) {
  RUN_TEST_CASE(inst_vlmaccr, vlmaccr_s8xs8_basic);
  RUN_TEST_CASE(inst_vlmaccr, vlmaccr_s16xs16_basic);
  RUN_TEST_CASE(inst_vlmaccr, vlmaccr_s32xs32_basic);
  RUN_TEST_CASE(inst_vlmaccr, vlmaccr_s8xs8_patterns);
  RUN_TEST_CASE(inst_vlmaccr, vlmaccr_s16xs16_patterns);
  RUN_TEST_CASE(inst_vlmaccr, vlmaccr_s32xs32_patterns);
}

TEST_GROUP(inst_vlmaccr);
TEST_SETUP(inst_vlmaccr) { fflush(stdout); }
TEST_TEAR_DOWN(inst_vlmaccr) {}

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

void VLMACCR_0(
    int8_t * vD_out, int8_t * vR_out, 
    const int8_t * vD_in, const int8_t * vR_in, 
    const int8_t * vC_in, const int8_t * Mem_in, int mode);

void VLMACCR_1(
    int8_t * vD_out, int8_t * vR_out, 
    const int8_t * vD_in, const int8_t * vR_in, 
    const int8_t * vC_in, const int8_t * Mem_in, int mode);

void VLMACCR_2(
    int8_t * vD_out, int8_t * vR_out, 
    const int8_t * vD_in, const int8_t * vR_in, 
    const int8_t * vC_in, const int8_t * Mem_in, int mode);

void VLMACCR_3(
    int8_t * vD_out, int8_t * vR_out, 
    const int8_t * vD_in, const int8_t * vR_in, 
    const int8_t * vC_in, const int8_t * Mem_in, int mode);


void ref_vlmaccr_s8(int32_t * accu, int8_t * vC, int8_t * mem){
    
    int64_t sum = accu[15];

    for(int i=0;i<32;i++){
        int16_t v = (int16_t)(vC[i]) * (int16_t)(mem[i]);
        sum += v;
    }

    for(int i=15;i>=1;i--)
        accu[i] = accu[i-1];

    if (sum > INT32_MAX) accu[0] = INT32_MAX;
    else if (sum < INT32_MIN) accu[0] = INT32_MIN;
    else accu[0] = sum;

}

void ref_vlmaccr_s16(int32_t * accu, int16_t * vC, int16_t * mem){
    int64_t sum = accu[15];
    
    for(int i=0;i<16;i++){
        sum += (int32_t)(vC[i]) * (int32_t)(mem[i]);
    }

    for(int i=15;i>=1;i--)
        accu[i] = accu[i-1];

    if (sum > INT32_MAX) accu[0] = INT32_MAX;
    else if (sum < INT32_MIN) accu[0] = INT32_MIN;
    else accu[0] = sum;

}
void ref_vlmaccr_s16xs8(int32_t * accu, int16_t * vC, int8_t * mem){
    int64_t sum = accu[15];
    
    for(int i=0;i<16;i++){
        sum += (int32_t)(vC[i]) * (int32_t)(mem[2*i]);
    }

    for(int i=15;i>=1;i--)
        accu[i] = accu[i-1];

    if (sum > INT32_MAX) accu[0] = INT32_MAX;
    else if (sum < INT32_MIN) accu[0] = INT32_MIN;
    else accu[0] = sum;

}

void ref_vlmaccr_s32(int64_t * accu, int32_t * vC, int32_t * mem){
    __int128 sum = accu[7];
    
    for(int i=0;i<8;i++){
        sum += (__int128)(vC[i]) * (__int128)(mem[i]);
    }

    for(int i=7;i>=1;i--)
        accu[i] = accu[i-1];

    if (sum > INT64_MAX) accu[0] = INT64_MAX;
    else if (sum < INT64_MIN) accu[0] = INT64_MIN;
    else accu[0] = sum;
}

#define S8_LEN 32
#define S16_LEN 16
#define S32_LEN 8


TEST(inst_vlmaccr, vlmaccr_s8xs8_basic)
{
    
    unsigned seed = SEED_FROM_FUNC_NAME();

    int32_t accu_ref[S16_LEN]; 
    int16_t vR[S16_LEN];    
    int16_t vD[S16_LEN];

    int8_t vC[S8_LEN];
    int8_t Mem[S8_LEN];

    memset(accu_ref, 0, sizeof(accu_ref));
    memset(vR, 0, sizeof(vR));
    memset(vD, 0, sizeof(vD));

    for(int v = 0; v < REPS; v++){

        setExtraInfo_R(v);

        for(int i = 0; i < S8_LEN; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            vC[i] = pseudo_rand_int8(&seed) >> shr;
        }

        for(int i = 0; i < S8_LEN; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            Mem[i] = pseudo_rand_int8(&seed) >> shr;
        }

        ref_vlmaccr_s8(accu_ref, vC, Mem);

        VLMACCR_0(vD, vR, vD, vR, Mem, vC, 8);

        for(int i = 0; i < S16_LEN; i++){
            int32_t accu = 0;
            ((int16_t*)&accu)[0] = ((int16_t*)vR)[i];
            ((int16_t*)&accu)[1] = ((int16_t*)vD)[i];
            TEST_ASSERT_EQUAL_INT32(accu, accu_ref[i]);
        }

    }
}

TEST(inst_vlmaccr, vlmaccr_s16xs16_basic)
{
     
    unsigned seed = SEED_FROM_FUNC_NAME();

    int32_t accu_ref[S16_LEN]; 
    int16_t vR[S16_LEN];    
    int16_t vD[S16_LEN];

    int16_t vC[S16_LEN];
    int16_t Mem[S16_LEN];

    memset(accu_ref, 0, sizeof(accu_ref));
    memset(vR, 0, sizeof(vR));
    memset(vD, 0, sizeof(vD));

    for(int v = 0; v < REPS; v++){

        setExtraInfo_R(v);

        for(int i = 0; i < S16_LEN; i++){
            int32_t accu = accu_ref[i];
            ((int16_t*)vR)[i] = ((int16_t*)&accu)[0];
            ((int16_t*)vD)[i] = ((int16_t*)&accu)[1];
        }

        for(int i = 0; i < S16_LEN; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            vC[i] = pseudo_rand_int16(&seed) >> shr;
        }

        for(int i = 0; i < S16_LEN; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            Mem[i] = pseudo_rand_int16(&seed) >> shr;
        }

        ref_vlmaccr_s16(accu_ref, vC, Mem);

        VLMACCR_0(vD, vR, vD, vR, Mem, vC, 16);
        VLMACCR_1(vD, vR, vD, vR, Mem, vC, 16);

        for(int i = 0; i < S16_LEN; i++){
            int32_t accu = 0;
            ((int16_t*)&accu)[0] = ((int16_t*)vR)[i];
            ((int16_t*)&accu)[1] = ((int16_t*)vD)[i];
            TEST_ASSERT_INT32_WITHIN(1, accu, accu_ref[i]);
        }

    }
}

TEST(inst_vlmaccr, vlmaccr_s32xs32_basic)
{
    unsigned seed = SEED_FROM_FUNC_NAME();

    int64_t accu_ref[S32_LEN]; 
    int32_t vR[S32_LEN];    
    int32_t vD[S32_LEN];

    int32_t vC[S32_LEN];
    int32_t Mem[S32_LEN];

    memset(accu_ref, 0, sizeof(accu_ref));
    memset(vR, 0, sizeof(vR));
    memset(vD, 0, sizeof(vD));

    for(int v = 0; v < REPS; v++){

        setExtraInfo_R(v);

        for(int i = 0; i < S32_LEN; i++){
            int64_t accu = accu_ref[i];
            ((int32_t*)vR)[i] = ((int32_t*)&accu)[0];
            ((int32_t*)vD)[i] = ((int32_t*)&accu)[1];
        }

        for(int i = 0; i < S32_LEN; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            vC[i] = pseudo_rand_int32(&seed) >> shr;
        }

        for(int i = 0; i < S32_LEN; i++){
            unsigned shr = pseudo_rand_uint32(&seed) % 5;
            Mem[i] = pseudo_rand_int32(&seed) >> shr;
        }

        ref_vlmaccr_s32(accu_ref, vC, Mem);

        VLMACCR_0(vD, vR, vD, vR, Mem, vC, 32);
        VLMACCR_1(vD, vR, vD, vR, Mem, vC, 32);
        VLMACCR_2(vD, vR, vD, vR, Mem, vC, 32);
        VLMACCR_3(vD, vR, vD, vR, Mem, vC, 32);

        for(int i = 0; i < S32_LEN; i++){
            int64_t accu = 0;
            ((int32_t*)&accu)[0] = ((int32_t*)vR)[i];
            ((int32_t*)&accu)[1] = ((int32_t*)vD)[i];
            TEST_ASSERT_INT64_WITHIN(1, accu, accu_ref[i]);
        }

    }
}


TEST(inst_vlmaccr, vlmaccr_s8xs8_patterns)
{
    
    unsigned seed = SEED_FROM_FUNC_NAME();

    int32_t accu_ref[S16_LEN]; 
    int16_t vR[S16_LEN];    
    int16_t vD[S16_LEN];

    int8_t vC[S8_LEN];
    int8_t Mem[S8_LEN];

#define S8_PAT_COUNT 7
    int8_t patterns[S8_PAT_COUNT] = {
        0, 1, -1, INT8_MAX, INT8_MIN, INT8_MAX-1, INT8_MIN+1
    }; 
    memset(accu_ref, 0, sizeof(accu_ref));
    memset(vR, 0, sizeof(vR));
    memset(vD, 0, sizeof(vD));

    for (int pc=0;pc<S8_PAT_COUNT;pc++){
        for (int pm=0;pm<S8_PAT_COUNT;pm++){
            for(int v = 0; v < 32; v++){

                setExtraInfo_R(v);
        
                for(int i = 0; i < S8_LEN; i++){
                    vC[i] = patterns[pc];
                    Mem[i] = patterns[pm];
                }

                ref_vlmaccr_s8(accu_ref, vC, Mem);

                VLMACCR_0(vD, vR, vD, vR, Mem, vC, 8);

                for(int i = 0; i < S16_LEN; i++){
                    int32_t accu = 0;
                    ((int16_t*)&accu)[0] = ((int16_t*)vR)[i];
                    ((int16_t*)&accu)[1] = ((int16_t*)vD)[i];
                    TEST_ASSERT_EQUAL_INT32(accu, accu_ref[i]);
                }
            }
        }
    }
}

TEST(inst_vlmaccr, vlmaccr_s16xs16_patterns)
{
     
    unsigned seed = SEED_FROM_FUNC_NAME();

    int32_t accu_ref[S16_LEN]; 
    int16_t vR[S16_LEN];    
    int16_t vD[S16_LEN];

    int16_t vC[S16_LEN];
    int16_t Mem[S16_LEN];

#define S16_PAT_COUNT 7
    int16_t patterns[S16_PAT_COUNT] = {
        0, 1, -1, INT16_MAX, INT16_MIN, INT16_MAX-1, INT16_MIN+1
    }; 
    memset(accu_ref, 0, sizeof(accu_ref));
    memset(vR, 0, sizeof(vR));
    memset(vD, 0, sizeof(vD));

    for (int pc=0;pc<S16_PAT_COUNT;pc++){
        for (int pm=0;pm<S16_PAT_COUNT;pm++){
            for(int v = 0; v < 32; v++){

                setExtraInfo_R(v);

                for(int i = 0; i < S16_LEN; i++){
                    int32_t accu = accu_ref[i];
                    ((int16_t*)vR)[i] = ((int16_t*)&accu)[0];
                    ((int16_t*)vD)[i] = ((int16_t*)&accu)[1];
                }

        
                for(int i = 0; i < S16_LEN; i++){
                    vC[i] = patterns[pc];
                    Mem[i] = patterns[pm];
                }


                for(int i = 0; i < S16_LEN; i++){
                    unsigned shr = pseudo_rand_uint32(&seed) % 5;
                    vC[i] = pseudo_rand_int16(&seed) >> shr;
                }

                for(int i = 0; i < S16_LEN; i++){
                    unsigned shr = pseudo_rand_uint32(&seed) % 5;
                    Mem[i] = pseudo_rand_int16(&seed) >> shr;
                }

                ref_vlmaccr_s16(accu_ref, vC, Mem);

                VLMACCR_0(vD, vR, vD, vR, Mem, vC, 16);
                VLMACCR_1(vD, vR, vD, vR, Mem, vC, 16);

                for(int i = 0; i < S16_LEN; i++){
                    int32_t accu = 0;
                    ((int16_t*)&accu)[0] = ((int16_t*)vR)[i];
                    ((int16_t*)&accu)[1] = ((int16_t*)vD)[i];
                    TEST_ASSERT_INT32_WITHIN(1, accu, accu_ref[i]);
                }

            }
        }
    }
}

TEST(inst_vlmaccr, vlmaccr_s32xs32_patterns)
{
    unsigned seed = SEED_FROM_FUNC_NAME();

    int64_t accu_ref[S32_LEN]; 
    int32_t vR[S32_LEN];    
    int32_t vD[S32_LEN];

    int32_t vC[S32_LEN];
    int32_t Mem[S32_LEN];

#define S32_PAT_COUNT 7
    int32_t patterns[S32_PAT_COUNT] = {
        0, 1, -1, INT32_MAX, INT32_MIN, INT32_MAX-1, INT32_MIN+1
    }; 
    memset(accu_ref, 0, sizeof(accu_ref));
    memset(vR, 0, sizeof(vR));
    memset(vD, 0, sizeof(vD));

    for (int pc=0;pc<S32_PAT_COUNT;pc++){
        for (int pm=0;pm<S32_PAT_COUNT;pm++){
            for(int v = 0; v < 32; v++){

                setExtraInfo_R(v);

                for(int i = 0; i < S32_LEN; i++){
                    int64_t accu = accu_ref[i];
                    ((int32_t*)vR)[i] = ((int32_t*)&accu)[0];
                    ((int32_t*)vD)[i] = ((int32_t*)&accu)[1];
                }
        
                for(int i = 0; i < S32_LEN; i++){
                    vC[i] = patterns[pc];
                    Mem[i] = patterns[pm];
                }

                for(int i = 0; i < S32_LEN; i++){
                    unsigned shr = pseudo_rand_uint32(&seed) % 5;
                    vC[i] = pseudo_rand_int32(&seed) >> shr;
                }

                for(int i = 0; i < S32_LEN; i++){
                    unsigned shr = pseudo_rand_uint32(&seed) % 5;
                    Mem[i] = pseudo_rand_int32(&seed) >> shr;
                }

                ref_vlmaccr_s32(accu_ref, vC, Mem);

                VLMACCR_0(vD, vR, vD, vR, Mem, vC, 32);
                VLMACCR_1(vD, vR, vD, vR, Mem, vC, 32);
                VLMACCR_2(vD, vR, vD, vR, Mem, vC, 32);
                VLMACCR_3(vD, vR, vD, vR, Mem, vC, 32);

                for(int i = 0; i < S32_LEN; i++){
                    int64_t accu = 0;
                    ((int32_t*)&accu)[0] = ((int32_t*)vR)[i];
                    ((int32_t*)&accu)[1] = ((int32_t*)vD)[i];
                    TEST_ASSERT_INT64_WITHIN(1, accu, accu_ref[i]);
                }
            }
        }
    }
}

