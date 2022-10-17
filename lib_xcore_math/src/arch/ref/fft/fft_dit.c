// Copyright 2020-2022 XMOS LIMITED.
// This Software is subject to the terms of the XMOS Public Licence: Version 1.

#include <stdint.h>
#include <stdio.h>

#include "xmath/xmath.h"
#include "../../../vect/vpu_helper.h"
#include "xmath_fft_lut.h"

#include "emu.h"

//load 4 complex 32-bit values into a buffer
static void load_vec(
    complex_s32_t dst[], 
    const complex_s32_t src[])
{
    for(int i = 0; i < 4; i++){
        dst[i] = src[i];
    }
}

#include<string.h>
static void load_vec_16(
    complex_s16_t dst[], 
    const complex_s16_t src[])
{
    memcpy(dst, src, sizeof(complex_s16_t) * 8);
}

#include<string.h>
static void load_vec_16_half(
    complex_s16_t dst[], 
    const complex_s16_t src[])
{
    memcpy(dst, src, sizeof(complex_s16_t) * 4);
}

void fft_dit_forward_quake_s16 (
    complex_s16_t x[], 
    const unsigned N, 
    headroom_t* hr, 
    exponent_t* exp)
{
    const unsigned FFT_N_LOG2 = 31 - CLS_S32(N);

    const complex_s16_t* W = xmath_dit_fft_lut_16;

    exponent_t exp_modifier = 0;

    right_shift_t shift_mode = 0;

    complex_s16_t vD[8] = {{0}}, vR[8] = {{0}}, vC[8] = {{0}};

    shift_mode = (*hr == 3)? 0 : (*hr < 3)? 1 : -1;
    exp_modifier += shift_mode;

    for(int j = 0; j < (N>>3); j++){
        load_vec_16(vD, &x[8*j]);
        VFTTF(vD, 16, shift_mode);
        load_vec_16(&x[8*j], vD);
    }

    // print_vector_16("x: ", x);

    // this always happens
    if (N >= 8){
        //do the special vadsb

        const headroom_t cur_hr = vect_complex_s16_headroom_quake(x, N);
        
        shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
        exp_modifier += shift_mode;

        load_vec_16(vC, W);
        W = &W[8];

        for(int k = 0; k<N; k+=8){
            load_vec_16(vD, &x[k]);

            // vect_complex_s16_mul_quake(vD, vD, vC, 8, 0, 0);
            VCMR_0(vR, vD, vC, 16);
            VCMR_1(vR, vD, vC, vR, 16);

            VCMI_0(vR, vD, vC, 16);
            VCMI_1(vR, vD, vC, vR, 16);

            complex_s16_t T[24];
            load_vec_16(T, vR);
            load_vec_16(vR, &T[4]);
            VLADSB(vD, vR, T, vR, 16, shift_mode);
            load_vec_16(&x[k], vR);
            load_vec_16_half(&x[k+4], vD);
        }

    }

    if(N > 8){

        // int a = N >> 3;

        for(int n = 0; n < FFT_N_LOG2-3; n++){
            
            int b = 1<<(n+3);
            int a = 1<<((FFT_N_LOG2-4)-n);

            headroom_t cur_hr = vect_complex_s16_headroom_quake(x, N);

            shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
            exp_modifier += shift_mode;

            for(int k = b-8; k >= 0; k -= 8){

                int s = k;

                load_vec_16(vC, W);
                W = &W[8];

                for(int j = 0; j < a; j++){
                    load_vec_16(vD, &x[s+b]);

                    // vect_complex_s16_mul_quake(vR, vD, vC, 8, 0, 0);
                    VCMR_0(vR, vD, vC, 16);
                    VCMR_1(vR, vD, vC, vR, 16);

                    VCMI_0(vR, vD, vC, 16);
                    VCMI_1(vR, vD, vC, vR, 16);

                    for(int i = 0; i < 8; i++){
                        vD[i].re = ASHR(16)(((int32_t)x[s+i].re) - vR[i].re, shift_mode);
                        vD[i].im = ASHR(16)(((int32_t)x[s+i].im) - vR[i].im, shift_mode);
                        vR[i].re = ASHR(16)(((int32_t)x[s+i].re) + vR[i].re, shift_mode);
                        vR[i].im = ASHR(16)(((int32_t)x[s+i].im) + vR[i].im, shift_mode);
                    }

                    load_vec_16(&x[s], vR);
                    load_vec_16(&x[s+b], vD);
                    s += 2*b;
                };
            }
        }
    }

    *hr = vect_complex_s16_headroom_quake(x, N);
    *exp = *exp + exp_modifier;
}

void fft_dit_forward (
    complex_s32_t x[], 
    const unsigned N, 
    headroom_t* hr, 
    exponent_t* exp)
{
    const unsigned FFT_N_LOG2 = 31 - CLS_S32(N);

    const complex_s32_t* W = xmath_dit_fft_lut_32;

    exponent_t exp_modifier = 0;

    right_shift_t shift_mode = 0;

    complex_s32_t vD[4] = {{0}}, vR[4] = {{0}}, vC[4] = {{0}};

    shift_mode = (*hr == 3)? 0 : (*hr < 3)? 1 : -1;
    exp_modifier += shift_mode;


    for(int j = 0; j < (N>>2); j++){
        load_vec(vD, &x[4*j]);
        VFTTF(vD, 32, shift_mode);
        load_vec(&x[4*j], vD);
    }

    if(N != 4){

        // int a = N >> 3;

        for(int n = 0; n < FFT_N_LOG2-2; n++){
            
            int b = 1<<(n+2);
            int a = 1<<((FFT_N_LOG2-3)-n);

            headroom_t cur_hr = vect_complex_s32_headroom(x, N);

            shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
            exp_modifier += shift_mode;

            for(int k = b-4; k >= 0; k -= 4){

                int s = k;

                load_vec(vC, W);
                W = &W[4];

                for(int j = 0; j < a; j++){
                    load_vec(vD, &x[s+b]);

                    // vect_complex_s32_mul(vR, vD, vC, 4, 0, 0);
                    VCMR_0(vR, vD, vC, 32);
                    VCMR_1(vR, vD, vC, vR, 32);
                    VCMR_2(vR, vD, vC, vR, 32);
                    VCMR_3(vR, vD, vC, vR, 32);

                    VCMI_0(vR, vD, vC, 32);
                    VCMI_1(vR, vD, vC, vR, 32);
                    VCMI_2(vR, vD, vC, vR, 32);
                    VCMI_3(vR, vD, vC, vR, 32);

                    for(int i = 0; i < 4; i++){
                        vD[i].re = ASHR(32)(((int64_t)x[s+i].re) - vR[i].re, shift_mode);
                        vD[i].im = ASHR(32)(((int64_t)x[s+i].im) - vR[i].im, shift_mode);
                        vR[i].re = ASHR(32)(((int64_t)x[s+i].re) + vR[i].re, shift_mode);
                        vR[i].im = ASHR(32)(((int64_t)x[s+i].im) + vR[i].im, shift_mode);
                    }

                    load_vec(&x[s], vR);
                    load_vec(&x[s+b], vD);
                    s += 2*b;
                };
            }
        }
    }

    *hr = vect_complex_s32_headroom(x, N);
    *exp = *exp + exp_modifier;
}

void fft_dit_inverse (
    complex_s32_t x[], 
    const unsigned N, 
    headroom_t* hr, 
    exponent_t* exp)
{
    const unsigned FFT_N_LOG2 = 31 - CLS_S32(N);

    const complex_s32_t* W = xmath_dit_fft_lut_32;

    exponent_t exp_modifier = 0;

    right_shift_t shift_mode = 0;

    complex_s32_t vD[4] = {{0}}, vR[4] = {{0}}, vC[4] = {{0}};

    shift_mode = (*hr == 3)? 0 : (*hr < 3)? 1 : -1;
    exp_modifier += shift_mode;
    exp_modifier += -2;

    for(int j = 0; j < (N>>2); j++){
        load_vec(vD, &x[4*j]);
        VFTTB(vD, 32, shift_mode);
        load_vec(&x[4*j], vD);
    }

    if(N != 4){

        for(int n = 0; n < FFT_N_LOG2-2; n++){
            
            int b = 1<<(n+2);
            int a = 1<<((FFT_N_LOG2-3)-n);

            headroom_t cur_hr = vect_complex_s32_headroom(x, N);

            shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
            exp_modifier += shift_mode;
            exp_modifier += -1;

            for(int k = b-4; k >= 0; k -= 4){

                int s = k;

                load_vec(vC, W);
                W = &W[4];

                for(int j = 0; j < a; j++){
                    load_vec(vD, &x[s+b]);

                    // vect_complex_s32_conj_mul(vR, vD, vC, 4, 0, 0);
                    VCMCR_0(vR, vD, vC, 32);
                    VCMCR_1(vR, vD, vC, vR, 32);
                    VCMCR_2(vR, vD, vC, vR, 32);
                    VCMCR_3(vR, vD, vC, vR, 32);

                    VCMCI_0(vR, vD, vC, 32);
                    VCMCI_1(vR, vD, vC, vR, 32);
                    VCMCI_2(vR, vD, vC, vR, 32);
                    VCMCI_3(vR, vD, vC, vR, 32);


                    for(int i = 0; i < 4; i++){
                        vD[i].re = ASHR(32)(((int64_t)x[s+i].re) - vR[i].re, shift_mode);
                        vD[i].im = ASHR(32)(((int64_t)x[s+i].im) - vR[i].im, shift_mode);
                        vR[i].re = ASHR(32)(((int64_t)x[s+i].re) + vR[i].re, shift_mode);
                        vR[i].im = ASHR(32)(((int64_t)x[s+i].im) + vR[i].im, shift_mode);
                    }

                    load_vec(&x[s], vR);
                    load_vec(&x[s+b], vD);
                    s += 2*b;
                };
            }
        }
    }

    *hr = vect_complex_s32_headroom(x, N);
    *exp = *exp + exp_modifier;
}

#include<stdlib.h>

void fft_dit_inverse_quake_s16 (
    complex_s16_t x[], 
    const unsigned N, 
    headroom_t* hr, 
    exponent_t* exp)
{
    const unsigned FFT_N_LOG2 = 31 - CLS_S32(N);

    const complex_s16_t* W = xmath_dit_fft_lut_16;

    exponent_t exp_modifier = 0;

    right_shift_t shift_mode = 0;

    complex_s16_t vD[8] = {{0}}, vR[8] = {{0}}, vC[8] = {{0}};

    shift_mode = (*hr == 3)? 0 : (*hr < 3)? 1 : -1;
    exp_modifier += shift_mode;
    exp_modifier += -2;

    for(int j = 0; j < (N>>3); j++){
        load_vec_16(vD, &x[8*j]);
        VFTTB(vD, 16, shift_mode);
        load_vec_16(&x[8*j], vD);
    }
    
    // this always happens
    if (N >= 8){
        //do the special vadsb

        const headroom_t cur_hr = vect_complex_s16_headroom_quake(x, N);
        
        shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
        exp_modifier += shift_mode;
        exp_modifier += -1;

        load_vec_16(vC, W);
        W = &W[8];

        for(int k = 0; k<N; k+=8){
            load_vec_16(vD, &x[k]);

            // vect_complex_s16_conj_mul_quake(vD, vD, vC, 8, 0, 0);
            VCMCR_0(vR, vD, vC, 16);
            VCMCR_1(vR, vD, vC, vR, 16);

            VCMCI_0(vR, vD, vC, 16);
            VCMCI_1(vR, vD, vC, vR, 16);

            // VSBAD((int8_t *)vR, (const int8_t *)vR, shift_mode);
            // load_vec_16(&x[k], vR);

            //store vR in T
            complex_s16_t T[24];
            load_vec_16(T, vR);
            load_vec_16(vR, &T[4]);
            VLADSB(vD, vR, T, vR, 16, shift_mode);
            load_vec_16(&x[k], vR);
            load_vec_16_half(&x[k+4], vD);

            // printf("\n");
            // // load_vec_16(vR, &x[k]);
            // for (int i=0;i<8;i++){
            //     printf("%d %d %d\n", i, vD[i].re, vD[i].im);
            // }
            // printf("\n");
            // for (int i=0;i<8;i++){
            //     printf("%d %d %d\n", i, vR[i].re, vR[i].im);
            // }
            // exit(1);



        }
    } 

    if(N  > 8 ){

        for(int n = 0; n < FFT_N_LOG2-3; n++){
            
            int b = 1<<(n+3);
            int a = 1<<((FFT_N_LOG2-4)-n);

            headroom_t cur_hr = vect_complex_s16_headroom_quake(x, N);

            shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
            exp_modifier += shift_mode;
            exp_modifier += -1;

            for(int k = b-8; k >= 0; k -= 8){

                int s = k;

                load_vec_16(vC, W);
                W = &W[8];

                for(int j = 0; j < a; j++){
                    load_vec_16(vD, &x[s+b]);

                    // vect_complex_s16_conj_mul_quake(vR, vD, vC, 8, 0, 0);
                    VCMCR_0(vR, vD, vC, 16);
                    VCMCR_1(vR, vD, vC, vR, 16);

                    VCMCI_0(vR, vD, vC, 16);
                    VCMCI_1(vR, vD, vC, vR, 16);

                    for(int i = 0; i < 8; i++){
                        vD[i].re = ASHR(16)(((int32_t)x[s+i].re) - vR[i].re, shift_mode);
                        vD[i].im = ASHR(16)(((int32_t)x[s+i].im) - vR[i].im, shift_mode);
                        vR[i].re = ASHR(16)(((int32_t)x[s+i].re) + vR[i].re, shift_mode);
                        vR[i].im = ASHR(16)(((int32_t)x[s+i].im) + vR[i].im, shift_mode);
                    }

                    load_vec_16(&x[s], vR);
                    load_vec_16(&x[s+b], vD);
                    s += 2*b;
                };
            }
        }
    }

    *hr = vect_complex_s16_headroom_quake(x, N);
    *exp = *exp + exp_modifier;
}

