// Copyright 2020-2022 XMOS LIMITED.
// This Software is subject to the terms of the XMOS Public Licence: Version 1.

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

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

void fft_dif_forward (
    complex_s32_t x[], 
    const unsigned N, 
    headroom_t* hr, 
    exponent_t* exp)
{
    const unsigned FFT_N_LOG2 = 31 - CLS_S32(N);

    const complex_s32_t* W = XMATH_DIF_FFT_LUT(N);

    exponent_t exp_modifier = 0;
    right_shift_t shift_mode = 0;

    complex_s32_t vD[4] = {{0}}, vR[4] = {{0}}, vC[4] = {{0}};

    shift_mode = (*hr == 3)? 0 : (*hr < 3)? 1 : -1;
    exp_modifier += shift_mode;

    if(N != 4){

        for(int n = 0; n < FFT_N_LOG2-2; n++){
            
            const int b = 1<<(FFT_N_LOG2-1-n);
            const int a = 1<<(2+n);

            for(int k = b-4; k >= 0; k -= 4){
                
                load_vec(vC, W);
                W = &W[4];

                for(int j = 0; j < a/4; j+=1){

                    const int s = 2*j*b+k;

                    load_vec(vR, &x[s]);

                    VLADSB(vD, vR, &(x[s+b]), vR, 32, shift_mode);

                    load_vec(&x[s], vR);

                    // vect_complex_s32_mul(vR, vD, vC, 4, 0, 0);
                    VCMR_0(vR, vD, vC, 32);
                    VCMR_1(vR, vD, vC, vR, 32);
                    VCMR_2(vR, vD, vC, vR, 32);
                    VCMR_3(vR, vD, vC, vR, 32);

                    VCMI_0(vR, vD, vC, 32);
                    VCMI_1(vR, vD, vC, vR, 32);
                    VCMI_2(vR, vD, vC, vR, 32);
                    VCMI_3(vR, vD, vC, vR, 32);

                    load_vec(&x[s+b], vR);
                    
                };
            }
            
            const headroom_t cur_hr = vect_complex_s32_headroom(x, N);
            
            shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
            exp_modifier += shift_mode;
            
        }
    }


    for(int j = 0; j < (N>>2); j++){
        load_vec(vR, &x[4*j]);
        VFTFF(vR, 32, shift_mode);
        load_vec(&x[4*j], vR);
    }

    *hr = vect_complex_s32_headroom(x, N);
    *exp = *exp + exp_modifier;
}

void fft_dif_forward_quake_s16 (
    complex_s16_t x[], 
    const unsigned N, 
    headroom_t* hr, 
    exponent_t* exp)
{
    const unsigned FFT_N_LOG2 = 31 - CLS_S32(N);

    const complex_s16_t* W = XMATH_DIF_FFT_LUT_16(N);

    exponent_t exp_modifier = 0;
    right_shift_t shift_mode = 0;

    complex_s16_t vD[8] = {{0}}, vR[8] = {{0}}, vC[8] = {{0}};

    shift_mode = (*hr == 3)? 0 : (*hr < 3)? 1 : -1;
    exp_modifier += shift_mode;

    if(N > 8){

        for(int n = 0; n < FFT_N_LOG2-3; n++){

            const int b = 1<<(FFT_N_LOG2-1-n);
            const int a = 1<<(3+n);

            for(int k = b-8; k >= 0; k -= 8){

                load_vec_16(vC, W);
                W = &W[8];

                for(int j = 0; j < a/8; j+=1){

                    const int s = 2*j*b+k;         

                    load_vec_16(vR, &x[s]);

                    VLADSB(vD, vR, &(x[s+b]), vR, 16, shift_mode);

                    load_vec_16(&x[s], vR);
                    // vect_complex_s16_mul_quake(vR, vD, vC, 8, 0, 0);
                    VCMR_0(vR, vD, vC, 16);
                    VCMR_1(vR, vD, vC, vR, 16);

                    VCMI_0(vR, vD, vC, 16);
                    VCMI_1(vR, vD, vC, vR, 16);
                    
                    load_vec_16(&x[s+b], vR);

                };
            }

            const headroom_t cur_hr = vect_complex_s16_headroom_quake(x, N);
            
            shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
            exp_modifier += shift_mode;
            
        } 
    }

    // this always happens
    if (N >= 8){
        //do the special vadsb

        load_vec_16(vC, W);

        for(int k = 0; k<N; k+=8){

            // load_vec_16(vD, &x[k]);
            // VADSB((int8_t *)vD, (const int8_t *)vD, shift_mode);

            load_vec_16(vR, &x[k]);
            complex_s16_t T[24];
            load_vec_16(T, vR);
            VLADSB(vD, vR, &T[4], vR, 16, shift_mode);
            load_vec_16(T, vR);
            load_vec_16(&T[4], vD);
            load_vec_16(vD, T);

            VCMR_0(vR, vD, vC, 16);
            VCMR_1(vR, vD, vC, vR, 16);

            VCMI_0(vR, vD, vC, 16);
            VCMI_1(vR, vD, vC, vR, 16);
            load_vec_16(&x[k], vR);

        }
        
        const headroom_t cur_hr = vect_complex_s16_headroom_quake(x, N);
        
        shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
        exp_modifier += shift_mode;
    }

    for(int j = 0; j < (N>>3); j++){
        load_vec_16(vR, &x[8*j]);
        VFTFF(vR, 16, shift_mode);
        load_vec_16(&x[8*j], vR);
    }

    *hr = vect_complex_s16_headroom_quake(x, N);
    *exp = *exp + exp_modifier;
}

void fft_dif_inverse_quake_s16 (
    complex_s16_t x[], 
    const unsigned N, 
    headroom_t* hr, 
    exponent_t* exp)
{
    const unsigned FFT_N_LOG2 = 31 - CLS_S32(N);

    const complex_s16_t* W = XMATH_DIF_FFT_LUT_16(N);

    exponent_t exp_modifier = -FFT_N_LOG2;
    right_shift_t shift_mode = 0;

    complex_s16_t vD[8] = {{0}}, vR[8] = {{0}}, vC[8] = {{0}};

    shift_mode = (*hr == 3)? 0 : (*hr < 3)? 1 : -1;
    exp_modifier += shift_mode;

    if(N > 8){

        for(int n = 0; n < FFT_N_LOG2-3; n++){
            
            const int b = 1<<(FFT_N_LOG2-1-n);
            const int a = 1<<(3+n);

            for(int k = b-8; k >= 0; k -= 8){
                
                load_vec_16(vC, W);
                W = &W[8];

                for(int j = 0; j < a/8; j+=1){

                    const int s = 2*j*b+k;

                    load_vec_16(vR, &x[s]);

                    VLADSB(vD, vR, &(x[s+b]), vR, 16, shift_mode);

                    load_vec_16(&x[s], vR);

                    VCMCR_0(vR, vD, vC, 16);
                    VCMCR_1(vR, vD, vC, vR, 16);

                    VCMCI_0(vR, vD, vC, 16);
                    VCMCI_1(vR, vD, vC, vR, 16);

                    load_vec_16(&x[s+b], vR);
                    
                };
            }
            
            const headroom_t cur_hr = vect_complex_s16_headroom_quake(x, N);
            
            shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
            exp_modifier += shift_mode;
            
        }
    }
    
    // this always happens
    if (N >= 8){
        //do the special vadsb

        load_vec_16(vC, W);

        for(int k = 0; k<N; k+=8){

            // load_vec_16(vD, &x[k]);

            // VADSB((int8_t *)vD, (const int8_t *)vD, shift_mode);
            load_vec_16(vR, &x[k]);
            VLADSB(vD, vR, &x[k+4], vR, 16, shift_mode);
            
            complex_s16_t T[12];
            load_vec_16(T, vR);
            load_vec_16(&T[4], vD);
            load_vec_16(vD, T);

            VCMCR_0(vR, vD, vC, 16);
            VCMCR_1(vR, vD, vC, vR, 16);

            VCMCI_0(vR, vD, vC, 16);
            VCMCI_1(vR, vD, vC, vR, 16);

            load_vec_16(&x[k], vR);

        }
        
        const headroom_t cur_hr = vect_complex_s16_headroom_quake(x, N);
        
        shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
        exp_modifier += shift_mode;
    }

    for(int j = 0; j < (N>>3); j++){
        load_vec_16(vR, &x[8*j]);
        // vftfb_quake_s16(vR, shift_mode);
        VFTFB(vR, 16, shift_mode);
        load_vec_16(&x[8*j], vR);
    }


    *hr = vect_complex_s16_headroom_quake(x, N);
    *exp = *exp + exp_modifier;
}

void fft_dif_inverse (
    complex_s32_t x[], 
    const unsigned N, 
    headroom_t* hr, 
    exponent_t* exp)
{
    const unsigned FFT_N_LOG2 = 31 - CLS_S32(N);

    const complex_s32_t* W = XMATH_DIF_FFT_LUT(N);

    exponent_t exp_modifier = -FFT_N_LOG2;
    right_shift_t shift_mode = 0;

    complex_s32_t vD[4] = {{0}}, vR[4] = {{0}}, vC[4] = {{0}};

    shift_mode = (*hr == 3)? 0 : (*hr < 3)? 1 : -1;
    exp_modifier += shift_mode;

    if(N != 4){

        for(int n = 0; n < FFT_N_LOG2-2; n++){
            
            const int b = 1<<(FFT_N_LOG2-1-n);
            const int a = 1<<(2+n);

            for(int k = b-4; k >= 0; k -= 4){
                
                load_vec(vC, W);
                W = &W[4];

                for(int j = 0; j < a/4; j+=1){

                    const int s = 2*j*b+k;

                    load_vec(vR, &x[s]);

                    VLADSB(vD, vR, &(x[s+b]), vR, 32, shift_mode);

                    load_vec(&x[s], vR);

                    memset(vR, 0, 32);

                    VCMCR_0(vR, vD, vC, 32);
                    VCMCR_1(vR, vD, vC, vR, 32);
                    VCMCR_2(vR, vD, vC, vR, 32);
                    VCMCR_3(vR, vD, vC, vR, 32);

                    VCMCI_0(vR, vD, vC, 32);
                    VCMCI_1(vR, vD, vC, vR, 32);
                    VCMCI_2(vR, vD, vC, vR, 32);
                    VCMCI_3(vR, vD, vC, vR, 32);

                    load_vec(&x[s+b], vR);
                    
                };
            }
            
            const headroom_t cur_hr = vect_complex_s32_headroom(x, N);
            
            shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
            exp_modifier += shift_mode;
            
        }
    }
    

    for(int j = 0; j < (N>>2); j++){
        load_vec(vR, &x[4*j]);
        // vftfb(vR, shift_mode);
        VFTFB(vR, 32, shift_mode);
        load_vec(&x[4*j], vR);
    }

    *hr = vect_complex_s32_headroom(x, N);
    *exp = *exp + exp_modifier;
}