// Copyright 2020-2022 XMOS LIMITED.
// This Software is subject to the terms of the XMOS Public Licence: Version 1.

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "xmath/xmath.h"
#include "../../../vect/vpu_helper.h"
#include "xmath_fft_lut.h"

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

static void print_vector_16(const char * name, const complex_s16_t x[]){
    printf("%s ", name);
    for (int i=0;i<8;i++){
        printf("(%.2f+%.2fj) ", (float)x[i].re/(1<<14) ,  (float)x[i].im/(1<<14));
    }
    printf("\n");
}

static void vftff(
    complex_s32_t vR[],
    const right_shift_t shift_mode)
{
    struct {
        int64_t re;
        int64_t im;
    } s[4];

    s[0].re = vR[0].re + vR[2].re;
    s[0].im = vR[0].im + vR[2].im;
    s[1].re = vR[1].re + vR[3].re;
    s[1].im = vR[1].im + vR[3].im;
    s[2].re = vR[0].re - vR[2].re;
    s[2].im = vR[0].im - vR[2].im;
    s[3].re = vR[1].im - vR[3].im;
    s[3].im = vR[3].re - vR[1].re;

    vR[0].re = (int32_t) ASHR(32)(s[0].re + s[1].re, shift_mode);
    vR[0].im = (int32_t) ASHR(32)(s[0].im + s[1].im, shift_mode);
    vR[1].re = (int32_t) ASHR(32)(s[0].re - s[1].re, shift_mode);
    vR[1].im = (int32_t) ASHR(32)(s[0].im - s[1].im, shift_mode);
    vR[2].re = (int32_t) ASHR(32)(s[2].re + s[3].re, shift_mode);
    vR[2].im = (int32_t) ASHR(32)(s[2].im + s[3].im, shift_mode);
    vR[3].re = (int32_t) ASHR(32)(s[2].re - s[3].re, shift_mode);
    vR[3].im = (int32_t) ASHR(32)(s[2].im - s[3].im, shift_mode);
}

static void vftff_quake_s16(
    complex_s16_t vR[],
    const right_shift_t shift_mode)
{
    struct {
        int32_t re;
        int32_t im;
    } s[8];

    for (int i=0; i<8; i+=4){
        s[0+i].re = vR[0+i].re + vR[2+i].re;
        s[0+i].im = vR[0+i].im + vR[2+i].im;
        s[1+i].re = vR[1+i].re + vR[3+i].re;
        s[1+i].im = vR[1+i].im + vR[3+i].im;
        s[2+i].re = vR[0+i].re - vR[2+i].re;
        s[2+i].im = vR[0+i].im - vR[2+i].im;
        s[3+i].re = vR[1+i].im - vR[3+i].im;
        s[3+i].im = vR[3+i].re - vR[1+i].re;
    }

    for (int i=0; i<8; i+=4){
        vR[0+i].re = (int16_t) ASHR(16)(s[0+i].re + s[1+i].re, shift_mode);
        vR[0+i].im = (int16_t) ASHR(16)(s[0+i].im + s[1+i].im, shift_mode);
        vR[1+i].re = (int16_t) ASHR(16)(s[0+i].re - s[1+i].re, shift_mode);
        vR[1+i].im = (int16_t) ASHR(16)(s[0+i].im - s[1+i].im, shift_mode);
        vR[2+i].re = (int16_t) ASHR(16)(s[2+i].re + s[3+i].re, shift_mode);
        vR[2+i].im = (int16_t) ASHR(16)(s[2+i].im + s[3+i].im, shift_mode);
        vR[3+i].re = (int16_t) ASHR(16)(s[2+i].re - s[3+i].re, shift_mode);
        vR[3+i].im = (int16_t) ASHR(16)(s[2+i].im - s[3+i].im, shift_mode);
    }
}
static void vftfb(
    complex_s32_t vR[],
    const right_shift_t shift_mode)
{
    struct {
        int64_t re;
        int64_t im;
    } s[4];

    s[0].re = vR[0].re + vR[2].re;
    s[0].im = vR[0].im + vR[2].im;
    s[1].re = vR[1].re + vR[3].re;
    s[1].im = vR[1].im + vR[3].im;
    s[2].re = vR[0].re - vR[2].re;
    s[2].im = vR[0].im - vR[2].im;
    s[3].re = vR[3].im - vR[1].im;
    s[3].im = vR[1].re - vR[3].re;

    vR[0].re = (int32_t) ASHR(32)(s[0].re + s[1].re, shift_mode);
    vR[0].im = (int32_t) ASHR(32)(s[0].im + s[1].im, shift_mode);
    vR[1].re = (int32_t) ASHR(32)(s[0].re - s[1].re, shift_mode);
    vR[1].im = (int32_t) ASHR(32)(s[0].im - s[1].im, shift_mode);
    vR[2].re = (int32_t) ASHR(32)(s[2].re + s[3].re, shift_mode);
    vR[2].im = (int32_t) ASHR(32)(s[2].im + s[3].im, shift_mode);
    vR[3].re = (int32_t) ASHR(32)(s[2].re - s[3].re, shift_mode);
    vR[3].im = (int32_t) ASHR(32)(s[2].im - s[3].im, shift_mode);
}

static void vftfb_quake_s16(
    complex_s16_t vR[],
    const right_shift_t shift_mode)
{
    struct {
        int64_t re;
        int64_t im;
    } s[8];

    for (int i=0; i<8; i+=4){
        s[0+i].re = vR[0+i].re + vR[2+i].re;
        s[0+i].im = vR[0+i].im + vR[2+i].im;
        s[1+i].re = vR[1+i].re + vR[3+i].re;
        s[1+i].im = vR[1+i].im + vR[3+i].im;
        s[2+i].re = vR[0+i].re - vR[2+i].re;
        s[2+i].im = vR[0+i].im - vR[2+i].im;
        s[3+i].re = vR[3+i].im - vR[1+i].im;
        s[3+i].im = vR[1+i].re - vR[3+i].re;
    }

    for (int i=0; i<8; i+=4){
        vR[0+i].re = (int16_t) ASHR(16)(s[0+i].re + s[1+i].re, shift_mode);
        vR[0+i].im = (int16_t) ASHR(16)(s[0+i].im + s[1+i].im, shift_mode);
        vR[1+i].re = (int16_t) ASHR(16)(s[0+i].re - s[1+i].re, shift_mode);
        vR[1+i].im = (int16_t) ASHR(16)(s[0+i].im - s[1+i].im, shift_mode);
        vR[2+i].re = (int16_t) ASHR(16)(s[2+i].re + s[3+i].re, shift_mode);
        vR[2+i].im = (int16_t) ASHR(16)(s[2+i].im + s[3+i].im, shift_mode);
        vR[3+i].re = (int16_t) ASHR(16)(s[2+i].re - s[3+i].re, shift_mode);
        vR[3+i].im = (int16_t) ASHR(16)(s[2+i].im - s[3+i].im, shift_mode);
    }
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

                    for(int i = 0; i < 4; i++){
                        vD[i].re = ASHR(32)(((int64_t)x[s+b+i].re) - vR[i].re, shift_mode);
                        vD[i].im = ASHR(32)(((int64_t)x[s+b+i].im) - vR[i].im, shift_mode);
                        vR[i].re = ASHR(32)(((int64_t)x[s+b+i].re) + vR[i].re, shift_mode);
                        vR[i].im = ASHR(32)(((int64_t)x[s+b+i].im) + vR[i].im, shift_mode);
                    }

                    load_vec(&x[s], vR);

                    vect_complex_s32_mul(vR, vD, vC, 4, 0, 0);

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
        vftff(vR, shift_mode);
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

                    for(int i = 0; i < 8; i++){
                        vD[i].re = ASHR(16)(((int32_t)x[s+b+i].re) - vR[i].re, shift_mode);
                        vD[i].im = ASHR(16)(((int32_t)x[s+b+i].im) - vR[i].im, shift_mode);
                        vR[i].re = ASHR(16)(((int32_t)x[s+b+i].re) + vR[i].re, shift_mode);
                        vR[i].im = ASHR(16)(((int32_t)x[s+b+i].im) + vR[i].im, shift_mode);
                    }

                    load_vec_16(&x[s], vR);
                    vect_complex_s16_mul_quake(vR, vD, vC, 8, 0, 0);
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

            load_vec_16(vD, &x[k]);

            for(int i = 0; i < 4; i++){
                int16_t t_re = ASHR(16)(((int32_t)vD[i+4].re) + vD[i].re, shift_mode);
                int16_t t_im = ASHR(16)(((int32_t)vD[i+4].im) + vD[i].im, shift_mode);
                vD[i+4].re   = ASHR(16)(((int32_t)vD[i+4].re) - vD[i].re, shift_mode);
                vD[i+4].im   = ASHR(16)(((int32_t)vD[i+4].im) - vD[i].im, shift_mode);
                vD[i].re = t_re;
                vD[i].im = t_im;
            }

            vect_complex_s16_mul_quake(vR, vD, vC, 8, 0, 0);
            load_vec_16(&x[k], vR);

        }
        
        const headroom_t cur_hr = vect_complex_s16_headroom_quake(x, N);
        
        shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
        exp_modifier += shift_mode;
    }

    for(int j = 0; j < (N>>3); j++){
        load_vec_16(vR, &x[8*j]);
        vftff_quake_s16(vR, shift_mode);
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

                    for(int i = 0; i < 8; i++){
                        vD[i].re = ASHR(16)(((int32_t)x[s+b+i].re) - vR[i].re, shift_mode);
                        vD[i].im = ASHR(16)(((int32_t)x[s+b+i].im) - vR[i].im, shift_mode);
                        vR[i].re = ASHR(16)(((int32_t)x[s+b+i].re) + vR[i].re, shift_mode);
                        vR[i].im = ASHR(16)(((int32_t)x[s+b+i].im) + vR[i].im, shift_mode);
                    }

                    load_vec_16(&x[s], vR);

                    vect_complex_s16_conj_mul_quake(vR, vD, vC, 8, 0, 0);

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

            load_vec_16(vD, &x[k]);

            for(int i = 0; i < 4; i++){
                int16_t t_re = ASHR(16)(((int32_t)vD[i+4].re) + vD[i].re, shift_mode);
                int16_t t_im = ASHR(16)(((int32_t)vD[i+4].im) + vD[i].im, shift_mode);
                vD[i+4].re   = ASHR(16)(((int32_t)vD[i+4].re) - vD[i].re, shift_mode);
                vD[i+4].im   = ASHR(16)(((int32_t)vD[i+4].im) - vD[i].im, shift_mode);
                vD[i].re = t_re;
                vD[i].im = t_im;
            }

            vect_complex_s16_conj_mul_quake(vR, vD, vC, 8, 0, 0);
            load_vec_16(&x[k], vR);

        }
        
        const headroom_t cur_hr = vect_complex_s16_headroom_quake(x, N);
        
        shift_mode = (cur_hr == 3)? 0 : (cur_hr < 3)? 1 : -1;
        exp_modifier += shift_mode;
    }

    for(int j = 0; j < (N>>3); j++){
        load_vec_16(vR, &x[8*j]);
        vftfb_quake_s16(vR, shift_mode);
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

                    for(int i = 0; i < 4; i++){
                        vD[i].re = ASHR(32)(((int64_t)x[s+b+i].re) - vR[i].re, shift_mode);
                        vD[i].im = ASHR(32)(((int64_t)x[s+b+i].im) - vR[i].im, shift_mode);
                        vR[i].re = ASHR(32)(((int64_t)x[s+b+i].re) + vR[i].re, shift_mode);
                        vR[i].im = ASHR(32)(((int64_t)x[s+b+i].im) + vR[i].im, shift_mode);
                    }

                    load_vec(&x[s], vR);

                    vect_complex_s32_conj_mul(vR, vD, vC, 4, 0, 0);

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
        vftfb(vR, shift_mode);
        load_vec(&x[4*j], vR);
    }

    *hr = vect_complex_s32_headroom(x, N);
    *exp = *exp + exp_modifier;
}