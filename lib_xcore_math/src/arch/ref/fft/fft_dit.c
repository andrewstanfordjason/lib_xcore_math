// Copyright 2020-2022 XMOS LIMITED.
// This Software is subject to the terms of the XMOS Public Licence: Version 1.

#include <stdint.h>
#include <stdio.h>

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

static void vfttf(
    complex_s32_t vD[],
    const right_shift_t shift_mode)
{
    struct {
        int64_t re;
        int64_t im;
    } s[4];

    s[0].re = vD[0].re + vD[1].re;
    s[0].im = vD[0].im + vD[1].im;
    s[1].re = vD[0].re - vD[1].re;
    s[1].im = vD[0].im - vD[1].im;
    s[2].re = vD[2].re + vD[3].re;
    s[2].im = vD[2].im + vD[3].im;
    s[3].re = vD[2].im - vD[3].im;
    s[3].im = vD[3].re - vD[2].re;

    vD[0].re = (int32_t) ASHR(32)(s[0].re + s[2].re, shift_mode);
    vD[0].im = (int32_t) ASHR(32)(s[0].im + s[2].im, shift_mode);
    vD[1].re = (int32_t) ASHR(32)(s[1].re + s[3].re, shift_mode);
    vD[1].im = (int32_t) ASHR(32)(s[1].im + s[3].im, shift_mode);
    vD[2].re = (int32_t) ASHR(32)(s[0].re - s[2].re, shift_mode);
    vD[2].im = (int32_t) ASHR(32)(s[0].im - s[2].im, shift_mode);
    vD[3].re = (int32_t) ASHR(32)(s[1].re - s[3].re, shift_mode);
    vD[3].im = (int32_t) ASHR(32)(s[1].im - s[3].im, shift_mode);
}


static void vfttf_quake_s16(
    complex_s16_t vD[],
    const right_shift_t shift_mode)
{
    struct {
        int32_t re;
        int32_t im;
    } s[8];

    for (int i=0; i<8; i+=4){
        s[0+i].re = vD[0+i].re + vD[1+i].re;
        s[0+i].im = vD[0+i].im + vD[1+i].im;
        s[1+i].re = vD[0+i].re - vD[1+i].re;
        s[1+i].im = vD[0+i].im - vD[1+i].im;
        s[2+i].re = vD[2+i].re + vD[3+i].re;
        s[2+i].im = vD[2+i].im + vD[3+i].im;
        s[3+i].re = vD[2+i].im - vD[3+i].im;
        s[3+i].im = vD[3+i].re - vD[2+i].re;
    }
    for (int i=0; i<8; i+=4){
        vD[0+i].re = (int16_t) ASHR(16)(s[0+i].re + s[2+i].re, shift_mode);
        vD[0+i].im = (int16_t) ASHR(16)(s[0+i].im + s[2+i].im, shift_mode);
        vD[1+i].re = (int16_t) ASHR(16)(s[1+i].re + s[3+i].re, shift_mode);
        vD[1+i].im = (int16_t) ASHR(16)(s[1+i].im + s[3+i].im, shift_mode);
        vD[2+i].re = (int16_t) ASHR(16)(s[0+i].re - s[2+i].re, shift_mode);
        vD[2+i].im = (int16_t) ASHR(16)(s[0+i].im - s[2+i].im, shift_mode);
        vD[3+i].re = (int16_t) ASHR(16)(s[1+i].re - s[3+i].re, shift_mode);
        vD[3+i].im = (int16_t) ASHR(16)(s[1+i].im - s[3+i].im, shift_mode);
    }
}


static void vfttb(
    complex_s32_t vD[],
    const right_shift_t shift_mode)
{
    struct {
        int64_t re;
        int64_t im;
    } s[4];

    s[0].re = vD[0].re + vD[1].re;
    s[0].im = vD[0].im + vD[1].im;
    s[1].re = vD[0].re - vD[1].re;
    s[1].im = vD[0].im - vD[1].im;
    s[2].re = vD[2].re + vD[3].re;
    s[2].im = vD[2].im + vD[3].im;
    s[3].re = vD[3].im - vD[2].im;
    s[3].im = vD[2].re - vD[3].re;

    vD[0].re = (int32_t) ASHR(32)(s[0].re + s[2].re, shift_mode);
    vD[0].im = (int32_t) ASHR(32)(s[0].im + s[2].im, shift_mode);
    vD[1].re = (int32_t) ASHR(32)(s[1].re + s[3].re, shift_mode);
    vD[1].im = (int32_t) ASHR(32)(s[1].im + s[3].im, shift_mode);
    vD[2].re = (int32_t) ASHR(32)(s[0].re - s[2].re, shift_mode);
    vD[2].im = (int32_t) ASHR(32)(s[0].im - s[2].im, shift_mode);
    vD[3].re = (int32_t) ASHR(32)(s[1].re - s[3].re, shift_mode);
    vD[3].im = (int32_t) ASHR(32)(s[1].im - s[3].im, shift_mode);
}

static void vfttb_quake_s16(
    complex_s16_t vD[],
    const right_shift_t shift_mode)
{
    struct {
        int32_t re;
        int32_t im;
    } s[8];

    for (int i=0; i<8; i+=4){
        s[0+i].re = vD[0+i].re + vD[1+i].re;
        s[0+i].im = vD[0+i].im + vD[1+i].im;
        s[1+i].re = vD[0+i].re - vD[1+i].re;
        s[1+i].im = vD[0+i].im - vD[1+i].im;
        s[2+i].re = vD[2+i].re + vD[3+i].re;
        s[2+i].im = vD[2+i].im + vD[3+i].im;
        s[3+i].re = vD[3+i].im - vD[2+i].im;
        s[3+i].im = vD[2+i].re - vD[3+i].re;
    }

    for (int i=0; i<8; i+=4){
        vD[0+i].re = (int16_t) ASHR(16)(s[0+i].re + s[2+i].re, shift_mode);
        vD[0+i].im = (int16_t) ASHR(16)(s[0+i].im + s[2+i].im, shift_mode);
        vD[1+i].re = (int16_t) ASHR(16)(s[1+i].re + s[3+i].re, shift_mode);
        vD[1+i].im = (int16_t) ASHR(16)(s[1+i].im + s[3+i].im, shift_mode);
        vD[2+i].re = (int16_t) ASHR(16)(s[0+i].re - s[2+i].re, shift_mode);
        vD[2+i].im = (int16_t) ASHR(16)(s[0+i].im - s[2+i].im, shift_mode);
        vD[3+i].re = (int16_t) ASHR(16)(s[1+i].re - s[3+i].re, shift_mode);
        vD[3+i].im = (int16_t) ASHR(16)(s[1+i].im - s[3+i].im, shift_mode);
    }
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
        vfttf_quake_s16(vD, shift_mode);
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

            vect_complex_s16_mul_quake(vD, vD, vC, 8, 0, 0);

            //FIXME this is backwards from the dif
            for(int i = 0; i < 4; i++){
                int16_t t_re = ASHR(16)(((int32_t)vD[i+4].re) + vD[i].re, shift_mode);
                int16_t t_im = ASHR(16)(((int32_t)vD[i+4].im) + vD[i].im, shift_mode);
                vD[i+4].re   = -ASHR(16)(((int32_t)vD[i+4].re) - vD[i].re, shift_mode);
                vD[i+4].im   = -ASHR(16)(((int32_t)vD[i+4].im) - vD[i].im, shift_mode);
                vD[i].re = t_re;
                vD[i].im = t_im;
            }

            load_vec_16(&x[k], vD);

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

                    vect_complex_s16_mul_quake(vR, vD, vC, 8, 0, 0);
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
        vfttf(vD, shift_mode);
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

                    vect_complex_s32_mul(vR, vD, vC, 4, 0, 0);

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
        vfttb(vD, shift_mode);
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

                    vect_complex_s32_conj_mul(vR, vD, vC, 4, 0, 0);

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
        vfttb_quake_s16(vD, shift_mode);
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

            vect_complex_s16_conj_mul_quake(vD, vD, vC, 8, 0, 0);

            //FIXME this is backwards from the dif
            for(int i = 0; i < 4; i++){
                int16_t t_re = ASHR(16)(((int32_t)vD[i+4].re) + vD[i].re, shift_mode);
                int16_t t_im = ASHR(16)(((int32_t)vD[i+4].im) + vD[i].im, shift_mode);
                vD[i+4].re   = -ASHR(16)(((int32_t)vD[i+4].re) - vD[i].re, shift_mode);
                vD[i+4].im   = -ASHR(16)(((int32_t)vD[i+4].im) - vD[i].im, shift_mode);
                vD[i].re = t_re;
                vD[i].im = t_im;
            }

            load_vec_16(&x[k], vD);
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

                    vect_complex_s16_conj_mul_quake(vR, vD, vC, 8, 0, 0);

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

