// Copyright 2020-2022 XMOS LIMITED.
// This Software is subject to the terms of the XMOS Public Licence: Version 1.

#include <stdint.h>
#include <stdio.h>

#include "xmath/xmath.h"
#include "../../../vect/vpu_helper.h"
#include "xmath/xs3/vpu_scalar_ops.h"


////////////////////////////////////////
//      16-Bit Multiplication         //
////////////////////////////////////////

// complex vector multiplied by real vector

headroom_t vect_complex_s16_real_mul(
    int16_t a_real[],
    int16_t a_imag[],
    const int16_t b_real[],
    const int16_t b_imag[],
    const int16_t c[],
    const unsigned length,
    const right_shift_t sat)
{
    for(int k = 0; k < length; k++){
        
        complex_s32_t B = {b_real[k], b_imag[k]};
        int32_t C = c[k];

        struct{ int64_t re; int64_t im; } P = {
            ((int64_t)B.re * C),
            ((int64_t)B.im * C),
        };

        P.re = ROUND_SHR(P.re, sat);
        P.im = ROUND_SHR(P.im, sat);

        a_real[k] = SAT(16)(P.re);
        a_imag[k] = SAT(16)(P.im);
    }

    return vect_complex_s16_headroom(a_real, a_imag, length);
}


// complex vector multiplied by complex vector

headroom_t vect_complex_s16_mul(
    int16_t a_real[],
    int16_t a_imag[],
    const int16_t b_real[],
    const int16_t b_imag[],
    const int16_t c_real[],
    const int16_t c_imag[],
    const unsigned length,
    const right_shift_t sat)
{
    for(int k = 0; k < length; k++){
        
        complex_s32_t B = {b_real[k], b_imag[k]};
        complex_s32_t C = {c_real[k], c_imag[k]};

        struct{ int64_t re; int64_t im; } P = {
            ((int64_t)B.re * C.re) - B.im * C.im,
            ((int64_t)B.re * C.im) + B.im * C.re,
        };

        P.re = ROUND_SHR(P.re, sat);
        P.im = ROUND_SHR(P.im, sat);

        a_real[k] = SAT(16)(P.re);
        a_imag[k] = SAT(16)(P.im);
    }

    return vect_complex_s16_headroom(a_real, a_imag, length);
}


// complex vector (conjugate) multiplied by complex vector

headroom_t vect_complex_s16_conj_mul(
    int16_t a_real[],
    int16_t a_imag[],
    const int16_t b_real[],
    const int16_t b_imag[],
    const int16_t c_real[],
    const int16_t c_imag[],
    const unsigned length,
    const right_shift_t sat)
{
    for(int k = 0; k < length; k++){
        
        complex_s32_t B = {b_real[k], b_imag[k]};
        complex_s32_t C = {c_real[k], c_imag[k]};

        struct{ int64_t re; int64_t im; } P = {
            ((int64_t)B.re * C.re) + B.im * C.im,
            ((int64_t)B.im * C.re) - B.re * C.im,
        };

        P.re = ROUND_SHR(P.re, sat);
        P.im = ROUND_SHR(P.im, sat);

        a_real[k] = SAT(16)(P.re);
        a_imag[k] = SAT(16)(P.im);
    }

    return vect_complex_s16_headroom(a_real, a_imag, length);
}


// complex vector multiplied by complex vector

headroom_t vect_complex_s16_scale(
    int16_t a_real[],
    int16_t a_imag[],
    const int16_t b_real[],
    const int16_t b_imag[],
    const int16_t c_real,
    const int16_t c_imag,
    const unsigned length,
    const right_shift_t sat)
{

    for(int k = 0; k < length; k++){
        
        complex_s32_t B = {b_real[k], b_imag[k]};

        struct{ int64_t re; int64_t im; } P = {
            ((int64_t)B.re * c_real) - B.im * c_imag,
            ((int64_t)B.re * c_imag) + B.im * c_real,
        };

        P.re = ROUND_SHR(P.re, sat);
        P.im = ROUND_SHR(P.im, sat);

        a_real[k] = SAT(16)(P.re);
        a_imag[k] = SAT(16)(P.im);
    }

    return vect_complex_s16_headroom(a_real, a_imag, length);
}


////////////////////////////////////////
//      16-Bit Multiplication  quake        //
////////////////////////////////////////


// complex vector multiplied by complex vector

headroom_t vect_complex_s16_mul_quake_original(
    complex_s16_t a[],
    const complex_s16_t b[],
    const complex_s16_t c[],
    const unsigned length,
    const right_shift_t b_shr,
    const right_shift_t c_shr)
{
    for(int k = 0; k < length; k++){
        
        complex_s16_t B = {
            ASHR(16)(b[k].re, b_shr), 
            ASHR(16)(b[k].im, b_shr),
        };
        complex_s16_t C = {
            ASHR(16)(c[k].re, c_shr), 
            ASHR(16)(c[k].im, c_shr),
        };

        int32_t q1 = ROUND_SHR( ((int32_t)B.re) * C.re, 14 );
        int32_t q2 = ROUND_SHR( ((int32_t)B.im) * C.im, 14 );
        int32_t q3 = ROUND_SHR( ((int32_t)B.re) * C.im, 14 );
        int32_t q4 = ROUND_SHR( ((int32_t)B.im) * C.re, 14 );

        a[k].re = SAT(16)(q1 - q2);
        a[k].im = SAT(16)(q3 + q4);
    }

    return vect_complex_s16_headroom_quake(a, length);
}

headroom_t vect_complex_s16_mul_quake(
    complex_s16_t a[],
    const complex_s16_t b[],
    const complex_s16_t c[],
    const unsigned length,
    const right_shift_t b_shr,
    const right_shift_t c_shr)
{
    for(int k = 0; k < length; k++){

        complex_s16_t B = {
            ASHR(16)(b[k].re, b_shr), 
            ASHR(16)(b[k].im, b_shr),
        };
        complex_s16_t C = {
            ASHR(16)(c[k].re, c_shr), 
            ASHR(16)(c[k].im, c_shr),
        };
        
        uint8_t b_re_l = (uint8_t)B.re;
        uint8_t b_im_l = (uint8_t)B.im;

        int32_t q1_0 = (int32_t)b_re_l * C.re;
        int32_t q2_0 = (int32_t)b_im_l * C.im;

        int32_t q3_0 = (int32_t)b_re_l * C.im;
        int32_t q4_0 = (int32_t)b_im_l * C.re;

        int32_t t_re = q1_0 - q2_0;
        int32_t t_im = q3_0 + q4_0;

        int16_t inter_re = ROUND_SHR(t_re, 9);
        int16_t inter_im = ROUND_SHR(t_im, 9);
        
        //now the second half

        int8_t b_re_h = (int8_t)(B.re>>8);
        int8_t b_im_h = (int8_t)(B.im>>8);

        int32_t q1_1 = (int32_t)b_re_h * C.re;
        int32_t q2_1 = (int32_t)b_im_h * C.im;

        int32_t q3_1 = (int32_t)b_re_h * C.im;
        int32_t q4_1 = (int32_t)b_im_h * C.re;

        int32_t a_re = ROUND_SHR(q1_1 - q2_1 + ((int32_t)inter_re << 1), 6);
        int32_t a_im = ROUND_SHR(q3_1 + q4_1 + ((int32_t)inter_im << 1), 6);

        a[k].re = SAT(16)(a_re);
        a[k].im = SAT(16)(a_im);
    }

    return vect_complex_s16_headroom_quake(a, length);
}


// complex vector (conjugate) multiplied by complex vector

headroom_t vect_complex_s16_conj_mul_quake_original(
    complex_s16_t a[],
    const complex_s16_t b[],
    const complex_s16_t c[],
    const unsigned length,
    const right_shift_t b_shr,
    const right_shift_t c_shr)
{
    for(int k = 0; k < length; k++){
        
        complex_s16_t B = {
            ASHR(16)(b[k].re, b_shr), 
            ASHR(16)(b[k].im, b_shr),
        };

        complex_s16_t C = {
            ASHR(16)(c[k].re, c_shr), 
            ASHR(16)(c[k].im, c_shr),
        };

        int32_t q1 = ROUND_SHR( ((int32_t)B.re) * C.re, 14 );
        int32_t q2 = ROUND_SHR( ((int32_t)B.im) * C.im, 14 );
        int32_t q3 = ROUND_SHR( ((int32_t)B.re) * C.im, 14 );
        int32_t q4 = ROUND_SHR( ((int32_t)B.im) * C.re, 14 );

        a[k].re = SAT(16)(q1 + q2);
        a[k].im = SAT(16)(q4 - q3);
    }

    return vect_complex_s16_headroom_quake(a, length);
}

headroom_t vect_complex_s16_conj_mul_quake(
    complex_s16_t a[],
    const complex_s16_t b[],
    const complex_s16_t c[],
    const unsigned length,
    const right_shift_t b_shr,
    const right_shift_t c_shr)
{
    for(int k = 0; k < length; k++){
        
        complex_s16_t B = {
            ASHR(16)(b[k].re, b_shr), 
            ASHR(16)(b[k].im, b_shr),
        };

        complex_s16_t C = {
            ASHR(16)(c[k].re, c_shr), 
            ASHR(16)(c[k].im, c_shr),
        };
        uint8_t b_re_l = (uint8_t)B.re;
        uint8_t b_im_l = (uint8_t)B.im;

        int32_t q1_0 = (int32_t)b_re_l * C.re;
        int32_t q2_0 = (int32_t)b_im_l * C.im;

        int32_t q3_0 = (int32_t)b_re_l * C.im;
        int32_t q4_0 = (int32_t)b_im_l * C.re;

        int32_t t_re = q1_0 + q2_0;
        int32_t t_im = q4_0 - q3_0;

        int16_t inter_re = ROUND_SHR(t_re, 9);
        int16_t inter_im = ROUND_SHR(t_im, 9);
        
        //now the second half

        int8_t b_re_h = (int8_t)(B.re>>8);
        int8_t b_im_h = (int8_t)(B.im>>8);

        int32_t q1_1 = (int32_t)b_re_h * C.re;
        int32_t q2_1 = (int32_t)b_im_h * C.im;

        int32_t q3_1 = (int32_t)b_re_h * C.im;
        int32_t q4_1 = (int32_t)b_im_h * C.re;

        int32_t a_re = ROUND_SHR(q1_1 + q2_1 + ((int32_t)inter_re << 1), 6);
        int32_t a_im = ROUND_SHR(q4_1 - q3_1 + ((int32_t)inter_im << 1), 6);

        a[k].re = SAT(16)(a_re);
        a[k].im = SAT(16)(a_im);
    }

    return vect_complex_s16_headroom_quake(a, length);
}





////////////////////////////////////////
//      32-Bit Multiplication         //
////////////////////////////////////////


// complex vector multiplied by real vector

headroom_t vect_complex_s32_real_mul(
    complex_s32_t a[],
    const complex_s32_t b[],
    const int32_t c[],
    const unsigned length,
    const right_shift_t b_shr,
    const right_shift_t c_shr)
{
    for(int k = 0; k < length; k++){
        
        complex_s32_t B = {
            ASHR(32)(b[k].re, b_shr), 
            ASHR(32)(b[k].im, b_shr),
        };

        int32_t C = ASHR(32)(c[k], c_shr);

        a[k].re = SAT(32)(ROUND_SHR(((int64_t)B.re) * C, 30));
        a[k].im = SAT(32)(ROUND_SHR(((int64_t)B.im) * C, 30));
    }


    return vect_complex_s32_headroom(a, length);
}


// complex vector multiplied by complex vector

headroom_t vect_complex_s32_mul(
    complex_s32_t a[],
    const complex_s32_t b[],
    const complex_s32_t c[],
    const unsigned length,
    const right_shift_t b_shr,
    const right_shift_t c_shr)
{
    for(int k = 0; k < length; k++){
        
        complex_s32_t B = {
            ASHR(32)(b[k].re, b_shr), 
            ASHR(32)(b[k].im, b_shr),
        };
        complex_s32_t C = {
            ASHR(32)(c[k].re, c_shr), 
            ASHR(32)(c[k].im, c_shr),
        };

        int64_t q1 = ROUND_SHR( ((int64_t)B.re) * C.re, 30 );
        int64_t q2 = ROUND_SHR( ((int64_t)B.im) * C.im, 30 );
        int64_t q3 = ROUND_SHR( ((int64_t)B.re) * C.im, 30 );
        int64_t q4 = ROUND_SHR( ((int64_t)B.im) * C.re, 30 );

        a[k].re = SAT(32)(q1 - q2);
        a[k].im = SAT(32)(q3 + q4);
    }

    return vect_complex_s32_headroom(a, length);
}


// complex vector (conjugate) multiplied by complex vector

headroom_t vect_complex_s32_conj_mul(
    complex_s32_t a[],
    const complex_s32_t b[],
    const complex_s32_t c[],
    const unsigned length,
    const right_shift_t b_shr,
    const right_shift_t c_shr)
{
    for(int k = 0; k < length; k++){
        
        complex_s32_t B = {
            ASHR(32)(b[k].re, b_shr), 
            ASHR(32)(b[k].im, b_shr),
        };

        complex_s32_t C = {
            ASHR(32)(c[k].re, c_shr), 
            ASHR(32)(c[k].im, c_shr),
        };

        int64_t q1 = ROUND_SHR( ((int64_t)B.re) * C.re, 30 );
        int64_t q2 = ROUND_SHR( ((int64_t)B.im) * C.im, 30 );
        int64_t q3 = ROUND_SHR( ((int64_t)B.re) * C.im, 30 );
        int64_t q4 = ROUND_SHR( ((int64_t)B.im) * C.re, 30 );

        a[k].re = SAT(32)(q1 + q2);
        a[k].im = SAT(32)(q4 - q3);
    }

    return vect_complex_s32_headroom(a, length);
}



// complex vector multiplied by complex scalar

headroom_t vect_complex_s32_scale(
    complex_s32_t a[],
    const complex_s32_t b[],
    const int32_t c_real,
    const int32_t c_imag,
    const unsigned length,
    const right_shift_t b_shr,
    const right_shift_t c_shr)
{

    const complex_s32_t C = {
        vlashr32(c_real, c_shr),
        vlashr32(c_imag, c_shr),
    };

    for(int k = 0; k < length; k++){
        
        const complex_s32_t B = {
            vlashr32(b[k].re, b_shr),
            vlashr32(b[k].im, b_shr),
        };

        a[k].re = vcmr32( B, C );
        a[k].im = vcmi32( B, C );

    }

    return vect_complex_s32_headroom(a, length);
}



