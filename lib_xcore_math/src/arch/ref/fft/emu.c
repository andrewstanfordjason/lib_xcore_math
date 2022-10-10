#include <stdint.h>
#include <assert.h>
#include <string.h>

void add_sub(int16_t * output_reg, const int16_t * input_reg, int shift_mode, int invert){

    assert(shift_mode <= 1);
    assert(shift_mode >= -1);

    for(int i = 0; i < 8; i++){
        int32_t sum = (int32_t)input_reg[i+8] + (int32_t)input_reg[i];
        int32_t dif = (int32_t)input_reg[i+8] - (int32_t)input_reg[i];

        if (invert)
            dif = -dif;

        if (shift_mode>=0){
            sum >>= shift_mode;
            dif >>= shift_mode;
        } else {
            sum <<= (-shift_mode);;
            dif <<= (-shift_mode);
        }
        if (sum > INT16_MAX) sum = INT16_MAX;
        if (sum < INT16_MIN) sum = INT16_MIN;
        if (dif > INT16_MAX) dif = INT16_MAX;
        if (dif < INT16_MIN) dif = INT16_MIN;
        
        output_reg[i] = sum;
        output_reg[i+8] = dif;
    }
}
void VADSB(int8_t * output_reg, const int8_t * input_reg, int shift_mode){
    add_sub((int16_t * )output_reg, (const int16_t * )input_reg, shift_mode, 0);
}

void VSBAD(int8_t * output_reg, const int8_t * input_reg, int shift_mode){
    add_sub((int16_t * )output_reg, (const int16_t * )input_reg, shift_mode, 1);
}
#include <stdio.h>
void complex_mul_32(int8_t * vR_out, 
    const int8_t * vR_in, 
    const int8_t * vD_in, 
    const int8_t * vC_in, int mode, int swap_C, int subword, int neg) {

    assert (subword >= 0);
    assert (subword <= 3);
    assert (swap_C >= 0);
    assert (swap_C <= 1);
    assert (neg >= 0);
    assert (neg <= 1);

    for (int i=0; i<4; i++){
        const int32_t * vR_in_32 = (const int32_t* )vR_in;
        const int32_t * vD_32 = (const int32_t* )vD_in;
        const int32_t * vC_32 = (const int32_t* )vC_in;
        int32_t * vR_out_32 = (int32_t* )vR_out;

        int re_idx = 2*i + swap_C;
        int im_idx = 2*i + 1 - swap_C;

        uint8_t vC_re_subword = (uint8_t)(vC_32[re_idx]>>(8*subword));
        uint8_t vC_im_subword = (uint8_t)(vC_32[im_idx]>>(8*subword));

        int64_t a = (int64_t)vC_re_subword * (int64_t)vD_32[2*i];
        int64_t b = (int64_t)vC_im_subword * (int64_t)vD_32[2*i + 1];

        int64_t carry = (int64_t)(vR_in_32[re_idx]);
        int64_t sum = (carry<<1) + (int64_t)a ;

        if (neg){
            sum -= (int64_t)b;
        } else {
            sum += (int64_t)b;
        }
        
        int64_t round = (sum + (1<<8)) >> 9;

        if (round > INT32_MAX) round = INT32_MAX;
        if (round < INT32_MIN) round = INT32_MIN;

        vR_out_32[re_idx] = round;
    }
}

void complex_mul_32_final(int8_t * vR_out, 
    const int8_t * vR_in, 
    const int8_t * vD_in, 
    const int8_t * vC_in, int mode, int swap_C, int neg) {

    for (int i=0; i<4; i++){
        const int32_t * vR_in_32 = (const int32_t* )vR_in;
        const int32_t * vD_32 = (const int32_t* )vD_in;
        const int32_t * vC_32 = (const int32_t* )vC_in;
        int32_t * vR_out_32 = (int32_t* )vR_out;

        int re_idx = 2*i + swap_C;
        int im_idx = 2*i + 1 - swap_C;

        int8_t vC_re_subword = (int8_t)(vC_32[re_idx]>>24);
        int8_t vC_im_subword = (int8_t)(vC_32[im_idx]>>24);

        int64_t a = (int64_t)vC_re_subword * (int64_t)vD_32[2*i];
        int64_t b = (int64_t)vC_im_subword * (int64_t)vD_32[2*i + 1];

        int64_t carry = (int64_t)vR_in_32[re_idx];
        int64_t sum = (carry<<1) + (int64_t)a ;
        if (neg){
            sum -= (int64_t)b;
        } else {
            sum += (int64_t)b;
        }
        
        int64_t round = (sum + (1<<5)) >> 6;

        if (round > INT32_MAX) round = INT32_MAX;
        if (round < INT32_MIN) round = INT32_MIN;

        vR_out_32[re_idx] = round;
    }
}

void complex_mul_16(int8_t * vR_out, 
    const int8_t * vR_in, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode, int swap_C, int subword, int neg) 
{
    for (int i=0;i<8;i++){

        const int16_t * vR_in_32 = (const int16_t* )vR_in;
        const int16_t * vD_32 = (const int16_t* )vD_in;
        const int16_t * vC_32 = (const int16_t* )vC_in;
        int16_t * vR_out_32 = (int16_t* )vR_out;

        int re_idx = 8*i + swap_C;
        int im_idx = 8*i + 1 - swap_C;

        uint8_t vC_re_subword = (uint8_t)(vC_32[re_idx]>>(8*subword));
        uint8_t vC_im_subword = (uint8_t)(vC_32[im_idx]>>(8*subword));

        int64_t a = (int64_t)vC_re_subword * (int64_t)vD_32[re_idx];
        int64_t b = (int64_t)vC_im_subword * (int64_t)vD_32[im_idx];

        int64_t carry = (int64_t)vR_in_32[re_idx];
        int64_t sum = (carry<<1) + (int64_t)a ;
        if (neg){
            sum -= (int64_t)b;
        } else {
            sum += (int64_t)b;
        }
        
        int64_t round = (sum + (1<<8)) >> 9;

        if (round > INT16_MAX) round = INT16_MAX;
        if (round < INT16_MIN) round = INT16_MIN;

        vR_out_32[re_idx] = round;
    }
}

void complex_mul_16_final(int8_t * vR_out, 
    const int8_t * vR_in, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode, int swap_C, int neg) 
{
    for (int i=0;i<8;i++){

        const int16_t * vR_in_32 = (const int16_t* )vR_in;
        const int16_t * vD_32 = (const int16_t* )vD_in;
        const int16_t * vC_32 = (const int16_t* )vC_in;
        int16_t * vR_out_32 = (int16_t* )vR_out;

        int re_idx = 8*i + swap_C;
        int im_idx = 8*i + 1 - swap_C;

        int8_t vC_re_subword = (int8_t)(vC_32[re_idx]>>8);
        int8_t vC_im_subword = (int8_t)(vC_32[im_idx]>>8);

        int64_t a = (int64_t)vC_re_subword * (int64_t)vD_32[re_idx];
        int64_t b = (int64_t)vC_im_subword * (int64_t)vD_32[im_idx];

        int64_t carry = (int64_t)vR_in_32[re_idx];
        int64_t sum = (carry<<1) + (int64_t)a ;
        if (neg){
            sum -= (int64_t)b;
        } else {
            sum += (int64_t)b;
        }
        
        int64_t round = (sum + (1<<8)) >> 9;

        if (round > INT16_MAX) round = INT16_MAX;
        if (round < INT16_MIN) round = INT16_MIN;

        vR_out_32[re_idx] = round;
    }
}

/////////////////////////////////////////////////////////////////////////

void VCMR_0(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode) {
    int8_t zero[32];
    memset(zero, 0, 32);
    if (mode == 16){
        // complex_mul_16(vR_out, zero, vD_in, vC_in, mode, 0, 0, 1);
    } else if (mode == 32){
        complex_mul_32(vR_out, zero, vD_in, vC_in, mode, 0, 0, 1);
    } else {
        assert(0);
    }
}

void VCMR_1(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in, 
    const int8_t * vR_in, int mode) {
    if (mode == 16){
        // complex_mul_16_final(vR_out, vR_in, vC_in, vD_in, mode, 0, 1, 1);
    } else if (mode == 32){
        complex_mul_32(vR_out, vR_in, vD_in, vC_in, mode, 0, 1, 1);
    } else {
        assert(0);
    }
}

void VCMR_2(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        complex_mul_32(vR_out, vR_in, vD_in, vC_in, mode, 0, 2, 1);
    } else {
        assert(0);
    }
}

void VCMR_3(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        complex_mul_32_final(vR_out, vR_in, vD_in, vC_in, mode, 0, 1);
    } else {
        assert(0);
    }
}

/////////////////////////////////////////////////////////////////////////


void VCMI_0(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode) {
    int8_t zero[32];
    memset(zero, 0, 32);
    if (mode == 16){
        // complex_mul_16(vR_out, zero, vC_in, vD_in, mode, 0, 0);
    } else if (mode == 32){
        complex_mul_32(vR_out, zero, vD_in, vC_in, mode, 1, 0, 0);
    } else {
        assert(0);
    }
}

void VCMI_1(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in, 
    const int8_t * vR_in, int mode) {
    if (mode == 16){
        // complex_mul_16_final(vR_out, vR_in, vC_in, vD_in, mode, 1, 0);
    } else if (mode == 32){
        complex_mul_32(vR_out, vR_in, vD_in, vC_in, mode, 1, 1, 0);
    } else {
        assert(0);
    }
}

void VCMI_2(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        complex_mul_32(vR_out, vR_in, vD_in, vC_in, mode, 1, 2, 0);
    } else {
        assert(0);
    }
}

void VCMI_3(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        complex_mul_32_final(vR_out, vR_in, vD_in, vC_in, mode, 1, 0);
    } else {
        assert(0);
    }
}


/////////////////////////////////////////////////////////////////////////

void VCMCR_0(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode) {
    int8_t zero[32];
    memset(zero, 0, 32);
    if (mode == 16){
        // complex_mul_16(vR_out, zero, vC_in, vD_in, mode, 0, 0);
    } else if (mode == 32){
        complex_mul_32(vR_out, zero, vD_in, vC_in, mode, 0, 0, 0);
    } else {
        assert(0);
    }
}

void VCMCR_1(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 16){
        // complex_mul_16_final(vR_out, vR_in, vC_in, vD_in, mode, 1, 0);
    } else if (mode == 32){
        complex_mul_32(vR_out, vR_in, vD_in, vC_in, mode, 0, 1, 0);
    } else {
        assert(0);
    }
}

void VCMCR_2(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        complex_mul_32(vR_out, vR_in, vD_in, vC_in, mode, 0, 2, 0);
    } else {
        assert(0);
    }
}

void VCMCR_3(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        complex_mul_32_final(vR_out, vR_in, vD_in, vC_in, mode, 0, 0);
    } else {
        assert(0);
    }
}

/////////////////////////////////////////////////////////////////////////


void VCMCI_0(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode) {
    int8_t zero[32];
    memset(zero, 0, 32);
    if (mode == 16){
        // complex_mul_16(vR_out, zero, vC_in, vD_in, mode, 0, 0);
    } else if (mode == 32){
        complex_mul_32(vR_out, zero, vD_in, vC_in, mode, 1, 0, 1);
    } else {
        assert(0);
    }
}

void VCMCI_1(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 16){
        // complex_mul_16_final(vR_out, vR_in, vC_in, vD_in, mode, 1, 0);
    } else if (mode == 32){
        complex_mul_32(vR_out, vR_in, vD_in, vC_in, mode, 1, 1, 1);
    } else {
        assert(0);
    }
}

void VCMCI_2(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        complex_mul_32(vR_out, vR_in, vD_in, vC_in, mode, 1, 2, 1);
    } else {
        assert(0);
    }
}

void VCMCI_3(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        complex_mul_32_final(vR_out, vR_in, vD_in, vC_in, mode, 1, 1);
    } else {
        assert(0);
    }
}


