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

void VCMR_0(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode) {

    if (mode == 16){
        for (int i=0;i<8;i++){
            int64_t a = (int64_t)(uint8_t)vC_in[4*i  ] * (int64_t)((int16_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[4*i+2] * (int64_t)((int16_t* )vD_in)[2*i+1];
            ((int16_t* )vR_out)[2*i] = (a - b + (1<<8)) >> 9;
        }
    } else if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i  ] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i+4] * (int64_t)((int32_t* )vD_in)[2*i+1];
            ((int32_t* )vR_out)[2*i] = (a - b + (1<<8)) >> 9;
        }
    } else {
        assert(0);
    }
}

void VCMR_1(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in, 
    const int8_t * vR_in, int mode) {
    if (mode == 16){
        for (int i=0;i<8;i++){
            int64_t a = (int64_t)(int8_t)vC_in[4*i+1] * (int64_t)((int16_t* )vD_in)[2*i];
            int64_t b = (int64_t)(int8_t)vC_in[4*i+3] * (int64_t)((int16_t* )vD_in)[2*i + 1];
            
            int64_t carry = (int64_t)((int16_t*)vR_in)[2*i];
            int64_t round = (a - b + (carry<<1) + (1<<5)) >> 6;

            if (round > INT16_MAX) round = INT16_MAX;
            if (round < INT16_MIN) round = INT16_MIN;

            ((int16_t* )vR_out)[2*i] = (int16_t)round;

        }
    } else if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i+1] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i+5] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i];
            ((int32_t* )vR_out)[2*i] = (a - b + (carry<<1) + (1<<8)) >> 9;
            //FIXME - write a test for overflow due to rounding
        }
    } else {
        assert(0);
    }
}

void VCMR_2(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i+2] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i+6] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i];
            ((int32_t* )vR_out)[2*i] = (a - b + (carry<<1) + (1<<8)) >> 9;
            //FIXME - write a test for overflow due to rounding
        }
    } else {
        assert(0);
    }
}

void VCMR_3(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(int8_t)vC_in[8*i+3] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(int8_t)vC_in[8*i+7] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i];
            int64_t round = (a - b + (carry<<1) + (1<<5)) >> 6;

            if (round > INT32_MAX) round = INT32_MAX;
            if (round < INT32_MIN) round = INT32_MIN;

            ((int32_t* )vR_out)[2*i] = round;

        }
    } else {
        assert(0);
    }
}

/////////////////////////////////////////////////////////////////////////

void VCMI_0(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode) {

    if (mode == 16){
        for (int i=0;i<8;i++){
            int64_t a = (int64_t)(uint8_t)vC_in[4*i+2] * (int64_t)((int16_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[4*i  ] * (int64_t)((int16_t* )vD_in)[2*i+1];
            ((int16_t* )vR_out)[2*i+1] = (a + b + (1<<8)) >> 9;
        }
    } else if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i+4] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i  ] * (int64_t)((int32_t* )vD_in)[2*i+1];
            ((int32_t* )vR_out)[2*i+1] = (a + b + (1<<8)) >> 9;
        }
    } else {
        assert(0);
    }
}

void VCMI_1(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in, 
    const int8_t * vR_in, int mode) {
    if (mode == 16){
        for (int i=0;i<8;i++){
            int64_t a = (int64_t)(int8_t)vC_in[4*i+3] * (int64_t)((int16_t* )vD_in)[2*i];
            int64_t b = (int64_t)(int8_t)vC_in[4*i+1] * (int64_t)((int16_t* )vD_in)[2*i + 1];
            
            int64_t carry = (int64_t)((int16_t*)vR_in)[2*i+1];
            int64_t round = (a + b + (carry<<1) + (1<<5)) >> 6;

            if (round > INT16_MAX) round = INT16_MAX;
            if (round < INT16_MIN) round = INT16_MIN;

            ((int16_t* )vR_out)[2*i+1] = (int16_t)round;

        }
    } else if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i+5] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i+1] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i+1];
            ((int32_t* )vR_out)[2*i+1] = (a + b + (carry<<1) + (1<<8)) >> 9;
            //FIXME - write a test for overflow due to rounding
        }
    } else {
        assert(0);
    }
}

void VCMI_2(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i+6] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i+2] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i+1];
            ((int32_t* )vR_out)[2*i+1] = (a + b + (carry<<1) + (1<<8)) >> 9;
            //FIXME - write a test for overflow due to rounding
        }
    } else {
        assert(0);
    }
}

void VCMI_3(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(int8_t)vC_in[8*i+7] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(int8_t)vC_in[8*i+3] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i+1];
            int64_t round = (a + b + (carry<<1) + (1<<5)) >> 6;

            if (round > INT32_MAX) round = INT32_MAX;
            if (round < INT32_MIN) round = INT32_MIN;

            ((int32_t* )vR_out)[2*i+1] = round;

        }
    } else {
        assert(0);
    }
}

/////////////////////////////////////////////////////////////////////////

void VCMCR_0(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode) {

    if (mode == 16){
        for (int i=0;i<8;i++){
            int64_t a = (int64_t)(uint8_t)vC_in[4*i  ] * (int64_t)((int16_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[4*i+2] * (int64_t)((int16_t* )vD_in)[2*i+1];
            ((int16_t* )vR_out)[2*i] = (a + b + (1<<8)) >> 9;
        }
    } else if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i  ] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i+4] * (int64_t)((int32_t* )vD_in)[2*i+1];
            ((int32_t* )vR_out)[2*i] = (a + b + (1<<8)) >> 9;
        }
    } else {
        assert(0);
    }
}

void VCMCR_1(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in, 
    const int8_t * vR_in, int mode) {
    if (mode == 16){
        for (int i=0;i<8;i++){
            int64_t a = (int64_t)(int8_t)vC_in[4*i+1] * (int64_t)((int16_t* )vD_in)[2*i];
            int64_t b = (int64_t)(int8_t)vC_in[4*i+3] * (int64_t)((int16_t* )vD_in)[2*i + 1];
            
            int64_t carry = (int64_t)((int16_t*)vR_in)[2*i];
            int64_t round = (a + b + (carry<<1) + (1<<5)) >> 6;

            if (round > INT16_MAX) round = INT16_MAX;
            if (round < INT16_MIN) round = INT16_MIN;

            ((int16_t* )vR_out)[2*i] = (int16_t)round;

        }
    } else if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i+1] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i+5] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i];
            ((int32_t* )vR_out)[2*i] = (a + b + (carry<<1) + (1<<8)) >> 9;
            //FIXME - write a test for overflow due to rounding
        }
    } else {
        assert(0);
    }
}

void VCMCR_2(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i+2] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i+6] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i];
            ((int32_t* )vR_out)[2*i] = (a + b + (carry<<1) + (1<<8)) >> 9;
            //FIXME - write a test for overflow due to rounding
        }
    } else {
        assert(0);
    }
}

void VCMCR_3(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(int8_t)vC_in[8*i+3] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(int8_t)vC_in[8*i+7] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i];
            int64_t round = (a + b + (carry<<1) + (1<<5)) >> 6;

            if (round > INT32_MAX) round = INT32_MAX;
            if (round < INT32_MIN) round = INT32_MIN;

            ((int32_t* )vR_out)[2*i] = round;

        }
    } else {
        assert(0);
    }
}

/////////////////////////////////////////////////////////////////////////

void VCMCI_0(int8_t * vR_out, 
    const int8_t * vC_in, 
    const int8_t * vD_in, int mode) {

    if (mode == 16){
        for (int i=0;i<8;i++){
            int64_t a = (int64_t)(uint8_t)vC_in[4*i+2] * (int64_t)((int16_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[4*i  ] * (int64_t)((int16_t* )vD_in)[2*i+1];
            ((int16_t* )vR_out)[2*i+1] = (a - b + (1<<8)) >> 9;
        }
    } else if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i+4] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i  ] * (int64_t)((int32_t* )vD_in)[2*i+1];
            ((int32_t* )vR_out)[2*i+1] = (a - b + (1<<8)) >> 9;
        }
    } else {
        assert(0);
    }
}

void VCMCI_1(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in, 
    const int8_t * vR_in, int mode) {
    if (mode == 16){
        for (int i=0;i<8;i++){
            int64_t a = (int64_t)(int8_t)vC_in[4*i+3] * (int64_t)((int16_t* )vD_in)[2*i];
            int64_t b = (int64_t)(int8_t)vC_in[4*i+1] * (int64_t)((int16_t* )vD_in)[2*i + 1];
            
            int64_t carry = (int64_t)((int16_t*)vR_in)[2*i+1];
            int64_t round = (a - b + (carry<<1) + (1<<5)) >> 6;

            if (round > INT16_MAX) round = INT16_MAX;
            if (round < INT16_MIN) round = INT16_MIN;

            ((int16_t* )vR_out)[2*i+1] = (int16_t)round;

        }
    } else if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i+5] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i+1] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i+1];
            ((int32_t* )vR_out)[2*i+1] = (a - b + (carry<<1) + (1<<8)) >> 9;
            //FIXME - write a test for overflow due to rounding
        }
    } else {
        assert(0);
    }
}

void VCMCI_2(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(uint8_t)vC_in[8*i+6] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(uint8_t)vC_in[8*i+2] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i+1];
            ((int32_t* )vR_out)[2*i+1] = (a - b + (carry<<1) + (1<<8)) >> 9;
            //FIXME - write a test for overflow due to rounding
        }
    } else {
        assert(0);
    }
}

void VCMCI_3(int8_t * vR_out,
    const int8_t * vC_in, 
    const int8_t * vD_in,
    const int8_t * vR_in, int mode) {
    if (mode == 32){
        for (int i=0; i<4; i++){
            int64_t a = (int64_t)(int8_t)vC_in[8*i+7] * (int64_t)((int32_t* )vD_in)[2*i  ];
            int64_t b = (int64_t)(int8_t)vC_in[8*i+3] * (int64_t)((int32_t* )vD_in)[2*i+1];
            int64_t carry = (int64_t)((int32_t*)vR_in)[2*i+1];
            int64_t round = (a - b + (carry<<1) + (1<<5)) >> 6;

            if (round > INT32_MAX) round = INT32_MAX;
            if (round < INT32_MIN) round = INT32_MIN;

            ((int32_t* )vR_out)[2*i+1] = round;

        }
    } else {
        assert(0);
    }
}


void VLADSB(int8_t * vD_out, int8_t * vR_out, const int8_t * Mem_in, const int8_t * vR_in, int mode,  int shift_mode){
    if (mode == 16){
        const int16_t * vR_in_16 = (const int16_t* )vR_in;
        const int16_t * Mem_16 = (const int16_t* )Mem_in;
        int16_t * vD_out_16 = (int16_t* )vD_out;
        int16_t * vR_out_16 = (int16_t* )vR_out;

        for(int i = 0; i < 16; i++){
            int32_t sum = (int32_t)Mem_16[i] + (int32_t)vR_in_16[i];
            int32_t dif = (int32_t)Mem_16[i] - (int32_t)vR_in_16[i];

            if (shift_mode>=0){
                sum >>= shift_mode;
                dif >>= shift_mode;
            } else {
                sum <<= (-shift_mode);
                dif <<= (-shift_mode);
            }
            if (sum > INT16_MAX) sum = INT16_MAX;
            if (sum < INT16_MIN) sum = INT16_MIN;
            if (dif > INT16_MAX) dif = INT16_MAX;
            if (dif < INT16_MIN) dif = INT16_MIN;
            
            vR_out_16[i] = sum;
            vD_out_16[i] = dif;
        }
    } else if (mode == 32){
        const int32_t * vR_in_32 = (const int32_t* )vR_in;
        const int32_t * Mem_32 = (const int32_t* )Mem_in;
        int32_t * vD_out_32 = (int32_t* )vD_out;
        int32_t * vR_out_32 = (int32_t* )vR_out;

        for(int i = 0; i < 8; i++){
            int64_t sum = (int64_t)Mem_32[i] + (int64_t)vR_in_32[i];
            int64_t dif = (int64_t)Mem_32[i] - (int64_t)vR_in_32[i];

            if (shift_mode>=0){
                sum >>= shift_mode;
                dif >>= shift_mode;
            } else {
                sum <<= (-shift_mode);
                dif <<= (-shift_mode);
            }
            if (sum > INT32_MAX) sum = INT32_MAX;
            if (sum < INT32_MIN) sum = INT32_MIN;
            if (dif > INT32_MAX) dif = INT32_MAX;
            if (dif < INT32_MIN) dif = INT32_MIN;
            
            vR_out_32[i] = sum;
            vD_out_32[i] = dif;
        }
    } else {
        assert(0);
    }
}

#include "xmath/xmath.h"
#define SAT16(VAL)      (((VAL) >= INT16_MAX)? INT16_MAX : (((VAL) <= INT16_MIN)? INT16_MIN : (VAL)))
#define SAT32(VAL)      (((VAL) >= INT32_MAX)? INT32_MAX : (((VAL) <= INT32_MIN)? INT32_MIN : (VAL)))

#define SHR(VAL, SHR)   ((VAL)>>SHR)
#define ASHR16(VAL, SHR_BITS)   (SAT16(((SHR_BITS) >= 0)? SHR((VAL),(SHR_BITS)) : (SAT16( ((int32_t)(VAL))<<(-(SHR_BITS)) ) )))
#define ASHR32(VAL, SHR_BITS)   (SAT32(((SHR_BITS) >= 0)? SHR((VAL),(SHR_BITS)) : (SAT32( ((int64_t)(VAL))<<(-(SHR_BITS)) ) )))
#define ASHR(BITS)  ASHR##BITS

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

void VFTFF(int8_t * vD, int mode, int shift_mode){
    if (mode == 16){
        vftff_quake_s16((complex_s16_t*)vD, shift_mode);
    } else if (mode == 32){
        vftff((complex_s32_t*)vD, shift_mode);
    } else{
        assert(0);
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


void VFTFB(int8_t * vD, int mode, int shift_mode){
    if (mode == 16){
        vftfb_quake_s16((complex_s16_t*)vD, shift_mode);
    } else if (mode == 32){
        vftfb((complex_s32_t*)vD, shift_mode);
    } else{
        assert(0);
    }
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


void VFTTF(int8_t * vD, int mode, int shift_mode){
    if (mode == 16){
        vfttf_quake_s16((complex_s16_t*)vD, shift_mode);
    } else if (mode == 32){
        vfttf((complex_s32_t*)vD, shift_mode);
    } else{
        assert(0);
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

void VFTTB(int8_t * vD, int mode, int shift_mode){
    if (mode == 16){
        vfttb_quake_s16((complex_s16_t*)vD, shift_mode);
    } else if (mode == 32){
        vfttb((complex_s32_t*)vD, shift_mode);
    } else{
        assert(0);
    }
}

