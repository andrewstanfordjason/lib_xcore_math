#include <stdint.h>
#include <assert.h>
#include <string.h>

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

#define SAT16(VAL)      (((VAL) >= INT16_MAX)? INT16_MAX : (((VAL) <= INT16_MIN)? INT16_MIN : (VAL)))
#define SAT32(VAL)      (((VAL) >= INT32_MAX)? INT32_MAX : (((VAL) <= INT32_MIN)? INT32_MIN : (VAL)))

#define SHR(VAL, SHR)   ((VAL)>>SHR)
#define ASHR16(VAL, SHR_BITS)   (SAT16(((SHR_BITS) >= 0)? SHR((VAL),(SHR_BITS)) : (SAT16( ((int32_t)(VAL))<<(-(SHR_BITS)) ) )))
#define ASHR32(VAL, SHR_BITS)   (SAT32(((SHR_BITS) >= 0)? SHR((VAL),(SHR_BITS)) : (SAT32( ((int64_t)(VAL))<<(-(SHR_BITS)) ) )))
#define ASHR(BITS)  ASHR##BITS

static void vftff(
    int8_t * v,
    const int shift_mode)
{

    int64_t s[8];

    s[0] =  (int64_t)((int32_t *)v)[0] + ((int32_t *)v)[4];
    s[1] =  (int64_t)((int32_t *)v)[1] + ((int32_t *)v)[5];
    s[2] =  (int64_t)((int32_t *)v)[2] + ((int32_t *)v)[6];
    s[3] =  (int64_t)((int32_t *)v)[3] + ((int32_t *)v)[7];
    s[4] =  (int64_t)((int32_t *)v)[0] - ((int32_t *)v)[4];
    s[5] =  (int64_t)((int32_t *)v)[1] - ((int32_t *)v)[5];
    s[6] =  (int64_t)((int32_t *)v)[3] - ((int32_t *)v)[7];
    s[7] =  (int64_t)((int32_t *)v)[6] - ((int32_t *)v)[2];

    ((int32_t *)v)[0] = (int32_t) ASHR(32)(s[0] + s[2], shift_mode);
    ((int32_t *)v)[1] = (int32_t) ASHR(32)(s[1] + s[3], shift_mode);
    ((int32_t *)v)[2] = (int32_t) ASHR(32)(s[0] - s[2], shift_mode);
    ((int32_t *)v)[3] = (int32_t) ASHR(32)(s[1] - s[3], shift_mode);
    ((int32_t *)v)[4] = (int32_t) ASHR(32)(s[4] + s[6], shift_mode);
    ((int32_t *)v)[5] = (int32_t) ASHR(32)(s[5] + s[7], shift_mode);
    ((int32_t *)v)[6] = (int32_t) ASHR(32)(s[4] - s[6], shift_mode);
    ((int32_t *)v)[7] = (int32_t) ASHR(32)(s[5] - s[7], shift_mode);
}

static void vftff_quake_s16(
    int8_t * v,
    const int shift_mode)
{

    for (int i=0; i<8; i+=4){

        int32_t s[8];

        s[0] =  (int32_t)((int16_t *)v)[2*i+0] + ((int16_t *)v)[2*i+4];
        s[1] =  (int32_t)((int16_t *)v)[2*i+1] + ((int16_t *)v)[2*i+5];
        s[2] =  (int32_t)((int16_t *)v)[2*i+2] + ((int16_t *)v)[2*i+6];
        s[3] =  (int32_t)((int16_t *)v)[2*i+3] + ((int16_t *)v)[2*i+7];
        s[4] =  (int32_t)((int16_t *)v)[2*i+0] - ((int16_t *)v)[2*i+4];
        s[5] =  (int32_t)((int16_t *)v)[2*i+1] - ((int16_t *)v)[2*i+5];
        s[6] =  (int32_t)((int16_t *)v)[2*i+3] - ((int16_t *)v)[2*i+7];
        s[7] =  (int32_t)((int16_t *)v)[2*i+6] - ((int16_t *)v)[2*i+2];

        ((int16_t *)v)[2*i+0] = (int16_t) ASHR(16)(s[0] + s[2], shift_mode);
        ((int16_t *)v)[2*i+1] = (int16_t) ASHR(16)(s[1] + s[3], shift_mode);
        ((int16_t *)v)[2*i+2] = (int16_t) ASHR(16)(s[0] - s[2], shift_mode);
        ((int16_t *)v)[2*i+3] = (int16_t) ASHR(16)(s[1] - s[3], shift_mode);
        ((int16_t *)v)[2*i+4] = (int16_t) ASHR(16)(s[4] + s[6], shift_mode);
        ((int16_t *)v)[2*i+5] = (int16_t) ASHR(16)(s[5] + s[7], shift_mode);
        ((int16_t *)v)[2*i+6] = (int16_t) ASHR(16)(s[4] - s[6], shift_mode);
        ((int16_t *)v)[2*i+7] = (int16_t) ASHR(16)(s[5] - s[7], shift_mode);
    }
}

void VFTFF(int8_t * vD, int mode, int shift_mode){
    if (mode == 16){
        vftff_quake_s16(vD, shift_mode);
    } else if (mode == 32){
        vftff(vD, shift_mode);
    } else{
        assert(0);
    }
}

static void vftfb(
    int8_t * v,
    const int shift_mode)
{

    int64_t s[8];

    s[0] = (int64_t)((int32_t *)v)[0] + (int64_t)((int32_t *)v)[4];
    s[1] = (int64_t)((int32_t *)v)[1] + (int64_t)((int32_t *)v)[5];
    s[2] = (int64_t)((int32_t *)v)[2] + (int64_t)((int32_t *)v)[6];
    s[3] = (int64_t)((int32_t *)v)[3] + (int64_t)((int32_t *)v)[7];
    s[4] = (int64_t)((int32_t *)v)[0] - (int64_t)((int32_t *)v)[4];
    s[5] = (int64_t)((int32_t *)v)[1] - (int64_t)((int32_t *)v)[5];
    s[6] = (int64_t)((int32_t *)v)[7] - (int64_t)((int32_t *)v)[3];
    s[7] = (int64_t)((int32_t *)v)[2] - (int64_t)((int32_t *)v)[6];

    ((int32_t *)v)[0] = (int32_t) ASHR(32)(s[0] + s[2], shift_mode);
    ((int32_t *)v)[1] = (int32_t) ASHR(32)(s[1] + s[3], shift_mode);
    ((int32_t *)v)[2] = (int32_t) ASHR(32)(s[0] - s[2], shift_mode);
    ((int32_t *)v)[3] = (int32_t) ASHR(32)(s[1] - s[3], shift_mode);
    ((int32_t *)v)[4] = (int32_t) ASHR(32)(s[4] + s[6], shift_mode);
    ((int32_t *)v)[5] = (int32_t) ASHR(32)(s[5] + s[7], shift_mode);
    ((int32_t *)v)[6] = (int32_t) ASHR(32)(s[4] - s[6], shift_mode);
    ((int32_t *)v)[7] = (int32_t) ASHR(32)(s[5] - s[7], shift_mode);
}

static void vftfb_quake_s16(
    int8_t * v,
    const int shift_mode)
{

    for (int i=0; i<8; i+=4){

        int32_t s[8];

        s[0] = (int32_t)((int16_t *)v)[2*i+0] + (int32_t)((int16_t *)v)[2*i+4];
        s[1] = (int32_t)((int16_t *)v)[2*i+1] + (int32_t)((int16_t *)v)[2*i+5];
        s[2] = (int32_t)((int16_t *)v)[2*i+2] + (int32_t)((int16_t *)v)[2*i+6];
        s[3] = (int32_t)((int16_t *)v)[2*i+3] + (int32_t)((int16_t *)v)[2*i+7];
        s[4] = (int32_t)((int16_t *)v)[2*i+0] - (int32_t)((int16_t *)v)[2*i+4];
        s[5] = (int32_t)((int16_t *)v)[2*i+1] - (int32_t)((int16_t *)v)[2*i+5];
        s[6] = (int32_t)((int16_t *)v)[2*i+7] - (int32_t)((int16_t *)v)[2*i+3];
        s[7] = (int32_t)((int16_t *)v)[2*i+2] - (int32_t)((int16_t *)v)[2*i+6];

        ((int16_t *)v)[2*i+0] = (int16_t) ASHR(16)(s[0] + s[2], shift_mode);
        ((int16_t *)v)[2*i+1] = (int16_t) ASHR(16)(s[1] + s[3], shift_mode);
        ((int16_t *)v)[2*i+2] = (int16_t) ASHR(16)(s[0] - s[2], shift_mode);
        ((int16_t *)v)[2*i+3] = (int16_t) ASHR(16)(s[1] - s[3], shift_mode);
        ((int16_t *)v)[2*i+4] = (int16_t) ASHR(16)(s[4] + s[6], shift_mode);
        ((int16_t *)v)[2*i+5] = (int16_t) ASHR(16)(s[5] + s[7], shift_mode);
        ((int16_t *)v)[2*i+6] = (int16_t) ASHR(16)(s[4] - s[6], shift_mode);
        ((int16_t *)v)[2*i+7] = (int16_t) ASHR(16)(s[5] - s[7], shift_mode);
    }
}


void VFTFB(int8_t * vD, int mode, int shift_mode){
    if (mode == 16){
        vftfb_quake_s16(vD, shift_mode);
    } else if (mode == 32){
        vftfb(vD, shift_mode);
    } else{
        assert(0);
    }
}

static void vfttf(
    int8_t * v,
    const int shift_mode)
{

        int64_t s[8];

        s[0] = (int64_t)((int32_t *)v)[0] + (int64_t)((int32_t *)v)[2];
        s[1] = (int64_t)((int32_t *)v)[1] + (int64_t)((int32_t *)v)[3];
        s[2] = (int64_t)((int32_t *)v)[0] - (int64_t)((int32_t *)v)[2];
        s[3] = (int64_t)((int32_t *)v)[1] - (int64_t)((int32_t *)v)[3];
        s[4] = (int64_t)((int32_t *)v)[4] + (int64_t)((int32_t *)v)[6];
        s[5] = (int64_t)((int32_t *)v)[5] + (int64_t)((int32_t *)v)[7];
        s[6] = (int64_t)((int32_t *)v)[5] - (int64_t)((int32_t *)v)[7];
        s[7] = (int64_t)((int32_t *)v)[6] - (int64_t)((int32_t *)v)[4];

        ((int32_t *)v)[0] = (int32_t) ASHR(32)(s[0] + s[4], shift_mode);
        ((int32_t *)v)[1] = (int32_t) ASHR(32)(s[1] + s[5], shift_mode);
        ((int32_t *)v)[2] = (int32_t) ASHR(32)(s[2] + s[6], shift_mode);
        ((int32_t *)v)[3] = (int32_t) ASHR(32)(s[3] + s[7], shift_mode);
        ((int32_t *)v)[4] = (int32_t) ASHR(32)(s[0] - s[4], shift_mode);
        ((int32_t *)v)[5] = (int32_t) ASHR(32)(s[1] - s[5], shift_mode);
        ((int32_t *)v)[6] = (int32_t) ASHR(32)(s[2] - s[6], shift_mode);
        ((int32_t *)v)[7] = (int32_t) ASHR(32)(s[3] - s[7], shift_mode);
}


static void vfttf_quake_s16(
    int8_t * v,
    const int shift_mode)
{
    for (int i=0; i<8; i+=4){

        int32_t s[8];

        s[0] = (int32_t)((int16_t *)v)[2*i+0] + (int32_t)((int16_t *)v)[2*i+2];
        s[1] = (int32_t)((int16_t *)v)[2*i+1] + (int32_t)((int16_t *)v)[2*i+3];
        s[2] = (int32_t)((int16_t *)v)[2*i+0] - (int32_t)((int16_t *)v)[2*i+2];
        s[3] = (int32_t)((int16_t *)v)[2*i+1] - (int32_t)((int16_t *)v)[2*i+3];
        s[4] = (int32_t)((int16_t *)v)[2*i+4] + (int32_t)((int16_t *)v)[2*i+6];
        s[5] = (int32_t)((int16_t *)v)[2*i+5] + (int32_t)((int16_t *)v)[2*i+7];
        s[6] = (int32_t)((int16_t *)v)[2*i+5] - (int32_t)((int16_t *)v)[2*i+7];
        s[7] = (int32_t)((int16_t *)v)[2*i+6] - (int32_t)((int16_t *)v)[2*i+4];

        ((int16_t *)v)[2*i+0] = (int16_t) ASHR(16)(s[0] + s[4], shift_mode);
        ((int16_t *)v)[2*i+1] = (int16_t) ASHR(16)(s[1] + s[5], shift_mode);
        ((int16_t *)v)[2*i+2] = (int16_t) ASHR(16)(s[2] + s[6], shift_mode);
        ((int16_t *)v)[2*i+3] = (int16_t) ASHR(16)(s[3] + s[7], shift_mode);
        ((int16_t *)v)[2*i+4] = (int16_t) ASHR(16)(s[0] - s[4], shift_mode);
        ((int16_t *)v)[2*i+5] = (int16_t) ASHR(16)(s[1] - s[5], shift_mode);
        ((int16_t *)v)[2*i+6] = (int16_t) ASHR(16)(s[2] - s[6], shift_mode);
        ((int16_t *)v)[2*i+7] = (int16_t) ASHR(16)(s[3] - s[7], shift_mode);
    }
}


void VFTTF(int8_t * vD, int mode, int shift_mode){
    if (mode == 16){
        vfttf_quake_s16(vD, shift_mode);
    } else if (mode == 32){
        vfttf(vD, shift_mode);
    } else{
        assert(0);
    }
}

static void vfttb(
    int8_t * v,
    const int shift_mode)
{
    int64_t s[8];

    s[0] = (int64_t)((int32_t *)v)[0] + (int64_t)((int32_t *)v)[2];
    s[1] = (int64_t)((int32_t *)v)[1] + (int64_t)((int32_t *)v)[3];
    s[2] = (int64_t)((int32_t *)v)[0] - (int64_t)((int32_t *)v)[2];
    s[3] = (int64_t)((int32_t *)v)[1] - (int64_t)((int32_t *)v)[3];
    s[4] = (int64_t)((int32_t *)v)[4] + (int64_t)((int32_t *)v)[6];
    s[5] = (int64_t)((int32_t *)v)[5] + (int64_t)((int32_t *)v)[7];
    s[6] = (int64_t)((int32_t *)v)[7] - (int64_t)((int32_t *)v)[5];
    s[7] = (int64_t)((int32_t *)v)[4] - (int64_t)((int32_t *)v)[6];

    ((int32_t *)v)[0] = (int32_t) ASHR(32)(s[0] + s[4], shift_mode);
    ((int32_t *)v)[1] = (int32_t) ASHR(32)(s[1] + s[5], shift_mode);
    ((int32_t *)v)[2] = (int32_t) ASHR(32)(s[2] + s[6], shift_mode);
    ((int32_t *)v)[3] = (int32_t) ASHR(32)(s[3] + s[7], shift_mode);
    ((int32_t *)v)[4] = (int32_t) ASHR(32)(s[0] - s[4], shift_mode);
    ((int32_t *)v)[5] = (int32_t) ASHR(32)(s[1] - s[5], shift_mode);
    ((int32_t *)v)[6] = (int32_t) ASHR(32)(s[2] - s[6], shift_mode);
    ((int32_t *)v)[7] = (int32_t) ASHR(32)(s[3] - s[7], shift_mode);
}

static void vfttb_quake_s16(
    int8_t * v,
    const int shift_mode)
{

    for (int i=0; i<8; i+=4){

        int64_t s[8];

        s[0] = (int32_t)((int16_t *)v)[2*i+0] + (int32_t)((int16_t *)v)[2*i+2];
        s[1] = (int32_t)((int16_t *)v)[2*i+1] + (int32_t)((int16_t *)v)[2*i+3];
        s[2] = (int32_t)((int16_t *)v)[2*i+0] - (int32_t)((int16_t *)v)[2*i+2];
        s[3] = (int32_t)((int16_t *)v)[2*i+1] - (int32_t)((int16_t *)v)[2*i+3];
        s[4] = (int32_t)((int16_t *)v)[2*i+4] + (int32_t)((int16_t *)v)[2*i+6];
        s[5] = (int32_t)((int16_t *)v)[2*i+5] + (int32_t)((int16_t *)v)[2*i+7];
        s[6] = (int32_t)((int16_t *)v)[2*i+7] - (int32_t)((int16_t *)v)[2*i+5];
        s[7] = (int32_t)((int16_t *)v)[2*i+4] - (int32_t)((int16_t *)v)[2*i+6];

        ((int16_t *)v)[2*i+0] = (int16_t) ASHR(16)(s[0] + s[4], shift_mode);
        ((int16_t *)v)[2*i+1] = (int16_t) ASHR(16)(s[1] + s[5], shift_mode);
        ((int16_t *)v)[2*i+2] = (int16_t) ASHR(16)(s[2] + s[6], shift_mode);
        ((int16_t *)v)[2*i+3] = (int16_t) ASHR(16)(s[3] + s[7], shift_mode);
        ((int16_t *)v)[2*i+4] = (int16_t) ASHR(16)(s[0] - s[4], shift_mode);
        ((int16_t *)v)[2*i+5] = (int16_t) ASHR(16)(s[1] - s[5], shift_mode);
        ((int16_t *)v)[2*i+6] = (int16_t) ASHR(16)(s[2] - s[6], shift_mode);
        ((int16_t *)v)[2*i+7] = (int16_t) ASHR(16)(s[3] - s[7], shift_mode);
    }
}

void VFTTB(int8_t * vD, int mode, int shift_mode){
    if (mode == 16){
        vfttb_quake_s16(vD, shift_mode);
    } else if (mode == 32){
        vfttb(vD, shift_mode);
    } else{
        assert(0);
    }
}

