#include <stdint.h>
#include <assert.h>
#include <string.h>

void VLMUL_0(int8_t * vR_out, int8_t * vD_out, const int8_t * vR_in, const int8_t * Mem_in, int mode){

    if (mode == 8){
        for (int i=0;i<32;i++){
            int64_t t = (int64_t)vR_in[i] *  (int64_t)Mem_in[i];
            vR_out[i] =  (t + (1<<6)) >> 7;
        }
    } else if (mode == 16){
        int16_t * vD_out_16 = (int16_t *)vD_out;
        const int16_t * vR_in_16 =(const int16_t * )vR_in;
        const int16_t * Mem_in_16 = (const int16_t * )Mem_in;

        for (int i=0;i<16;i++){
            int64_t t = (int64_t)Mem_in_16[i] * (uint8_t)vR_in_16[i];
            vD_out_16[i] = t >> 8;
        }
    } else if (mode == 32){
        int32_t * vD_out_32 = (int32_t *)vD_out;
        const int32_t * vR_in_32 =(const int32_t * )vR_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;

        for (int i=0;i<8;i++){
            int64_t t = (int64_t)Mem_in_32[i] *  (uint8_t)vR_in_32[i];
            vD_out_32[i] = t >> 8;
        }
    } else{
        assert(0);
    }
}

void VLMUL_1(int8_t * vR_out, int8_t * vD_out, const int8_t * vD_in, const int8_t * vR_in, const int8_t * Mem_in, int mode){

    //Note does it cost any more to have int8 * uint8 for this mode?
    if (mode == 16){
        int16_t * vR_out_16 = (int16_t *)vR_out;
        const int16_t * vR_in_16 =(const int16_t * )vR_in;
        const int16_t * vD_in_16 =(const int16_t * )vD_in;
        const int16_t * Mem_in_16 = (const int16_t * )Mem_in;

        for (int i=0;i<16;i++){
            int64_t t = (int64_t)Mem_in_16[i] *  (int8_t)(vR_in_16[i]>>8);
            t  += ((int64_t)vD_in_16[i]);
            vR_out_16[i] =  (t + (1<<6)) >> 7; //Note, cannot overflow
        }

    } else if (mode == 32){
        int32_t * vD_out_32 = (int32_t *)vD_out;
        const int32_t * vD_in_32 =(const int32_t * )vD_in;
        const int32_t * vR_in_32 =(const int32_t * )vR_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;

        for (int i=0;i<8;i++){
            int64_t t = (int64_t)Mem_in_32[i] *  (uint8_t)(vR_in_32[i]>>8);
            t += vD_in_32[i];
            vD_out_32[i] = t >> 8;
        }
    } else{
        assert(0);
    }
}

void VLMUL_2(int8_t * vR_out, int8_t * vD_out, const int8_t * vD_in, const int8_t * vR_in, const int8_t * Mem_in, int mode){

    //[asj:removed unless semantics are worked out]
    // if (mode == 16){
    //     // int16_t * vD_out_16 = (int16_t *)vD_out;
    //     int16_t * vR_out_16 = (int16_t * )vR_out;
    //     const int16_t * vR_in_16 =(const int16_t * )vR_in;
    //     const int16_t * Mem_in_16 = (const int16_t * )Mem_in;

    //     for (int i=0;i<16;i++){
    //         int64_t t = (int64_t)Mem_in_16[i] * (int8_t)vR_in_16[i];
    //         vR_out_16[i] =  (t + (1<<6)) >> 7;
    //     }

    // } else 
    if (mode == 32){
        int32_t * vD_out_32 = (int32_t *)vD_out;
        const int32_t * vD_in_32 =(const int32_t * )vD_in;
        const int32_t * vR_in_32 =(const int32_t * )vR_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;

        for (int i=0;i<8;i++){
            int64_t t = (int64_t)Mem_in_32[i] *  (uint8_t)(vR_in_32[i]>>16);
            t += vD_in_32[i];
            vD_out_32[i] = t >> 8;
        }
    } else{
        assert(0);
    }
}

void VLMUL_3(int8_t * vR_out, int8_t * vD_out, const int8_t * vD_in, const int8_t * vR_in, const int8_t * Mem_in, int mode){

    //[asj:removed unless semantics are worked out]
    // if (mode == 16){
    //     int16_t * vR_out_16 = (int16_t * )vR_out;
    //     const uint16_t * vR_in_16 =(const uint16_t * )vR_in;
    //     const int16_t * vD_in_16 =(const int16_t * )vD_in;
    //     const int16_t * Mem_in_16 = (const int16_t * )Mem_in;

    //     for (int i=0;i<16;i++){
    //         int64_t t = (int64_t)Mem_in_16[i] * (uint8_t)(vR_in_16[i]>>8);
    //         t += vD_in_16[i];
    //         vR_out_16[i] = (t + (1<<7)) >> 8;
    //     }
    // } else 
    if (mode == 32){
        int32_t * vR_out_32 = (int32_t *)vR_out;
        const int32_t * vD_in_32 =(const int32_t * )vD_in;
        const int32_t * vR_in_32 =(const int32_t * )vR_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;

        for (int i=0;i<8;i++){
            int64_t t = (int64_t)Mem_in_32[i] *  (int8_t)(vR_in_32[i]>>24);
            t += vD_in_32[i];
            vR_out_32[i] =  (t + (1<<6)) >> 7;//Note, cannot overflow
        }
    } else{
        assert(0);
    }
}

void VLMACC_0(int8_t * vD_out, int8_t * vR_out, int8_t * vD_in, int8_t * vR_in, int8_t * Mem_in, int8_t * vC_in, int mode){
    
    if (mode == 8){
        int16_t * vR_out_16 = (int16_t * )vR_out;
        int16_t * vD_out_16 = (int16_t * )vD_out;
        const int16_t * vR_in_16 = (const int16_t * )vR_in;
        const int16_t * vD_in_16 = (const int16_t * )vD_in;
        
        for (int i=0;i<16;i++){

            int32_t accu = 0;
            ((int16_t*)&accu)[0] = vR_in_16[i];
            ((int16_t*)&accu)[1] = vD_in_16[i];

            int64_t t = (int64_t)Mem_in[i] * (int8_t)vC_in[i];

            t += accu;

            if (t > INT32_MAX) accu = INT32_MAX;
            else if (t < INT32_MIN) accu = INT32_MIN;
            else accu = t;

            vR_out_16[i] = ((int16_t*)&accu)[0];
            vD_out_16[i] = ((int16_t*)&accu)[1];
        }

    } else if (mode == 16){
        const int16_t * vR_in_16 = (const int16_t * )vR_in;
        const int16_t * vD_in_16 = (const int16_t * )vD_in;
        const int16_t * vC_in_16 = (const int16_t * )vC_in;
        const int16_t * Mem_in_16 = (const int16_t * )Mem_in;

        int16_t * vR_out_16 = (int16_t * )vR_out;
        int16_t * vD_out_16 = (int16_t * )vD_out;
        
        for (int i=0;i<16;i++){

            int32_t accu = 0;
            ((int16_t*)&accu)[0] = vR_in_16[i];
            ((int16_t*)&accu)[1] = vD_in_16[i];

            int64_t t = (int64_t)Mem_in_16[i] * (uint8_t)vC_in_16[i];

            t += accu;
            accu = (t+0) >> 1;

            vR_out_16[i] = ((int16_t*)&accu)[0];
            vD_out_16[i] = ((int16_t*)&accu)[1];
        }

    } else if (mode == 32){
        const int32_t * vR_in_32 = (const int32_t * )vR_in;
        const int32_t * vD_in_32 = (const int32_t * )vD_in;
        const int32_t * vC_in_32 = (const int32_t * )vC_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;
        int32_t * vR_out_32 = (int32_t * )vR_out;
        int32_t * vD_out_32 = (int32_t * )vD_out;
        
        for (int i=0;i<8;i++){

            int64_t accu = 0;
            ((int32_t*)&accu)[0] = vR_in_32[i];
            ((int32_t*)&accu)[1] = vD_in_32[i];

            int64_t t = (int64_t)Mem_in_32[i] * (uint8_t)vC_in_32[i];

            #if 0
            accu >>= 1;
            accu += (t>>1);
            #else
            accu = ((__int128_t)accu + (__int128_t)t +1)>>1;
            #endif
            vR_out_32[i] = ((int32_t*)&accu)[0];
            vD_out_32[i] = ((int32_t*)&accu)[1];
        }
    } else{
        assert(0);
    }
}

void VLMACC_1(int8_t * vD_out, int8_t * vR_out, int8_t * vD_in, int8_t * vR_in, int8_t * Mem_in, int8_t * vC_in, int mode){
    if (mode == 16){
        const int16_t * vR_in_16 = (const int16_t * )vR_in;
        const int16_t * vD_in_16 = (const int16_t * )vD_in;
        const int16_t * vC_in_16 = (const int16_t * )vC_in;
        const int16_t * Mem_in_16 = (const int16_t * )Mem_in;

        int16_t * vR_out_16 = (int16_t * )vR_out;
        int16_t * vD_out_16 = (int16_t * )vD_out;
        
        for (int i=0;i<16;i++){

            int32_t accu = 0;
            ((int16_t*)&accu)[0] = vR_in_16[i];
            ((int16_t*)&accu)[1] = vD_in_16[i];

            int64_t t = (int64_t)Mem_in_16[i] * (int8_t)(vC_in_16[i]>>8);

            t <<= 8;
            t += (int64_t)accu << 1;

            if (t > INT32_MAX) accu = INT32_MAX;
            else if (t < INT32_MIN) accu = INT32_MIN;
            else accu = t;

            vR_out_16[i] = ((int16_t*)&accu)[0];
            vD_out_16[i] = ((int16_t*)&accu)[1];
        }

    } else if (mode == 32){
        const int32_t * vR_in_32 = (const int32_t * )vR_in;
        const int32_t * vD_in_32 = (const int32_t * )vD_in;
        const int32_t * vC_in_32 = (const int32_t * )vC_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;
        int32_t * vR_out_32 = (int32_t * )vR_out;
        int32_t * vD_out_32 = (int32_t * )vD_out;
        
        for (int i=0;i<8;i++){

            int64_t accu = 0;
            ((int32_t*)&accu)[0] = vR_in_32[i];
            ((int32_t*)&accu)[1] = vD_in_32[i];

            int64_t t = (int64_t)Mem_in_32[i] * (uint8_t)(vC_in_32[i]>>8);

            #if 0
            accu += (t<<(8 - 1));
            #else
            accu = ((((__int128_t)accu)<<1) + ((__int128_t)t<<8) +1)>>1;
            #endif
            vR_out_32[i] = ((int32_t*)&accu)[0];
            vD_out_32[i] = ((int32_t*)&accu)[1];
        }
    } else{
        assert(0);
    }
}


void VLMACC_2(int8_t * vD_out, int8_t * vR_out, int8_t * vD_in, int8_t * vR_in, int8_t * Mem_in, int8_t * vC_in, int mode){
    if (mode == 32){
        const int32_t * vR_in_32 = (const int32_t * )vR_in;
        const int32_t * vD_in_32 = (const int32_t * )vD_in;
        const int32_t * vC_in_32 = (const int32_t * )vC_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;
        int32_t * vR_out_32 = (int32_t * )vR_out;
        int32_t * vD_out_32 = (int32_t * )vD_out;
        
        for (int i=0;i<8;i++){

            int64_t accu = 0;
            ((int32_t*)&accu)[0] = vR_in_32[i];
            ((int32_t*)&accu)[1] = vD_in_32[i];

            int64_t t = (int64_t)Mem_in_32[i] * (uint8_t)(vC_in_32[i]>>16);

            #if 0
            accu += (t<<(16 - 1));
            #else
            accu = ((((__int128_t)accu)<<1) + ((__int128_t)t<<16) +1)>>1;
            #endif

            vR_out_32[i] = ((int32_t*)&accu)[0];
            vD_out_32[i] = ((int32_t*)&accu)[1];
        }
    } else{
        assert(0);
    }
}
void VLMACC_3(int8_t * vD_out, int8_t * vR_out, int8_t * vD_in, int8_t * vR_in, int8_t * Mem_in, int8_t * vC_in, int mode){
    if (mode == 32){
        const int32_t * vR_in_32 = (const int32_t * )vR_in;
        const int32_t * vD_in_32 = (const int32_t * )vD_in;
        const int32_t * vC_in_32 = (const int32_t * )vC_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;
        int32_t * vR_out_32 = (int32_t * )vR_out;
        int32_t * vD_out_32 = (int32_t * )vD_out;
        
        for (int i=0;i<8;i++){

            int64_t accu = 0;
            ((int32_t*)&accu)[0] = vR_in_32[i];
            ((int32_t*)&accu)[1] = vD_in_32[i];

            int64_t t = (int64_t)Mem_in_32[i] * (int8_t)(vC_in_32[i]>>24);
            
            #if 0
            accu += (t<<(24 - 1));
            if (accu > INT64_MAX/2) accu = INT64_MAX;
            else if (accu < INT64_MIN/2) accu = INT64_MIN;
            else accu <<= 1;
            #else
            __int128_t a = (((__int128_t)accu)<<1) + ((__int128_t)t<<24);
            if (a > INT64_MAX) accu = INT64_MAX;
            else if (a < INT64_MIN) accu = INT64_MIN;
            else accu = a;
            #endif

            vR_out_32[i] = ((int32_t*)&accu)[0];
            vD_out_32[i] = ((int32_t*)&accu)[1];
        }
    } else{
        assert(0);
    }
}














void VLMACCR_0(
    int8_t * vD_out, int8_t * vR_out, 
    const int8_t * vD_in, const int8_t * vR_in, 
    const int8_t * vC_in, const int8_t * Mem_in, int mode){

    if (mode == 8){
        int16_t * vR_out_16 = (int16_t * )vR_out;
        int16_t * vD_out_16 = (int16_t * )vD_out;
        const int16_t * vR_in_16 = (const int16_t * )vR_in;
        const int16_t * vD_in_16 = (const int16_t * )vD_in;
        
        int32_t accu = 0;
        ((int16_t*)&accu)[0] = vR_in_16[15];
        ((int16_t*)&accu)[1] = vD_in_16[15];

        int64_t a = accu;

        for (int i=0;i<32;i++)
            a += (int64_t)Mem_in[i] * (int8_t)vC_in[i];
        
        if (a > INT32_MAX) accu = INT32_MAX;
        else if (a < INT32_MIN) accu = INT32_MIN;
        else accu = a;

        // then rotate {vD:vR}
        for(int i=15;i>=1;i--){
            vR_out_16[i] = vR_in_16[i-1];
            vD_out_16[i] = vD_in_16[i-1];
        }
        vR_out_16[0] = ((int16_t*)&accu)[0];
        vD_out_16[0] = ((int16_t*)&accu)[1];

    } else if (mode == 16){
        const int16_t * vR_in_16 = (const int16_t * )vR_in;
        const int16_t * vD_in_16 = (const int16_t * )vD_in;
        const int16_t * vC_in_16 = (const int16_t * )vC_in;
        const int16_t * Mem_in_16 = (const int16_t * )Mem_in;

        int16_t * vR_out_16 = (int16_t * )vR_out;
        int16_t * vD_out_16 = (int16_t * )vD_out;
        
        int32_t accu = 0;
        ((int16_t*)&accu)[0] = vR_in_16[15];
        ((int16_t*)&accu)[1] = vD_in_16[15];

        int64_t t = 0;

        for (int i=0;i<16;i++){
            int64_t v = (int64_t)Mem_in_16[i] * (uint8_t)vC_in_16[i];
            t += v;
        }

        //t is a 16+8+4 = 28 bit number
        accu = ((__int128_t)accu + (__int128_t)t + 1)>>1;

        for(int i=15;i>=1;i--){
            vR_out_16[i] = vR_in_16[i-1];
            vD_out_16[i] = vD_in_16[i-1];
        }

        vR_out_16[0] = ((int16_t*)&accu)[0];
        vD_out_16[0] = ((int16_t*)&accu)[1];

    } else if (mode == 32){
        const int32_t * vR_in_32 = (const int32_t * )vR_in;
        const int32_t * vD_in_32 = (const int32_t * )vD_in;
        const int32_t * vC_in_32 = (const int32_t * )vC_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;
        int32_t * vR_out_32 = (int32_t * )vR_out;
        int32_t * vD_out_32 = (int32_t * )vD_out;
        
        int64_t accu = 0;
        ((int32_t*)&accu)[0] = vR_in_32[7];
        ((int32_t*)&accu)[1] = vD_in_32[7];

        int64_t t = 0;

        for (int i=0;i<8;i++){
            int64_t v = (int64_t)Mem_in_32[i] * (uint8_t)vC_in_32[i];
            t += v;
        }

        //t is a 32+8+4 = 44 bit number
        accu = ((__int128_t)accu + (__int128_t)t + 1)>>1;

        // then rotate {vD:vR}
        for(int i=7;i>=1;i--){
            vR_out_32[i] = vR_in_32[i-1];
            vD_out_32[i] = vD_in_32[i-1];
        }
        vR_out_32[0] = ((int32_t*)&accu)[0];
        vD_out_32[0] = ((int32_t*)&accu)[1];
    } else{
        assert(0);
    }
}

void VLMACCR_1(
    int8_t * vD_out, int8_t * vR_out, 
    const int8_t * vD_in, const int8_t * vR_in, 
    const int8_t * vC_in, const int8_t * Mem_in, int mode)
{
    if (mode == 16){
        const int16_t * vR_in_16 = (const int16_t * )vR_in;
        const int16_t * vD_in_16 = (const int16_t * )vD_in;
        const int16_t * vC_in_16 = (const int16_t * )vC_in;
        const int16_t * Mem_in_16 = (const int16_t * )Mem_in;

        int16_t * vR_out_16 = (int16_t * )vR_out;
        int16_t * vD_out_16 = (int16_t * )vD_out;
        
        int32_t accu = 0;
        ((int16_t*)&accu)[0] = vR_in_16[0];
        ((int16_t*)&accu)[1] = vD_in_16[0];

        int64_t t = 0;

        for (int i=0;i<16;i++){
            int64_t v = (int64_t)Mem_in_16[i] * (int8_t)(vC_in_16[i]>>8);
            t += v;
        }
        t <<= 8;
        //t is a 16+7 +8 + 3 = 34 bit number
        int64_t a = ((int64_t)accu << 1) + t;

        if (a > INT32_MAX) accu = INT32_MAX;
        else if (a < INT32_MIN) accu = INT32_MIN;
        else accu = a;

        vR_out_16[0] = ((int16_t*)&accu)[0];
        vD_out_16[0] = ((int16_t*)&accu)[1];

    } else if (mode == 32){
        const int32_t * vR_in_32 = (const int32_t * )vR_in;
        const int32_t * vD_in_32 = (const int32_t * )vD_in;
        const int32_t * vC_in_32 = (const int32_t * )vC_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;
        int32_t * vR_out_32 = (int32_t * )vR_out;
        int32_t * vD_out_32 = (int32_t * )vD_out;
        
        int64_t accu = 0;
        ((int32_t*)&accu)[0] = vR_in_32[0];
        ((int32_t*)&accu)[1] = vD_in_32[0];

        int64_t t = 0;

        for (int i=0;i<8;i++){
            int64_t v = (int64_t)Mem_in_32[i] * (uint8_t)(vC_in_32[i]>>8);
            t += v;
        }
        t <<= (8-1);
        //t is a 32+8 + 8 + 3 = 51 bit number
        accu += t; //no overflow

        vR_out_32[0] = ((int32_t*)&accu)[0];
        vD_out_32[0] = ((int32_t*)&accu)[1];
    } else{
        assert(0);
    }
}

void VLMACCR_2(
    int8_t * vD_out, int8_t * vR_out, 
    const int8_t * vD_in, const int8_t * vR_in, 
    const int8_t * vC_in, const int8_t * Mem_in, int mode)
{
    if (mode == 32){
        const int32_t * vR_in_32 = (const int32_t * )vR_in;
        const int32_t * vD_in_32 = (const int32_t * )vD_in;
        const int32_t * vC_in_32 = (const int32_t * )vC_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;
        int32_t * vR_out_32 = (int32_t * )vR_out;
        int32_t * vD_out_32 = (int32_t * )vD_out;
        
        int64_t accu = 0;
        ((int32_t*)&accu)[0] = vR_in_32[0];
        ((int32_t*)&accu)[1] = vD_in_32[0];

        int64_t t = 0;

        for (int i=0;i<8;i++){
            int64_t v = (int64_t)Mem_in_32[i] * (uint8_t)(vC_in_32[i]>>16);
            t += v;
        }
        t <<= (16-1);
        //t is a 32+8 + 16 + 3 = 59 bit number
        accu += t; //no overflow

        vR_out_32[0] = ((int32_t*)&accu)[0];
        vD_out_32[0] = ((int32_t*)&accu)[1];
    } else{
        assert(0);
    }
}

void VLMACCR_3(
    int8_t * vD_out, int8_t * vR_out, 
    const int8_t * vD_in, const int8_t * vR_in, 
    const int8_t * vC_in, const int8_t * Mem_in, int mode)
{
    if (mode == 32){
        const int32_t * vR_in_32 = (const int32_t * )vR_in;
        const int32_t * vD_in_32 = (const int32_t * )vD_in;
        const int32_t * vC_in_32 = (const int32_t * )vC_in;
        const int32_t * Mem_in_32 = (const int32_t * )Mem_in;
        int32_t * vR_out_32 = (int32_t * )vR_out;
        int32_t * vD_out_32 = (int32_t * )vD_out;
        
        int64_t accu = 0;
        ((int32_t*)&accu)[0] = vR_in_32[0];
        ((int32_t*)&accu)[1] = vD_in_32[0];

        __int128 t = 0;

        for (int i=0;i<8;i++){
            int64_t v = (int64_t)Mem_in_32[i] * (int8_t)(vC_in_32[i]>>24);
            t += v;
        }
        t <<= 24;
        //t is a 32 + 7 + 24 + 3 = 66 bit number

        __int128 a = (((__int128)accu) << 1) + t;

        if (a > INT64_MAX) accu = INT64_MAX;
        else if (a < INT64_MIN) accu = INT64_MIN;
        else accu = a;
        
        vR_out_32[0] = ((int32_t*)&accu)[0];
        vD_out_32[0] = ((int32_t*)&accu)[1];
    } else{
        assert(0);
    }
}