#include <stdint.h>

//operates on vD
void VADSB(int8_t * output_reg, const int8_t * input_reg, int shift_mode);
void VFTTF(int8_t * vD, int mode, int shift_mode);
void VFTTB(int8_t * vD, int mode, int shift_mode);

//operates on vR
void VSBAD(int8_t * output_reg, const int8_t * input_reg, int shift_mode);
void VFTFF(int8_t * vR, int mode, int shift_mode);
void VFTFB(int8_t * vR, int mode, int shift_mode);

void VCMR_0(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, int mode);
void VCMR_1(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);
void VCMR_2(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);
void VCMR_3(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);

void VCMI_0(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, int mode);
void VCMI_1(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);
void VCMI_2(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);
void VCMI_3(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);

void VCMCR_0(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, int mode);
void VCMCR_1(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);
void VCMCR_2(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);
void VCMCR_3(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);

void VCMCI_0(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, int mode);
void VCMCI_1(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);
void VCMCI_2(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);
void VCMCI_3(int8_t * vR_out, const int8_t * vD_in, const int8_t * vC_in, const int8_t * vR_in, int mode);

void VLADSB(int8_t * vD_out, int8_t * vR_out, const int8_t * Mem_in, const int8_t * vR_in, int mode, int shift_mode);

// void VLMUL(int8_t * vR_out, int8_t * vR_in, int8_t * Mem_in, int mode);
// void VLMACC(int8_t * vD_out, int8_t * vR_out, int8_t * vD_in, int8_t * vR_in,int8_t * Mem_in, int mode);
// void VLMACCR(int8_t * vD_out, int8_t * vR_out, int8_t * vD_in, int8_t * vR_in,int8_t * Mem_in, int mode);

