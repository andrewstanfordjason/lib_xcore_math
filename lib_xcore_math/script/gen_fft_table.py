# Copyright 2020-2022 XMOS LIMITED.
# This Software is subject to the terms of the XMOS Public Licence: Version 1.
from contextlib import suppress
import numpy as np
import argparse
import os


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--out_file",
        default="xmath_fft_lut",
        help="Filename to be used (with '.h' and '.c') for the generated files. (default: 'xmath_fft_lut')",
    )
    parser.add_argument(
        "--out_dir", default="./", help="Directory to output generated files to."
    )
    parser.add_argument(
        "--max_fft_log2",
        type=int,
        default=5,
        help="Log2 of the maximum FFT size supported. (default: 5)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Enable verbose mode.",
    )
    parser.add_argument(
        "--dit",
        action="store_true",
        default=False,
        help="Generate decimation-in-time tables.",
    )
    parser.add_argument(
        "--dif",
        action="store_true",
        default=False,
        help="Generate decimation-in-frequency tables.",
    )

    args = parser.parse_args()

    if not (args.dit or args.dif):
        print("Neither the --dit nor --dif flags were provided. No work to be done.")
        return

    header_filename = f"{args.out_file}.h"
    source_filename = f"{args.out_file}.c"

    if args.verbose:
        print(f"Header filename: {header_filename}")
        print(f"Source filename: {source_filename}")

    with open(os.path.join(args.out_dir, source_filename), "w+") as source_file:
        with open(os.path.join(args.out_dir, header_filename), "w+") as header_file:

            header_file.write(
                "// Copyright 2021 XMOS LIMITED. This Software is subject to the terms of the\n"
            )
            header_file.write("// XMOS Public License: Version 1\n")
            header_file.write("#pragma once\n")

            header_file.write('#include "xmath/xmath.h"\n\n')

            source_file.write(f'#include "{header_filename}"\n\n')
            
            for short in [
                True, 
                False
                ]:
                print('*'*80)
                if short:
                    output_dtype = np.int16
                    define_suffix = '_16' #this keeps it backwards compatable for now
                else:
                    output_dtype = np.int32
                    define_suffix =''
                    

                if args.dit:

                    if args.verbose:
                        print("Generating DIT tables..")

                    table_name = generate_DIT_FFT(
                        args.max_fft_log2, header_file, source_file, args, M=32, output_dtype = output_dtype
                    )

                    header_file.write(
                        "\n/** @brief Maximum FFT length (log2) that can be performed using decimation-in-time. */\n"
                    )
                    header_file.write(f"#define MAX_DIT_FFT_LOG2" + str(define_suffix) + " " + str(args.max_fft_log2) + "\n")
                    header_file.write(
                        "\n/** @brief Convenience macro to index into the decimation-in-time FFT look-up table. \n\n"
                    )
                    header_file.write(
                        "\tThis will return the address at which the coefficients for the final pass of the real DIT\n"
                    )
                    header_file.write("\tFFT algorithm begin. \n\n")
                    # header_file.write("\t\n")
                    header_file.write("\t@param N\tThe FFT length.\n*/\n")
                    header_file.write(
                        "#define XMATH_DIT_REAL_FFT_LUT" + str(define_suffix) + "(N) &" + str(table_name) + "[(N)-8]\n\n"
                    )

                if args.dif:

                    if args.verbose:
                        print("Generating DIF tables..")

                    table_name = generate_DIF_FFT(
                        args.max_fft_log2, header_file, source_file, args, M=32, output_dtype = output_dtype
                    )

                    header_file.write(
                        "\n/** @brief Maximum FFT length (log2) that can be performed using decimation-in-frequency. */\n"
                    )
                    header_file.write(f"#define MAX_DIF_FFT_LOG2" + str(define_suffix) + " " + str(args.max_fft_log2) + "\n")
                    header_file.write(
                        "\n/** @brief Convenience macro to index into the decimation-in-frequency FFT look-up table. \n\n"
                    )
                    header_file.write(
                        "\tUse this to retrieve the correct address for the DIF FFT look-up table when performing\n"
                    )
                    header_file.write(
                        "\tan FFT (or IFFT) using the DIF algorithm. (@see fft_dif_forward).\n\n"
                    )
                    header_file.write("\t@param N\tThe FFT length.\n*/\n")
                    header_file.write(
                        "#define XMATH_DIF_FFT_LUT" + str(define_suffix) + "(N) &" + str(table_name) + "[(1<<(MAX_DIF_FFT_LOG2)) - (N)]\n\n"
                    )

def generate_DIF_FFT(N, header_file, source_file, args, M = 8 * 4, output_dtype = np.int32, base_name = 'xmath_dif_fft_lut_'):
    print('dif')
    return generate_FFT(N, header_file, source_file, args, base_name, populate_twiddle_factors_dif, M = M, output_dtype = output_dtype)

def generate_DIT_FFT(N, header_file, source_file, args, M = 8 * 4, output_dtype = np.int32, base_name = 'xmath_dit_fft_lut_'):
    print('dit')
    return generate_FFT(N, header_file, source_file, args, base_name, populate_twiddle_factors_dit, M = M, output_dtype = output_dtype)

def generate_FFT(N, header_file, source_file, args, base_name, gen_fn, M = 8 * 4, output_dtype = np.int32):
    
    bytes_per_twiddle = output_dtype(0).nbytes*2
    twiddles_per_load = M // bytes_per_twiddle

    twiddle_factor_count = (2 ** N)//twiddles_per_load - 1

    if output_dtype == np.int16:
        twiddle_factor_count += 1
        #this one is for the half empty one

    twiddle_factors = np.zeros((twiddle_factor_count, twiddles_per_load), dtype=np.complex128)

    r = int(np.log2(twiddles_per_load))

    bits_per_byte = 8
    output_dtype_num_bits = (bytes_per_twiddle * bits_per_byte) // 2
    output_complex_type = 'complex_s' + str(output_dtype_num_bits) + '_t'

    table_name = base_name + str(output_dtype_num_bits)

    total_bytes = bytes_per_twiddle * twiddle_factor_count * twiddles_per_load

    header_file.write(
        f"extern const {output_complex_type} {table_name}[{twiddles_per_load*twiddle_factor_count}]; // {total_bytes} bytes\n"
    )

    gen_fn(N, M, twiddles_per_load, twiddle_factors, r)

    twiddle_exponent = -(output_dtype_num_bits - 2)
    r = twiddle_factors.view(np.float64).flatten()
    twiddle_factors_int = np.asarray(
        np.rint(r * 2 ** -twiddle_exponent), dtype=output_dtype
    )

    twiddle_factors_int = np.reshape(twiddle_factors_int, (twiddle_factor_count * twiddles_per_load, 2))

    source_file.write(
        f"const {output_complex_type} {table_name}[{twiddles_per_load*twiddle_factor_count}] = \n\t{{\n\t"
    )

    for t in range(len(twiddle_factors_int)):
        source_file.write(
            "{%11d, %11d}, " % (twiddle_factors_int[t][0], twiddle_factors_int[t][1])
        )
        if (t + 1) % twiddles_per_load == 0:
            source_file.write("\n\t")

    source_file.write("};\n\n")

    source_file.write(f"const uint32_t {table_name}_size = sizeof({table_name});\n\n")

    return table_name

def populate_twiddle_factors_dit(N, M, twiddles_per_load, twiddle_factors, r):
    P = int(np.log2(twiddles_per_load)) + 1
    a = M * 2 ** (N - P)
    b = M
    t = 0

    # insert special bit
    if P == 4: 
        half_tw = np.arange(twiddles_per_load//2)
        twiddle_stride = 2 ** (N - P + 1)
        q =  (twiddle_stride * half_tw)
        twiddle_factors[t, twiddles_per_load//2:] = np.exp(-2.0j * np.pi * q / (2 ** N))
        twiddle_factors[t, :twiddles_per_load//2] = 1

        # twiddle_factors[t, twiddles_per_load//2:] = twiddle_factors[t, :twiddles_per_load//2]
        t += 1

    for i in range(N - r):
        for k in range(b - M, -M, -M):
            twiddle_offset = (k / M) * (a / M) * twiddles_per_load
            twiddle_stride = 2 ** (N - P) / (b / M)

            tw = np.arange(twiddles_per_load)
            q = twiddle_offset + (twiddle_stride * tw)

            twiddle_factors[t, :] = np.exp(-2.0j * np.pi * q / (2 ** N))

            t += 1

        b *= 2
        a //= 2

def populate_twiddle_factors_dif(N, M, twiddles_per_load, twiddle_factors, r):
    P = int(np.log2(twiddles_per_load)) + 1
    a = M
    b = M * 2 ** (N - P)  # first two rounds are done by the fft inst

    t = 0
    for i in range(N - r):
        for k in range(b - M, -M, -M):
            twiddle_offset = (k / M) * (a / M) * twiddles_per_load
            twiddle_stride = 2 ** (N - P) / (b / M)

            tw = np.arange(twiddles_per_load)
            q = twiddle_offset + (twiddle_stride * tw)

            twiddle_factors[t, :] = -np.exp(-2.0j * np.pi * q / (2 ** N))
            t += 1

        b //= 2
        a *= 2

    # insert special bit
    if P == 4: 
        half_tw = np.arange(twiddles_per_load//2)
        twiddle_stride = 2 ** (N - P+1)
        q =  (twiddle_stride * half_tw)
        twiddle_factors[t, twiddles_per_load//2:] = -np.exp(-2.0j * np.pi * q / (2 ** N))
        twiddle_factors[t, :twiddles_per_load//2] = 1
        t += 1

if __name__ == "__main__":
    main()

