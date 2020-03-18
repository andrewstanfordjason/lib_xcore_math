
#ifndef XS3_BFP_H_
#define XS3_BFP_H_

#include "xs3_math_types.h"


#ifdef __XC__
extern "C" {
#endif

/*
    Saturation Logic:

    Where specified, the VPU-enhanced arithmetic functions apply symmetric saturation logic, which is applied
    where results would otherwise overflow.

    The resulting integer range for N-bit saturation is   -(2^(N-1))+1  to  (2^(N-1))-1.

    8-bit:      (-127, 127)
    16-bit:     (-65535, 65535)
    32-bit:     (-2147483647, 2147483647)
*/


/*
    Exponent Calculation / Headroom Minimization:

    In many of the arithmetic functions, an output exponent is chosen prior to actually performing element-by-element
    operations. The exponent is chosen to be the minimum exponent which guarantees no overflow/saturation will occur.**

    This is intended to minimize headroom, but in most situations cannot guarantee zero headroom in the result.

    ** Because of the symmetric saturation logic, a single LSb of saturation may occur, which is considered to be an
    acceptable trade-off for an extra bit of precision.

*/


/** Get the headroom of a 16-bit signed integer.
 *
 * \param a    Input value 
 * 
 * \returns    Number of leading sign bits of ``a``. 
 */
headroom_t xs3_headroom_s16(
    const int16_t a);

/** Get the headroom of a 32-bit signed integer.
 * 
 * \param a    Input value 
 * 
 * \returns    Number of leading sign bits of ``a``. 
 */
headroom_t xs3_headroom_s32(
    const int32_t a);

/** Get the headroom of a 64-bit signed integer.
 * 
 * IMPLEMENTATION_NOTE: Simplest way is just to check whether hte top word is 0 or -1,
 * then CLS whichever word that suggests needs it.
 * 
 * \param a    Input value 
 * 
 * \returns    Number of leading sign bits of ``a``. 
 */
headroom_t xs3_headroom_s64(
    const int64_t a);



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


/** Count the leading sign bits of a 16-bit BFP vector.
 * 
 * The number of leading sign bits of a sequence is the minimum of
 * that among its elements.
 *
 * \param a         input BFP vector
 * 
 * \returns    Number of leading sign bits of ``a``. 
 */
headroom_t bfp_headroom_vect_s16(
    bfp_s16_t* a);


/** Count the leading sign bits of a 32-bit BFP vector.
 * 
 * The number of leading sign bits of a sequence is the minimum of
 * that among its elements.
 *
 * \param a         input BFP vector
 * 
 * \returns    Number of leading sign bits of ``a``. 
 */
headroom_t bfp_headroom_vect_s32(
    bfp_s32_t* a);


/** Apply a left-shift to the elements of a 16-bit BFP vector.
 * 
 * Conceptually, the operation performed is:
 *      A[i] <- B[i] * 2^(shl)
 *        for each index i < length
 *        where A[] and B[] are integer vectors
 *              shift is an integer
 * 
 * This operation can be used to remove headroom from a BFP vector.
 * 
 * The operation saturates to 16-bit bounds.
 * 
 * \safe_in_place{a,b}
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector
 * \param shl   Number of bits to left shift
 */
void bfp_shl_vect_s16(
    bfp_s16_t* a,
    const bfp_s16_t* b,
    const left_shift_t shl);

/** Apply a left-shift to the elements of a 32-bit BFP vector.
 * 
 * Conceptually, the operation performed is:
 *      A[i] <- B[i] * 2^(shift)
 *        for each index i < length
 *        where A[] and B[] are integer vectors
 *              shift is an integer
 * 
 * This operation can be used to remove headroom from a BFP vector.
 * 
 * The operation saturates to 32-bit bounds.
 * 
 * \safe_in_place{a,b}
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector
 * \param shl   Number of bits to left shift
 */
void bfp_shl_vect_s32(
    bfp_s32_t* a,
    const bfp_s32_t* b,
    const left_shift_t shl);


/**
 * \brief Calculate the output exponent and input shifts needed to perform a
 * BFP vector addition.
 * 
 * This function is used to calcular the output exponent and operand shift parameter values 
 * required to compute the element-wise sum @f$\bar A@f$ of BFP vectors @f$\bar B@f$ and 
 * @f$\bar C@f$
 * 
 * Logic is identical for `bfp_s16_t`, `bfp_s32_t`, `bfp_complex_s16_t`, `bfp_complex_s32_t`.
 * 
 * The `allow_saturation` parameter is used to specify whether to check for the possibility of 
 * corner-case saturation. For an explanation of corner-case saturation, see \ref saturation.
 * Corner-case saturation is avoided if `allow_saturation` is zero.
 * 
 * For BFP vector addition, corner-case saturation may occur when BFP vectors `B` and `C` have
 * similar magnitudes (in the sense of their @f$L\infty@f$-norm). Specifically, for input BFP 
 * vectors @f$B@f$ and @f$C@f$, the corner case will be detected and avoided when 
 * @f$\left(B_{exp}-B_{hr}\right) = \left(C_{exp}-C_{hr}\right)@f$.
 * 
 * The outputs of this function `b_shr` and `c_shr` can be used with `xs3_add_vect_s16()`, 
 * `xs3_add_vect_s32()` as the shift values for that function's corresponding parameters. The output 
 * `a_exp` is the exponent associated with the result computed by those functions.
 * 
 * If a specific output exponent `desired_exp` is needed for the result, the `b_shr` and `c_shr` 
 * produced by this function can be adjusted according to the following:
 * \code{.c}
 *      exponent_t desired_exp = ...; // Value known a priori
 *      right_shift_t new_b_shr = b_shr + (desired_exp - a_exp);
 *      right_shift_t new_c_shr = c_shr + (desired_exp - a_exp);
 * \endcode
 * Note that using smaller values than necessary for `b_shr` and `c_shr` can result in saturation, 
 * and using larger values may result in unnecessary underflows.
 * 
 * \param[out] a_exp    Exponent @f$A_{exp}@f$ of the result vector @f$\bar A@f$
 * \param[out] b_shr    Right-shift to be applied to vector @f$\bar B@f$
 * \param[out] c_shr    Right-shift to be applied to vector @f$\bar C@f$
 * \param[in]  b_exp    Exponent @f$B_{exp}@f$ associated with @f$\bar B@f$
 * \param[in]  c_exp    Exponent @f$C_{exp}@f$ associated with @f$\bar C@f$
 * \param[in]  b_hr     Headroom @f$B_{hr}@f$ associated with @f$\bar B@f$
 * \param[in]  c_hr     Headroom @f$B_{hr}@f$ associated with @f$\bar B@f$
 * \param[in]  allow_saturation  Whether to avoid corner-case saturation.
 * 
 * \sa xs3_add_vect_s16
 * \sa xs3_add_vect_s32
 */
void bfp_add_vect_calc_params(
    exponent_t* a_exp,
    right_shift_t* b_shr,
    right_shift_t* c_shr,
    const exponent_t b_exp,
    const exponent_t c_exp,
    const headroom_t b_hr,
    const headroom_t c_hr,
    const unsigned allow_saturation);


/** Add two 16-bit BFP vectors together.
 * 
 * Conceptually, the operation performed is:
 *      A[] <- B[] + C[]
 *        where A[], B[] and C[] are BFP vectors
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * Any of ``a``, ``b`` and ``c`` may safely point at the same ``bfp_s16_t``.
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector 1
 * \param c     Input BFP vector 2
 */
void bfp_add_vect_s16(
    bfp_s16_t* a, 
    const bfp_s16_t* b, 
    const bfp_s16_t* c);


/** Add two 32-bit BFP vectors together.
 * 
 * Conceptually, the operation performed is:
 *      A[] <- B[] + C[]
 *        where A[], B[] and C[] are BFP vectors
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * Any of ``a``, ``b`` and ``c`` may safely point at the same ``bfp_s32_t``.
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector 1
 * \param c     Input BFP vector 2
 */
void bfp_add_vect_s32(
    bfp_s32_t* a, 
    const bfp_s32_t* b, 
    const bfp_s32_t* c);


/** Subtract one 16-bit BFP vector from another.
 * 
 * Conceptually, the operation performed is:
 *      A[] <- B[] - C[]
 *        where A[], B[] and C[] are BFP vectors
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * Any of ``a``, ``b`` and ``c`` may safely point at the same ``bfp_s16_t``.
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector 1
 * \param c     Input BFP vector 2
 */
void bfp_sub_vect_s16(
    bfp_s16_t* a, 
    const bfp_s16_t* b, 
    const bfp_s16_t* c);

/** Subtract one 32-bit BFP vector from another.
 * 
 * Conceptually, the operation performed is:
 *      A[] <- B[] - C[]
 *        where A[], B[] and C[] are BFP vectors
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * Any of ``a``, ``b`` and ``c`` may safely point at the same ``bfp_s32_t``.
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector 1
 * \param c     Input BFP vector 2
 */
void bfp_sub_vect_s32(
    bfp_s32_t* a, 
    const bfp_s32_t* b, 
    const bfp_s32_t* c);




/** Compute the element-by-element product two 16-bit BFP vectors.
 * 
 * Conceptually, the operation performed is:
 *      A[i] <- B[i] * C[i]
 *        for each index i in B[]
 *        where A[], B[] and C[] are BFP vectors
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * Any of ``a``, ``b`` and ``c`` may safely point at the same ``bfp_s16_t``.
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector 1
 * \param c     Input BFP vector 2
 */
void bfp_mul_vect_s16(
    bfp_s16_t* a, 
    const bfp_s16_t* b, 
    const bfp_s16_t* c);

/** Compute the element-by-element product of two 32-bit BFP vectors.
 * 
 * Conceptually, the operation performed is:
 *      A[i] <- B[i] * C[i]
 *        for each index i in B[]
 *        where A[], B[] and C[] are BFP vectors
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * Any of ``a``, ``b`` and ``c`` may safely point at the same ``bfp_s32_t``.
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector 1
 * \param c     Input BFP vector 2
 */
void bfp_mul_vect_s32(
    bfp_s32_t* a, 
    const bfp_s32_t* b, 
    const bfp_s32_t* c);


/** Multiply a 16-bit BFP vector by a scalar.
 * 
 * Conceptually, the operation performed is:
 *      A[] <- B[] * C
 *        where A[] and B[] are BFP vectors
 *              C is a scalar
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * It is safe to supply the same ``bfp_s16_t*`` for ``a`` and ``b``.
 * 
 * \param a         Output BFP vector
 * \param b         Input BFP vector
 * \param c         Scalar value
 */
void bfp_scalar_mul_vect_s16(
    bfp_s16_t* a, 
    const bfp_s16_t* b, 
    const float c);

/** Multiply a 32-bit BFP vector by a scalar.
 * 
 * Conceptually, the operation performed is:
 *      A[] <- B[] * C
 *        where A[] and B[] are BFP vectors
 *              C is a scalar
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * It is safe to supply the same ``bfp_s32_t*`` for ``a`` and ``b``.
 * 
 * \param a         Output BFP vector
 * \param b         Input BFP vector
 * \param c         Scalar Value
 */
void bfp_scalar_mul_vect_s32(
    bfp_s32_t* a, 
    const bfp_s32_t* b,
    const float c);

/** Take the absolute value of a 16-bit BFP vector, element-wise.
 * 
 * Conceptually, the operation performed is:
 *      A[i] <- abs( B[i] )
 *        for each index i in B[]
 *        where A[] and B[] are BFP vectors
 * 
 * This operation saturates to 16-bit bounds. (i.e. INT16_MIN becomes INT16_MAX)
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * It is safe to supply the same ``bfp_s16_t*`` for ``a`` and ``b``.
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector
 */
void bfp_abs_vect_s16(
    bfp_s16_t* a,
    const bfp_s16_t* b);

/** Take the absolute value of a 32-bit BFP vector, element-wise.
 * 
 * Conceptually, the operation performed is:
 *      A[i] <- abs( B[i] )
 *        for each index i in B[]
 *        where A[] and B[] are BFP vectors
 * 
 * This operation saturates to 32-bit bounds. (i.e. INT32_MIN becomes INT32_MAX)
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * It is safe to supply the same ``bfp_s32_t*`` for ``a`` and ``b``.
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector
 */
void bfp_abs_vect_s32(
    bfp_s32_t* a,
    const bfp_s32_t* b);

/** Compute the sum of the elements of a 16-bit BFP vector.
 * 
 * Conceptually, the operation performed is:
 *      A <- ( B[0] + B[1] + ... + B[N-1] }
 *        where A is a floating-point scalar
 *              B[] is a BFP vector
 *              N is the size of the vector B[]
 * 
 * \param b         Input BFP vector to be summed.
 * 
 * \returns         TODO
 */
float bfp_sum_s16(
    const bfp_s16_t* b);

/** Compute the sum of the elements of a 32-bit BFP vector.
 * 
 * Conceptually, the operation performed is:
 *      A <- ( B[0] + B[1] + ... + B[N-1] }
 *        where A is a floating-point scalar
 *              B[] is a BFP vector
 *              N is the size of the vector B[]
 * 
 * \param b         Input BFP vector to be summed.
 * 
 * \returns         TODO
 */
float bfp_sum_s32(
    const bfp_s32_t* b);


/** Compute the inner product of two 16-bit BFP vectors.
 * 
 * Conceptually, the operation performed is:
 *      A <- ( (B[0]*C[0]) + (B[1]*C[1]) + ... + (B[N-1]*C[N-1]) }
 *        where A is a floating-point scalar
 *              B[] is a BFP vector
 *              N is the size of the vector B[]
 * 
 * \param b         Input BFP vector 1
 * \param c         Input BFP vector 2
 */
float bfp_dot_s16(
    const bfp_s16_t* b, 
    const bfp_s16_t* c);

/** Compute the inner product of two 32-bit BFP vectors.
 * 
 * Conceptually, the operation performed is:
 *      A <- ( (B[0]*C[0]) + (B[1]*C[1]) + ... + (B[N-1]*C[N-1]) }
 *        where A is a floating-point scalar
 *              B[] is a BFP vector
 *              N is the size of the vector B[]
 * 
 * \param b         Input BFP vector 1
 * \param c         Input BFP vector 2
 */
float bfp_dot_s32(
    const bfp_s32_t* b, 
    const bfp_s32_t* c);


/** Clip the elements of a 16-bit BFP vector to a specified range.
 * 
 * Conceptually, the operation performed is:
 *      A[i] <-  min( U,  max( L,  B[i] ) )
 *        for each index i in B[]
 *        where A[] and B[] are BFP vectors
 *              U is a sclar; the range's upper bound
 *              L is a scalar; the range's lower bound
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * It is safe to supply the same ``bfp_s16_t*`` for ``a`` and ``b``.
 * 
 * \param a             Output BFP vector
 * \param b             Input BFP vector
 * \param lower_bound   Mantissa of the lower clipping bound
 * \param upper_bound   Mantissa of the upper clipping bound
 * \param bound_exp     Shared exponent of the clipping bounds
 */
void bfp_clip_vect_s16(
    bfp_s16_t* a, 
    const bfp_s16_t* b, 
    const int16_t lower_bound, 
    const int16_t upper_bound, 
    const int bound_exp);

/** Clip the elements of a 32-bit BFP vector to a specified range.
 * 
 * Conceptually, the operation performed is:
 *      A[i] <-  min( U,  max( L,  B[i] ) )
 *        for each index i in B[]
 *        where A[] and B[] are BFP vectors
 *              U is a sclar; the range's upper bound
 *              L is a scalar; the range's lower bound
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * It is safe to supply the same ``bfp_s32_t*`` for ``a`` and ``b``.
 * 
 * \param a             Output BFP vector
 * \param b             Input BFP vector
 * \param lower_bound   Mantissa of the lower clipping bound
 * \param upper_bound   Mantissa of the upper clipping bound
 * \param bound_exp     Shared exponent of the clipping bounds
 */
void bfp_clip_vect_s32(
    bfp_s32_t* a, 
    const bfp_s32_t* b, 
    const int32_t lower_bound, 
    const int32_t upper_bound, 
    const int bound_exp);


/** Rectify a 16-bit integer vector.
 * 
 * Conceptually, the operation performed is:
 *      A[i] <-  max( 0,  B[i] )
 *        for each index i < length
 *        where A[] and B[] are integer vectors
 * 
 * \param a     Output integer vector
 * \param b     Input integer vector
 */
void bfp_rect_vect_s16(
    bfp_s16_t* a,
    const bfp_s16_t* b);


/** Rectify the BFP vector ``b`` so that no elements are negative.
 * 
 *      a[i] = (b[i] < 0)? 0 : b[i];
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector
 */
void bfp_rect_vect_s32(
    bfp_s32_t* a, 
    const bfp_s32_t* b);

/** Convert a 32-bit BFP vector to a 16-bit BFP vector.
 * 
 * Conceptually, the operation performed is:
 *      A[] <-  B[]
 *        where A[] and B[] are BFP vectors
 * 
 * A conversion from 32-bit elements to 16-bit elements entails a possible loss 
 * of precision. Each element will be down-shifted an appropriate number of bits
 * (while taking into account headroom in ``b`` to minimize precision loss) and 
 * will be rounded to the nearest 16-bit value.
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector.
 **/
void bfp_s32_to_s16(
    bfp_s16_t* a,
    const bfp_s32_t* b);

/** Convert a 16-bit BFP vector to a 32-bit BFP vector.
 * 
 * Conceptually, the operation performed is:
 *      A[] <-  B[]
 *        where A[] and B[] are BFP vectors
 * 
 * A conversion from 16-bit elements to 32-bit elements entails no loss of precision. No
 * attempt will be made to minimize headroom in the output BFP vector (though the resulting
 * headroom will be recorded).
 * 
 * ``a->data`` must already be initialized to a valid memory buffer.
 * 
 * \param a     Output BFP vector
 * \param b     Input BFP vector.
 **/
void bfp_s16_to_s32(
    bfp_s32_t* a,
    const bfp_s16_t* b);


#ifdef __XC__
}   //extern "C"
#endif

#endif //XS3_BFP_H_