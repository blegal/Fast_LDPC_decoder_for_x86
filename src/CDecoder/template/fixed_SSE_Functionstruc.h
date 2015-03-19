/**
  Copyright (c) 2012-2015 "Bordeaux INP, Bertrand LE GAL"
  [http://legal.vvv.enseirb-matmeca.fr]

  This file is part of LDPC_C_Simulator.

  LDPC_C_Simulator is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __SSE_Functions__
#define __SSE_Functions__

#define llr_from_input(v)  ((2.0 * v)/(sigB * sigB))

#include <xmmintrin.h>

#define SSE_8S_LOAD(ptr)            (_mm_load_si128(ptr))
#define SSE_8S_UNCACHED_LOAD(ptr)   (_mm_stream_load_si128(ptr))
#define SSE_8S_STORE(ptr,v)         (_mm_store_si128(ptr,v))
#define SSE_8S_UNPACK_HIGH(a)       (_mm_srai_epi16(_mm_unpackhi_epi8(a,a),8))
#define SSE_8S_UNPACK_LOW(a)        (_mm_srai_epi16(_mm_unpacklo_epi8(a,a),8))
#define SSE_8S_PACK(hi,lo)          (_mm_packs_epi16(lo,hi))
#define SSE_8S_ADD(a,b)             (_mm_add_epi16(a,b))
#define SSE_8S_SUB(a,b)             (_mm_sub_epi16(a,b))
#define SSE_8S_ABS(a)               (_mm_abs_epi16(a))
#define SSE_8S_MAX(a,b)             (_mm_max_epi16(a,b))
#define SSE_8S_MIN(a,b)             (_mm_min_epi16(a,b))
#define SSE_8S_MUL(a,b)             (_mm_mullo_epi16(a,b))
#define SSE_8S_DIV(a,b)             (_mm_div_epi16(a,b))
#define SSE_8S_DIV32(a)             (_mm_srli_epi16(a,5))

const unsigned short sign_value = 0x8000;
const __m128i mask_sign = _mm_set_epi16(sign_value, sign_value, sign_value, sign_value, sign_value, sign_value, sign_value, sign_value);
inline __m128i SSE_8S_SIGN(__m128i a){
    __m128i b = _mm_and_si128 (a, mask_sign);            //  Invert all the bits
    __m128i f = _mm_xor_si128 (b, mask_sign);            //  Invert all the bits
    return f;
}

const unsigned short isign_value = 0xC000;
const __m128i mask_isign = _mm_set_epi16(isign_value, isign_value, isign_value, isign_value, isign_value, isign_value, isign_value, isign_value);
inline __m128i SSE_8S_invSIGN( __m128i val, __m128i sig ){
    __m128i x = _mm_xor_si128 (sig, mask_isign);      // PLUS LA VALEUR 64 CAR L'INSTRUCTION _mm_sign_epi16 A
    __m128i r = _mm_sign_epi16(val, x);     // PAS NULLE !!!
    return r;
}


//inline __m128i SSE_8S_MIN_1( __m128i a, __m128i min1){
//    return _mm_min_epi16(a, min1);
//}

#define SSE_8S_MIN_1(a,min1) \
    (SSE_8S_MIN(a,min1))

//inline __m128i SSE_8S_MIN_2( __m128i val, __m128i old_min1, __m128i min2){
//    return _mm_min_epi16(min2, _mm_max_epi16(val, old_min1));
//}

#define SSE_8S_MIN_2(val,old_min1,min2) \
    (SSE_8S_MIN(min2,SSE_8S_MAX(val,old_min1)))

/*
const __m128i min_var = _mm_set1_epi16( -127 );
const __m128i max_var = _mm_set1_epi16( +127 );
*/
#define SSE_8S_VAR_SATURATE(a) \
    (SSE_8S_MAX(SSE_8S_MIN(a, max_var), min_var))

#define SSE_8S_SATURATE(a, max, min) \
    (SSE_8S_MAX(SSE_8S_MIN(a, max), min))

//inline __m128i SSE_8S_VAR_SATURATE(__m128i a){
//    __m128i var_pos = _mm_set1_epi16(vSAT_POS_VAR);
//    __m128i var_neg = _mm_set1_epi16(vSAT_NEG_VAR);
//    __m128i b = _mm_min_epi16(a, var_pos);
//    __m128i c = _mm_max_epi16(b, var_neg);
//    return c;
//}
/*
inline __m128i SSE_8S_SUB_SATURATE_VAR(__m128i a, __m128i b){
    return SSE_8S_VAR_SATURATE(SSE_8S_SUB( a, b ));
}
*/
/*
inline __m128i SSE_8S_ADD_SATURATE_VAR(__m128i a, __m128i b){
    return SSE_8S_VAR_SATURATE(SSE_8S_ADD( a, b ));
}
*/
#define SSE_8S_SUB_AND_SATURATE_VAR(a,b,max,min) SSE_8S_SATURATE(SSE_8S_SUB(a,b),max,min)
#define SSE_8S_ADD_AND_SATURATE_VAR(a,b,max,min) SSE_8S_SATURATE(SSE_8S_ADD(a,b),max,min)

const unsigned short value = 0xFFFF;
const __m128i nn = _mm_set_epi16(value, value, value, value, value, value, value, value);
inline __m128i SSE_8S_CMOV( __m128i a, __m128i b, __m128i v1, __m128i v2){
    __m128i m1 = _mm_cmpeq_epi16 ( a, b );    // EQUALS
    __m128i m2 = _mm_xor_si128   ( m1, nn );  // NOT EQUALS
    __m128i m3 = _mm_and_si128   ( m1, v1 );
    __m128i m4 = _mm_and_si128   ( m2, v2 );
    return _mm_or_si128(m3, m4);
}

inline __m128i SSE_8S_XOR( __m128i a, __m128i b ){
    __m128i c = _mm_xor_si128(a, b); //  Invert all the bits
    return c;
}

inline int SSE_8S_XOR_REDUCE(__m128i reg){
    unsigned int *p;
    p = (unsigned int*)&reg;
    return ((p[0] | p[1] | p[2] | p[3]) != 0);
}

#endif
