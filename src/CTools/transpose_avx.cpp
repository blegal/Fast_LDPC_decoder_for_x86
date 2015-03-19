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

/* 
 * File:   Code_Analyzer.cpp
 * Author: cleroux
 * 
 * Created on September 17, 2013, 9:49 AM
 */

#include "transpose_avx.h"
#include "CTools.h"

#define LOAD_SIMD_FX    _mm256_load_si256
#define STORE_SIMD_FX   _mm256_store_si256

#define _mm256_unpacklo_epi128(a,b) (_mm256_permute2x128_si256(a,b,0x20))
#define _mm256_unpackhi_epi128(a,b) (_mm256_permute2x128_si256(a,b,0x31))

void uchar_transpose_avx(__m256i *src, __m256i *dst, int n)
{
    const int constN = n /32; // NOMBRE DE PAQUET (128 bits) PAR TRAME
    __m256i *p_input  = src;
    __m256i *p_output = dst;
    int loop = constN;

    while( loop-- ){
//IACA_START

        // STAGE 2 DU BUTTERFLY
        __m256i *copy_p_input  = p_input;

        // STAGE 2 DU BUTTERFLY
        __m256i a = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i b = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i c = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i d = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i e = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i f = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i g = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i h = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i i = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i j = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i k = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i l = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i m = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i n = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i o = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i p = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;

        __m256i a_0_15_b_0_15   = _mm256_unpacklo_epi8(a, b);
        __m256i a_16_31_b_16_31 = _mm256_unpackhi_epi8(a, b);
        __m256i c_0_15_d_0_15   = _mm256_unpacklo_epi8(c, d);
        __m256i c_16_31_d_16_31 = _mm256_unpackhi_epi8(c, d);
        __m256i e_0_15_f_0_15   = _mm256_unpacklo_epi8(e, f);
        __m256i e_16_31_f_16_31 = _mm256_unpackhi_epi8(e, f);
        __m256i g_0_15_h_0_15   = _mm256_unpacklo_epi8(g, h);
        __m256i g_16_31_h_16_31 = _mm256_unpackhi_epi8(g, h);
        __m256i i_0_15_j_0_15   = _mm256_unpacklo_epi8(i, j);
        __m256i i_16_31_j_16_31 = _mm256_unpackhi_epi8(i, j);
        __m256i k_0_15_l_0_15   = _mm256_unpacklo_epi8(k, l);
        __m256i k_16_31_l_16_31 = _mm256_unpackhi_epi8(k, l);
        __m256i m_0_15_n_0_15   = _mm256_unpacklo_epi8(m, n);
        __m256i m_16_31_n_16_31 = _mm256_unpackhi_epi8(m, n);
        __m256i o_0_15_p_0_15   = _mm256_unpacklo_epi8(o, p);
        __m256i o_16_31_p_16_31 = _mm256_unpackhi_epi8(o, p);

        __m256i a_0_7_d_0_7     = _mm256_unpacklo_epi16(a_0_15_b_0_15,   c_0_15_d_0_15);
        __m256i a_8_15_d_8_15   = _mm256_unpackhi_epi16(a_0_15_b_0_15,   c_0_15_d_0_15);
        __m256i a_16_23_d_16_23 = _mm256_unpacklo_epi16(a_16_31_b_16_31, c_16_31_d_16_31);
        __m256i a_24_31_d_24_31 = _mm256_unpackhi_epi16(a_16_31_b_16_31, c_16_31_d_16_31);
        __m256i e_0_7_h_0_7     = _mm256_unpacklo_epi16(e_0_15_f_0_15,   g_0_15_h_0_15);
        __m256i e_8_15_h_8_15   = _mm256_unpackhi_epi16(e_0_15_f_0_15,   g_0_15_h_0_15);
        __m256i e_16_23_h_16_23 = _mm256_unpacklo_epi16(e_16_31_f_16_31, g_16_31_h_16_31);
        __m256i e_24_31_h_24_31 = _mm256_unpackhi_epi16(e_16_31_f_16_31, g_16_31_h_16_31);
        __m256i i_0_7_l_0_7     = _mm256_unpacklo_epi16(i_0_15_j_0_15,   k_0_15_l_0_15);
        __m256i i_8_15_l_8_15   = _mm256_unpackhi_epi16(i_0_15_j_0_15,   k_0_15_l_0_15);
        __m256i i_16_23_l_16_23 = _mm256_unpacklo_epi16(i_16_31_j_16_31, k_16_31_l_16_31);
        __m256i i_24_31_l_24_31 = _mm256_unpackhi_epi16(i_16_31_j_16_31, k_16_31_l_16_31);
        __m256i m_0_7_p_0_7     = _mm256_unpacklo_epi16(m_0_15_n_0_15,   o_0_15_p_0_15);
        __m256i m_8_15_p_8_15   = _mm256_unpackhi_epi16(m_0_15_n_0_15,   o_0_15_p_0_15);
        __m256i m_16_23_p_16_23 = _mm256_unpacklo_epi16(m_16_31_n_16_31, o_16_31_p_16_31);
        __m256i m_24_31_p_24_31 = _mm256_unpackhi_epi16(m_16_31_n_16_31, o_16_31_p_16_31);

        //////////////////////////////////////////////////////////////////////////////////////////////////
        
        __m256i A               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i B               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i C               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i D               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i E               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i F               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i G               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i H               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i I               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i J               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i K               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i L               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i M               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i N               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i O               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        __m256i P               = LOAD_SIMD_FX(copy_p_input); copy_p_input += constN;
        
        __m256i A_0_15_B_0_15   = _mm256_unpacklo_epi8(A, B);
        __m256i A_16_31_B_16_31 = _mm256_unpackhi_epi8(A, B);
        __m256i C_0_15_D_0_15   = _mm256_unpacklo_epi8(C, D);
        __m256i C_16_31_D_16_31 = _mm256_unpackhi_epi8(C, D);
        __m256i E_0_15_F_0_15   = _mm256_unpacklo_epi8(E, F);
        __m256i E_16_31_F_16_31 = _mm256_unpackhi_epi8(E, F);
        __m256i G_0_15_H_0_15   = _mm256_unpacklo_epi8(G, H);
        __m256i G_16_31_H_16_31 = _mm256_unpackhi_epi8(G, H);
        __m256i I_0_15_J_0_15   = _mm256_unpacklo_epi8(I, J);
        __m256i I_16_31_J_16_31 = _mm256_unpackhi_epi8(I, J);
        __m256i K_0_15_L_0_15   = _mm256_unpacklo_epi8(K, L);
        __m256i K_16_31_L_16_31 = _mm256_unpackhi_epi8(K, L);
        __m256i M_0_15_N_0_15   = _mm256_unpacklo_epi8(M, N);
        __m256i M_16_31_N_16_31 = _mm256_unpackhi_epi8(M, N);
        __m256i O_0_15_P_0_15   = _mm256_unpacklo_epi8(O, P);
        __m256i O_16_31_P_16_31 = _mm256_unpackhi_epi8(O, P);

        __m256i A_0_7_D_0_7     = _mm256_unpacklo_epi16(A_0_15_B_0_15,   C_0_15_D_0_15);
        __m256i E_0_7_H_0_7     = _mm256_unpacklo_epi16(E_0_15_F_0_15,   G_0_15_H_0_15);
        __m256i I_0_7_L_0_7     = _mm256_unpacklo_epi16(I_0_15_J_0_15,   K_0_15_L_0_15);
        __m256i M_0_7_P_0_7     = _mm256_unpacklo_epi16(M_0_15_N_0_15,   O_0_15_P_0_15);

        __m256i A_8_15_D_8_15   = _mm256_unpackhi_epi16(A_0_15_B_0_15,   C_0_15_D_0_15);
        __m256i E_8_15_H_8_15   = _mm256_unpackhi_epi16(E_0_15_F_0_15,   G_0_15_H_0_15);
        __m256i I_8_15_L_8_15   = _mm256_unpackhi_epi16(I_0_15_J_0_15,   K_0_15_L_0_15);
        __m256i M_8_15_P_8_15   = _mm256_unpackhi_epi16(M_0_15_N_0_15,   O_0_15_P_0_15);

        __m256i A_16_23_D_16_23 = _mm256_unpacklo_epi16(A_16_31_B_16_31, C_16_31_D_16_31);
        __m256i E_16_23_H_16_23 = _mm256_unpacklo_epi16(E_16_31_F_16_31, G_16_31_H_16_31);
        __m256i I_16_23_L_16_23 = _mm256_unpacklo_epi16(I_16_31_J_16_31, K_16_31_L_16_31);
        __m256i M_16_23_P_16_23 = _mm256_unpacklo_epi16(M_16_31_N_16_31, O_16_31_P_16_31);

        __m256i A_24_31_D_24_31 = _mm256_unpackhi_epi16(A_16_31_B_16_31, C_16_31_D_16_31);
        __m256i E_24_31_H_24_31 = _mm256_unpackhi_epi16(E_16_31_F_16_31, G_16_31_H_16_31);
        __m256i I_24_31_L_24_31 = _mm256_unpackhi_epi16(I_16_31_J_16_31, K_16_31_L_16_31);
        __m256i M_24_31_P_24_31 = _mm256_unpackhi_epi16(M_16_31_N_16_31, O_16_31_P_16_31);

        //
        // STAGE 4 DU BUTTERFLY
        //
        __m256i a_0_3_h_0_3      = _mm256_unpacklo_epi32(a_0_7_d_0_7,     e_0_7_h_0_7); // LOW
        __m256i a_4_7_h_4_7      = _mm256_unpackhi_epi32(a_0_7_d_0_7,     e_0_7_h_0_7); // HIGH
        __m256i a_8_11_h_8_11    = _mm256_unpacklo_epi32(a_8_15_d_8_15,   e_8_15_h_8_15); // LOW
        __m256i a_12_15_h_12_15  = _mm256_unpackhi_epi32(a_8_15_d_8_15,   e_8_15_h_8_15); // HIGH
        __m256i a_16_19_h_16_19  = _mm256_unpacklo_epi32(a_16_23_d_16_23, e_16_23_h_16_23); // LOW
        __m256i a_20_23_h_20_23  = _mm256_unpackhi_epi32(a_16_23_d_16_23, e_16_23_h_16_23); // HIGH
        __m256i a_24_27_h_24_27  = _mm256_unpacklo_epi32(a_24_31_d_24_31, e_24_31_h_24_31); // LOW
        __m256i a_28_31_h_28_31  = _mm256_unpackhi_epi32(a_24_31_d_24_31, e_24_31_h_24_31); // HIGH
        __m256i i_0_3_p_0_3      = _mm256_unpacklo_epi32(i_0_7_l_0_7,     m_0_7_p_0_7); // LOW
        __m256i i_4_7_p_4_7      = _mm256_unpackhi_epi32(i_0_7_l_0_7,     m_0_7_p_0_7); // HIGH
        __m256i i_8_11_p_8_11    = _mm256_unpacklo_epi32(i_8_15_l_8_15,   m_8_15_p_8_15); // LOW
        __m256i i_12_15_p_12_15  = _mm256_unpackhi_epi32(i_8_15_l_8_15,   m_8_15_p_8_15); // HIGH
        __m256i i_16_19_p_16_19  = _mm256_unpacklo_epi32(i_16_23_l_16_23, m_16_23_p_16_23); // LOW
        __m256i i_20_23_p_20_23  = _mm256_unpackhi_epi32(i_16_23_l_16_23, m_16_23_p_16_23); // HIGH
        __m256i i_24_27_p_24_27  = _mm256_unpacklo_epi32(i_24_31_l_24_31, m_24_31_p_24_31); // LOW
        __m256i i_28_31_p_28_31  = _mm256_unpackhi_epi32(i_24_31_l_24_31, m_24_31_p_24_31); // HIGH

        __m256i A_0_3_H_0_3      = _mm256_unpacklo_epi32(A_0_7_D_0_7,     E_0_7_H_0_7);     // LOW
        __m256i A_4_7_H_4_7      = _mm256_unpackhi_epi32(A_0_7_D_0_7,     E_0_7_H_0_7);     // HIGH
        __m256i A_8_11_H_8_11    = _mm256_unpacklo_epi32(A_8_15_D_8_15,   E_8_15_H_8_15);     // LOW
        __m256i A_12_15_H_12_15  = _mm256_unpackhi_epi32(A_8_15_D_8_15,   E_8_15_H_8_15);     // HIGH
        __m256i A_16_19_H_16_19  = _mm256_unpacklo_epi32(A_16_23_D_16_23, E_16_23_H_16_23); // LOW
        __m256i A_20_23_H_20_23  = _mm256_unpackhi_epi32(A_16_23_D_16_23, E_16_23_H_16_23); // HIGH
        __m256i A_24_27_H_24_27  = _mm256_unpacklo_epi32(A_24_31_D_24_31, E_24_31_H_24_31); // LOW
        __m256i A_28_31_H_28_31  = _mm256_unpackhi_epi32(A_24_31_D_24_31, E_24_31_H_24_31); // HIGH
        __m256i I_0_3_P_0_3      = _mm256_unpacklo_epi32(I_0_7_L_0_7,     M_0_7_P_0_7);     // LOW
        __m256i I_4_7_P_4_7      = _mm256_unpackhi_epi32(I_0_7_L_0_7,     M_0_7_P_0_7);     // HIGH
        __m256i I_8_11_P_8_11    = _mm256_unpacklo_epi32(I_8_15_L_8_15,   M_8_15_P_8_15);     // LOW
        __m256i I_12_15_P_12_15  = _mm256_unpackhi_epi32(I_8_15_L_8_15,   M_8_15_P_8_15);     // HIGH
        __m256i I_16_19_P_16_19  = _mm256_unpacklo_epi32(I_16_23_L_16_23, M_16_23_P_16_23); // LOW
        __m256i I_20_23_P_20_23  = _mm256_unpackhi_epi32(I_16_23_L_16_23, M_16_23_P_16_23); // HIGH
        __m256i I_24_27_P_24_27  = _mm256_unpacklo_epi32(I_24_31_L_24_31, M_24_31_P_24_31); // LOW
        __m256i I_28_31_P_28_31  = _mm256_unpackhi_epi32(I_24_31_L_24_31, M_24_31_P_24_31); // HIGH

        //
        // STAGE 5 DU BUTTERFLY
        //
        __m256i a_0_1_p_0_1      = _mm256_unpacklo_epi64(a_0_3_h_0_3,     i_0_3_p_0_3);
        __m256i a_2_3_p_2_3      = _mm256_unpackhi_epi64(a_0_3_h_0_3,     i_0_3_p_0_3);
        __m256i a_4_5_p_4_5      = _mm256_unpacklo_epi64(a_4_7_h_4_7,     i_4_7_p_4_7);
        __m256i a_6_7_p_6_7      = _mm256_unpackhi_epi64(a_4_7_h_4_7,     i_4_7_p_4_7);
        __m256i a_8_9_p_8_9      = _mm256_unpacklo_epi64(a_8_11_h_8_11,   i_8_11_p_8_11);
        __m256i a_10_11_p_10_11  = _mm256_unpackhi_epi64(a_8_11_h_8_11,   i_8_11_p_8_11);
        __m256i a_12_13_p_12_13  = _mm256_unpacklo_epi64(a_12_15_h_12_15, i_12_15_p_12_15);
        __m256i a_14_15_p_14_15  = _mm256_unpackhi_epi64(a_12_15_h_12_15, i_12_15_p_12_15);
        __m256i a_16_17_p_16_17  = _mm256_unpacklo_epi64(a_16_19_h_16_19, i_16_19_p_16_19);
        __m256i a_18_19_p_18_19  = _mm256_unpackhi_epi64(a_16_19_h_16_19, i_16_19_p_16_19);
        __m256i a_20_21_p_20_21  = _mm256_unpacklo_epi64(a_20_23_h_20_23, i_20_23_p_20_23);
        __m256i a_22_23_p_22_23  = _mm256_unpackhi_epi64(a_20_23_h_20_23, i_20_23_p_20_23);
        __m256i a_24_25_p_24_25  = _mm256_unpacklo_epi64(a_24_27_h_24_27, i_24_27_p_24_27);
        __m256i a_26_27_p_26_27  = _mm256_unpackhi_epi64(a_24_27_h_24_27, i_24_27_p_24_27);
        __m256i a_28_29_p_28_29  = _mm256_unpacklo_epi64(a_28_31_h_28_31, i_28_31_p_28_31);
        __m256i a_30_31_p_30_31  = _mm256_unpackhi_epi64(a_28_31_h_28_31, i_28_31_p_28_31);

        __m256i A_0_1_P_0_1      = _mm256_unpacklo_epi64(A_0_3_H_0_3,     I_0_3_P_0_3);
        __m256i A_2_3_P_2_3      = _mm256_unpackhi_epi64(A_0_3_H_0_3,     I_0_3_P_0_3);
        __m256i A_4_5_P_4_5      = _mm256_unpacklo_epi64(A_4_7_H_4_7,     I_4_7_P_4_7);
        __m256i A_6_7_P_6_7      = _mm256_unpackhi_epi64(A_4_7_H_4_7,     I_4_7_P_4_7);
        __m256i A_8_9_P_8_9      = _mm256_unpacklo_epi64(A_8_11_H_8_11,   I_8_11_P_8_11);
        __m256i A_10_11_P_10_11  = _mm256_unpackhi_epi64(A_8_11_H_8_11,   I_8_11_P_8_11);
        __m256i A_12_13_P_12_13  = _mm256_unpacklo_epi64(A_12_15_H_12_15, I_12_15_P_12_15);
        __m256i A_14_15_P_14_15  = _mm256_unpackhi_epi64(A_12_15_H_12_15, I_12_15_P_12_15);
        __m256i A_16_17_P_16_17  = _mm256_unpacklo_epi64(A_16_19_H_16_19, I_16_19_P_16_19);
        __m256i A_18_19_P_18_19  = _mm256_unpackhi_epi64(A_16_19_H_16_19, I_16_19_P_16_19);
        __m256i A_20_21_P_20_21  = _mm256_unpacklo_epi64(A_20_23_H_20_23, I_20_23_P_20_23);
        __m256i A_22_23_P_22_23  = _mm256_unpackhi_epi64(A_20_23_H_20_23, I_20_23_P_20_23);
        __m256i A_24_25_P_24_25  = _mm256_unpacklo_epi64(A_24_27_H_24_27, I_24_27_P_24_27);
        __m256i A_26_27_P_26_27  = _mm256_unpackhi_epi64(A_24_27_H_24_27, I_24_27_P_24_27);
        __m256i A_28_29_P_28_29  = _mm256_unpacklo_epi64(A_28_31_H_28_31, I_28_31_P_28_31);
        __m256i A_30_31_P_30_31  = _mm256_unpackhi_epi64(A_28_31_H_28_31, I_28_31_P_28_31);
        
        
/*
 SELECT4(src1, src2, control){
    CASE(control[1:0])
        0:    tmp[127:0] := src1[127:0]
        1:    tmp[127:0] := src1[255:128]
        2:    tmp[127:0] := src2[127:0]
        3:    tmp[127:0] := src2[255:128]
    ESAC
    IF control[3]
        tmp[127:0] := 0
    ENDIF
    RETURN tmp[127:0]
}

 * __m256i _mm256_permute2x128_si256 (__m256i a, __m256i b, const int imm)
  dst[127:  0] := SELECT4(a[255:0], b[255:0], imm[3:0])
  dst[255:128] := SELECT4(b[255:0], b[255:0], imm[7:4])
  dst[MAX:256] := 0
 */        
        //
        // STAGE 6 DU BUTTERFLY
        //
        __m256i a_0_p_0_A_0_P_0     = _mm256_unpacklo_epi128(a_0_1_p_0_1,     A_0_1_P_0_1);
        __m256i a_1_p_1_A_1_P_1     = _mm256_unpackhi_epi128(a_0_1_p_0_1,     A_0_1_P_0_1);

        __m256i a_2_p_2_A_2_P_2     = _mm256_unpacklo_epi128(a_2_3_p_2_3,     A_2_3_P_2_3);
        __m256i a_3_p_3_A_3_P_3     = _mm256_unpackhi_epi128(a_2_3_p_2_3,     A_2_3_P_2_3);

        __m256i a_4_p_4_A_4_P_4     = _mm256_unpacklo_epi128(a_4_5_p_4_5,     A_4_5_P_4_5);
        __m256i a_5_p_5_A_5_P_5     = _mm256_unpackhi_epi128(a_4_5_p_4_5,     A_4_5_P_4_5);

        __m256i a_6_p_6_A_6_P_6     = _mm256_unpacklo_epi128(a_6_7_p_6_7,     A_6_7_P_6_7);
        __m256i a_7_p_7_A_7_P_7     = _mm256_unpackhi_epi128(a_6_7_p_6_7,     A_6_7_P_6_7);

        __m256i a_8_p_8_A_8_P_8     = _mm256_unpacklo_epi128(a_8_9_p_8_9,     A_8_9_P_8_9);
        __m256i a_9_p_9_A_9_P_9     = _mm256_unpackhi_epi128(a_8_9_p_8_9,     A_8_9_P_8_9);

        __m256i a_10_p_10_A_10_P_12 = _mm256_unpacklo_epi128(a_10_11_p_10_11, A_10_11_P_10_11);
        __m256i a_11_p_11_A_11_P_12 = _mm256_unpackhi_epi128(a_10_11_p_10_11, A_10_11_P_10_11);

        __m256i a_12_p_12_A_12_P_12 = _mm256_unpacklo_epi128(a_12_13_p_12_13, A_12_13_P_12_13);
        __m256i a_13_p_13_A_13_P_13 = _mm256_unpackhi_epi128(a_12_13_p_12_13, A_12_13_P_12_13);

        __m256i a_14_p_14_A_14_P_14 = _mm256_unpacklo_epi128(a_14_15_p_14_15, A_14_15_P_14_15);
        __m256i a_15_p_15_A_15_P_15 = _mm256_unpackhi_epi128(a_14_15_p_14_15, A_14_15_P_14_15);

        __m256i a_16_p_16_A_16_P_16 = _mm256_unpacklo_epi128(a_16_17_p_16_17, A_16_17_P_16_17);
        __m256i a_17_p_17_A_17_P_17 = _mm256_unpackhi_epi128(a_16_17_p_16_17, A_16_17_P_16_17);

        __m256i a_18_p_18_A_18_P_18 = _mm256_unpacklo_epi128(a_18_19_p_18_19, A_18_19_P_18_19);
        __m256i a_19_p_19_A_19_P_19 = _mm256_unpackhi_epi128(a_18_19_p_18_19, A_18_19_P_18_19);

        __m256i a_20_p_20_A_20_P_20 = _mm256_unpacklo_epi128(a_20_21_p_20_21, A_20_21_P_20_21);
        __m256i a_21_p_21_A_21_P_21 = _mm256_unpackhi_epi128(a_20_21_p_20_21, A_20_21_P_20_21);

        __m256i a_22_p_22_A_22_P_22 = _mm256_unpacklo_epi128(a_22_23_p_22_23, A_22_23_P_22_23);
        __m256i a_23_p_23_A_23_P_23 = _mm256_unpackhi_epi128(a_22_23_p_22_23, A_22_23_P_22_23);

        __m256i a_24_p_24_A_24_P_24 = _mm256_unpacklo_epi128(a_24_25_p_24_25, A_24_25_P_24_25);
        __m256i a_25_p_25_A_25_P_25 = _mm256_unpackhi_epi128(a_24_25_p_24_25, A_24_25_P_24_25);

        __m256i a_26_p_26_A_26_P_26 = _mm256_unpacklo_epi128(a_26_27_p_26_27, A_26_27_P_26_27);
        __m256i a_27_p_27_A_27_P_27 = _mm256_unpackhi_epi128(a_26_27_p_26_27, A_26_27_P_26_27);

        __m256i a_28_p_28_A_28_P_28 = _mm256_unpacklo_epi128(a_28_29_p_28_29, A_28_29_P_28_29);
        __m256i a_29_p_29_A_29_P_29 = _mm256_unpackhi_epi128(a_28_29_p_28_29, A_28_29_P_28_29);

        __m256i a_30_p_30_A_30_P_30 = _mm256_unpacklo_epi128(a_30_31_p_30_31, A_30_31_P_30_31);
        __m256i a_31_p_31_A_31_P_31 = _mm256_unpackhi_epi128(a_30_31_p_30_31, A_30_31_P_30_31);

        STORE_SIMD_FX(p_output++, a_0_p_0_A_0_P_0);
        STORE_SIMD_FX(p_output++, a_2_p_2_A_2_P_2);
        STORE_SIMD_FX(p_output++, a_4_p_4_A_4_P_4);
        STORE_SIMD_FX(p_output++, a_6_p_6_A_6_P_6);
        STORE_SIMD_FX(p_output++, a_8_p_8_A_8_P_8);
        STORE_SIMD_FX(p_output++, a_10_p_10_A_10_P_12);
        STORE_SIMD_FX(p_output++, a_12_p_12_A_12_P_12);
        STORE_SIMD_FX(p_output++, a_14_p_14_A_14_P_14);
        STORE_SIMD_FX(p_output++, a_16_p_16_A_16_P_16);
        STORE_SIMD_FX(p_output++, a_18_p_18_A_18_P_18);
        STORE_SIMD_FX(p_output++, a_20_p_20_A_20_P_20);
        STORE_SIMD_FX(p_output++, a_22_p_22_A_22_P_22);
        STORE_SIMD_FX(p_output++, a_24_p_24_A_24_P_24);
        STORE_SIMD_FX(p_output++, a_26_p_26_A_26_P_26);
        STORE_SIMD_FX(p_output++, a_28_p_28_A_28_P_28);
        STORE_SIMD_FX(p_output++, a_30_p_30_A_30_P_30);
        
        STORE_SIMD_FX(p_output++, a_1_p_1_A_1_P_1);
        STORE_SIMD_FX(p_output++, a_3_p_3_A_3_P_3);
        STORE_SIMD_FX(p_output++, a_5_p_5_A_5_P_5);
        STORE_SIMD_FX(p_output++, a_7_p_7_A_7_P_7);
        STORE_SIMD_FX(p_output++, a_9_p_9_A_9_P_9);
        STORE_SIMD_FX(p_output++, a_11_p_11_A_11_P_12);
        STORE_SIMD_FX(p_output++, a_13_p_13_A_13_P_13);
        STORE_SIMD_FX(p_output++, a_15_p_15_A_15_P_15);
        STORE_SIMD_FX(p_output++, a_17_p_17_A_17_P_17);
        STORE_SIMD_FX(p_output++, a_19_p_19_A_19_P_19);
        STORE_SIMD_FX(p_output++, a_21_p_21_A_21_P_21);
        STORE_SIMD_FX(p_output++, a_23_p_23_A_23_P_23);
        STORE_SIMD_FX(p_output++, a_25_p_25_A_25_P_25);
        STORE_SIMD_FX(p_output++, a_27_p_27_A_27_P_27);
        STORE_SIMD_FX(p_output++, a_29_p_29_A_29_P_29);
        STORE_SIMD_FX(p_output++, a_31_p_31_A_31_P_31);
        p_input  += 1;
    }    
//IACA_END
}

#define LOAD_AND_DECIDE(a)    _mm256_and_si256(_mm256_cmpgt_epi8(_mm256_load_si256(a),zero),mask)

void uchar_itranspose_avx(__m256i *src, __m256i *dst, int n)
{
    const int constN = n /32; // NOMBRE DE PAQUET (128 bits) PAR TRAME
    __m256i *p_input  = src;
    __m256i *p_output = dst;
    int loop = constN;
    const __m256i mask = _mm256_set1_epi8( 0x01 );
    const __m256i zero = _mm256_set1_epi8( 0x00 );

    while( loop-- ){
//IACA_START
        __m256i a = LOAD_AND_DECIDE(p_input +  0);
        __m256i b = LOAD_AND_DECIDE(p_input +  1);
        __m256i c = LOAD_AND_DECIDE(p_input +  2);
        __m256i d = LOAD_AND_DECIDE(p_input +  3);
        __m256i e = LOAD_AND_DECIDE(p_input +  4);
        __m256i f = LOAD_AND_DECIDE(p_input +  5);
        __m256i g = LOAD_AND_DECIDE(p_input +  6);
        __m256i h = LOAD_AND_DECIDE(p_input +  7);
        __m256i i = LOAD_AND_DECIDE(p_input +  8);
        __m256i j = LOAD_AND_DECIDE(p_input +  9);
        __m256i k = LOAD_AND_DECIDE(p_input + 10);
        __m256i l = LOAD_AND_DECIDE(p_input + 11);
        __m256i m = LOAD_AND_DECIDE(p_input + 12);
        __m256i n = LOAD_AND_DECIDE(p_input + 13);
        __m256i o = LOAD_AND_DECIDE(p_input + 14);
        __m256i p = LOAD_AND_DECIDE(p_input + 15);

        __m256i a_0_15_b_0_15   = _mm256_unpacklo_epi8(a, b);
        __m256i a_16_31_b_16_31 = _mm256_unpackhi_epi8(a, b);
        __m256i c_0_15_d_0_15   = _mm256_unpacklo_epi8(c, d);
        __m256i c_16_31_d_16_31 = _mm256_unpackhi_epi8(c, d);
        __m256i e_0_15_f_0_15   = _mm256_unpacklo_epi8(e, f);
        __m256i e_16_31_f_16_31 = _mm256_unpackhi_epi8(e, f);
        __m256i g_0_15_h_0_15   = _mm256_unpacklo_epi8(g, h);
        __m256i g_16_31_h_16_31 = _mm256_unpackhi_epi8(g, h);
        __m256i i_0_15_j_0_15   = _mm256_unpacklo_epi8(i, j);
        __m256i i_16_31_j_16_31 = _mm256_unpackhi_epi8(i, j);
        __m256i k_0_15_l_0_15   = _mm256_unpacklo_epi8(k, l);
        __m256i k_16_31_l_16_31 = _mm256_unpackhi_epi8(k, l);
        __m256i m_0_15_n_0_15   = _mm256_unpacklo_epi8(m, n);
        __m256i m_16_31_n_16_31 = _mm256_unpackhi_epi8(m, n);
        __m256i o_0_15_p_0_15   = _mm256_unpacklo_epi8(o, p);
        __m256i o_16_31_p_16_31 = _mm256_unpackhi_epi8(o, p);

        __m256i a_0_7_d_0_7     = _mm256_unpacklo_epi16(a_0_15_b_0_15,   c_0_15_d_0_15);
        __m256i a_8_15_d_8_15   = _mm256_unpackhi_epi16(a_0_15_b_0_15,   c_0_15_d_0_15);
        __m256i a_16_23_d_16_23 = _mm256_unpacklo_epi16(a_16_31_b_16_31, c_16_31_d_16_31);
        __m256i a_24_31_d_24_31 = _mm256_unpackhi_epi16(a_16_31_b_16_31, c_16_31_d_16_31);
        __m256i e_0_7_h_0_7     = _mm256_unpacklo_epi16(e_0_15_f_0_15,   g_0_15_h_0_15);
        __m256i e_8_15_h_8_15   = _mm256_unpackhi_epi16(e_0_15_f_0_15,   g_0_15_h_0_15);
        __m256i e_16_23_h_16_23 = _mm256_unpacklo_epi16(e_16_31_f_16_31, g_16_31_h_16_31);
        __m256i e_24_31_h_24_31 = _mm256_unpackhi_epi16(e_16_31_f_16_31, g_16_31_h_16_31);
        __m256i i_0_7_l_0_7     = _mm256_unpacklo_epi16(i_0_15_j_0_15,   k_0_15_l_0_15);
        __m256i i_8_15_l_8_15   = _mm256_unpackhi_epi16(i_0_15_j_0_15,   k_0_15_l_0_15);
        __m256i i_16_23_l_16_23 = _mm256_unpacklo_epi16(i_16_31_j_16_31, k_16_31_l_16_31);
        __m256i i_24_31_l_24_31 = _mm256_unpackhi_epi16(i_16_31_j_16_31, k_16_31_l_16_31);
        __m256i m_0_7_p_0_7     = _mm256_unpacklo_epi16(m_0_15_n_0_15,   o_0_15_p_0_15);
        __m256i m_8_15_p_8_15   = _mm256_unpackhi_epi16(m_0_15_n_0_15,   o_0_15_p_0_15);
        __m256i m_16_23_p_16_23 = _mm256_unpacklo_epi16(m_16_31_n_16_31, o_16_31_p_16_31);
        __m256i m_24_31_p_24_31 = _mm256_unpackhi_epi16(m_16_31_n_16_31, o_16_31_p_16_31);

        //////////////////////////////////////////////////////////////////////////////////////////////////
        
        __m256i A               = LOAD_AND_DECIDE(p_input + 16);
        __m256i B               = LOAD_AND_DECIDE(p_input + 17);
        __m256i C               = LOAD_AND_DECIDE(p_input + 18);
        __m256i D               = LOAD_AND_DECIDE(p_input + 19);
        __m256i E               = LOAD_AND_DECIDE(p_input + 20);
        __m256i F               = LOAD_AND_DECIDE(p_input + 21);
        __m256i G               = LOAD_AND_DECIDE(p_input + 22);
        __m256i H               = LOAD_AND_DECIDE(p_input + 23);
        __m256i I               = LOAD_AND_DECIDE(p_input + 24);
        __m256i J               = LOAD_AND_DECIDE(p_input + 25);
        __m256i K               = LOAD_AND_DECIDE(p_input + 26);
        __m256i L               = LOAD_AND_DECIDE(p_input + 27);
        __m256i M               = LOAD_AND_DECIDE(p_input + 28);
        __m256i N               = LOAD_AND_DECIDE(p_input + 29);
        __m256i O               = LOAD_AND_DECIDE(p_input + 30);
        __m256i P               = LOAD_AND_DECIDE(p_input + 31);
        
        __m256i A_0_15_B_0_15   = _mm256_unpacklo_epi8(A, B);
        __m256i A_16_31_B_16_31 = _mm256_unpackhi_epi8(A, B);
        __m256i C_0_15_D_0_15   = _mm256_unpacklo_epi8(C, D);
        __m256i C_16_31_D_16_31 = _mm256_unpackhi_epi8(C, D);
        __m256i E_0_15_F_0_15   = _mm256_unpacklo_epi8(E, F);
        __m256i E_16_31_F_16_31 = _mm256_unpackhi_epi8(E, F);
        __m256i G_0_15_H_0_15   = _mm256_unpacklo_epi8(G, H);
        __m256i G_16_31_H_16_31 = _mm256_unpackhi_epi8(G, H);
        __m256i I_0_15_J_0_15   = _mm256_unpacklo_epi8(I, J);
        __m256i I_16_31_J_16_31 = _mm256_unpackhi_epi8(I, J);
        __m256i K_0_15_L_0_15   = _mm256_unpacklo_epi8(K, L);
        __m256i K_16_31_L_16_31 = _mm256_unpackhi_epi8(K, L);
        __m256i M_0_15_N_0_15   = _mm256_unpacklo_epi8(M, N);
        __m256i M_16_31_N_16_31 = _mm256_unpackhi_epi8(M, N);
        __m256i O_0_15_P_0_15   = _mm256_unpacklo_epi8(O, P);
        __m256i O_16_31_P_16_31 = _mm256_unpackhi_epi8(O, P);

        __m256i A_0_7_D_0_7     = _mm256_unpacklo_epi16(A_0_15_B_0_15,   C_0_15_D_0_15);
        __m256i E_0_7_H_0_7     = _mm256_unpacklo_epi16(E_0_15_F_0_15,   G_0_15_H_0_15);
        __m256i I_0_7_L_0_7     = _mm256_unpacklo_epi16(I_0_15_J_0_15,   K_0_15_L_0_15);
        __m256i M_0_7_P_0_7     = _mm256_unpacklo_epi16(M_0_15_N_0_15,   O_0_15_P_0_15);

        __m256i A_8_15_D_8_15   = _mm256_unpackhi_epi16(A_0_15_B_0_15,   C_0_15_D_0_15);
        __m256i E_8_15_H_8_15   = _mm256_unpackhi_epi16(E_0_15_F_0_15,   G_0_15_H_0_15);
        __m256i I_8_15_L_8_15   = _mm256_unpackhi_epi16(I_0_15_J_0_15,   K_0_15_L_0_15);
        __m256i M_8_15_P_8_15   = _mm256_unpackhi_epi16(M_0_15_N_0_15,   O_0_15_P_0_15);

        __m256i A_16_23_D_16_23 = _mm256_unpacklo_epi16(A_16_31_B_16_31, C_16_31_D_16_31);
        __m256i E_16_23_H_16_23 = _mm256_unpacklo_epi16(E_16_31_F_16_31, G_16_31_H_16_31);
        __m256i I_16_23_L_16_23 = _mm256_unpacklo_epi16(I_16_31_J_16_31, K_16_31_L_16_31);
        __m256i M_16_23_P_16_23 = _mm256_unpacklo_epi16(M_16_31_N_16_31, O_16_31_P_16_31);

        __m256i A_24_31_D_24_31 = _mm256_unpackhi_epi16(A_16_31_B_16_31, C_16_31_D_16_31);
        __m256i E_24_31_H_24_31 = _mm256_unpackhi_epi16(E_16_31_F_16_31, G_16_31_H_16_31);
        __m256i I_24_31_L_24_31 = _mm256_unpackhi_epi16(I_16_31_J_16_31, K_16_31_L_16_31);
        __m256i M_24_31_P_24_31 = _mm256_unpackhi_epi16(M_16_31_N_16_31, O_16_31_P_16_31);

        //
        // STAGE 4 DU BUTTERFLY
        //
        __m256i a_0_3_h_0_3      = _mm256_unpacklo_epi32(a_0_7_d_0_7,     e_0_7_h_0_7); // LOW
        __m256i a_4_7_h_4_7      = _mm256_unpackhi_epi32(a_0_7_d_0_7,     e_0_7_h_0_7); // HIGH
        __m256i a_8_11_h_8_11    = _mm256_unpacklo_epi32(a_8_15_d_8_15,   e_8_15_h_8_15); // LOW
        __m256i a_12_15_h_12_15  = _mm256_unpackhi_epi32(a_8_15_d_8_15,   e_8_15_h_8_15); // HIGH
        __m256i a_16_19_h_16_19  = _mm256_unpacklo_epi32(a_16_23_d_16_23, e_16_23_h_16_23); // LOW
        __m256i a_20_23_h_20_23  = _mm256_unpackhi_epi32(a_16_23_d_16_23, e_16_23_h_16_23); // HIGH
        __m256i a_24_27_h_24_27  = _mm256_unpacklo_epi32(a_24_31_d_24_31, e_24_31_h_24_31); // LOW
        __m256i a_28_31_h_28_31  = _mm256_unpackhi_epi32(a_24_31_d_24_31, e_24_31_h_24_31); // HIGH
        __m256i i_0_3_p_0_3      = _mm256_unpacklo_epi32(i_0_7_l_0_7,     m_0_7_p_0_7); // LOW
        __m256i i_4_7_p_4_7      = _mm256_unpackhi_epi32(i_0_7_l_0_7,     m_0_7_p_0_7); // HIGH
        __m256i i_8_11_p_8_11    = _mm256_unpacklo_epi32(i_8_15_l_8_15,   m_8_15_p_8_15); // LOW
        __m256i i_12_15_p_12_15  = _mm256_unpackhi_epi32(i_8_15_l_8_15,   m_8_15_p_8_15); // HIGH
        __m256i i_16_19_p_16_19  = _mm256_unpacklo_epi32(i_16_23_l_16_23, m_16_23_p_16_23); // LOW
        __m256i i_20_23_p_20_23  = _mm256_unpackhi_epi32(i_16_23_l_16_23, m_16_23_p_16_23); // HIGH
        __m256i i_24_27_p_24_27  = _mm256_unpacklo_epi32(i_24_31_l_24_31, m_24_31_p_24_31); // LOW
        __m256i i_28_31_p_28_31  = _mm256_unpackhi_epi32(i_24_31_l_24_31, m_24_31_p_24_31); // HIGH

        __m256i A_0_3_H_0_3      = _mm256_unpacklo_epi32(A_0_7_D_0_7,     E_0_7_H_0_7);     // LOW
        __m256i A_4_7_H_4_7      = _mm256_unpackhi_epi32(A_0_7_D_0_7,     E_0_7_H_0_7);     // HIGH
        __m256i A_8_11_H_8_11    = _mm256_unpacklo_epi32(A_8_15_D_8_15,   E_8_15_H_8_15);     // LOW
        __m256i A_12_15_H_12_15  = _mm256_unpackhi_epi32(A_8_15_D_8_15,   E_8_15_H_8_15);     // HIGH
        __m256i A_16_19_H_16_19  = _mm256_unpacklo_epi32(A_16_23_D_16_23, E_16_23_H_16_23); // LOW
        __m256i A_20_23_H_20_23  = _mm256_unpackhi_epi32(A_16_23_D_16_23, E_16_23_H_16_23); // HIGH
        __m256i A_24_27_H_24_27  = _mm256_unpacklo_epi32(A_24_31_D_24_31, E_24_31_H_24_31); // LOW
        __m256i A_28_31_H_28_31  = _mm256_unpackhi_epi32(A_24_31_D_24_31, E_24_31_H_24_31); // HIGH
        __m256i I_0_3_P_0_3      = _mm256_unpacklo_epi32(I_0_7_L_0_7,     M_0_7_P_0_7);     // LOW
        __m256i I_4_7_P_4_7      = _mm256_unpackhi_epi32(I_0_7_L_0_7,     M_0_7_P_0_7);     // HIGH
        __m256i I_8_11_P_8_11    = _mm256_unpacklo_epi32(I_8_15_L_8_15,   M_8_15_P_8_15);     // LOW
        __m256i I_12_15_P_12_15  = _mm256_unpackhi_epi32(I_8_15_L_8_15,   M_8_15_P_8_15);     // HIGH
        __m256i I_16_19_P_16_19  = _mm256_unpacklo_epi32(I_16_23_L_16_23, M_16_23_P_16_23); // LOW
        __m256i I_20_23_P_20_23  = _mm256_unpackhi_epi32(I_16_23_L_16_23, M_16_23_P_16_23); // HIGH
        __m256i I_24_27_P_24_27  = _mm256_unpacklo_epi32(I_24_31_L_24_31, M_24_31_P_24_31); // LOW
        __m256i I_28_31_P_28_31  = _mm256_unpackhi_epi32(I_24_31_L_24_31, M_24_31_P_24_31); // HIGH

        //
        // STAGE 5 DU BUTTERFLY
        //
        __m256i a_0_1_p_0_1      = _mm256_unpacklo_epi64(a_0_3_h_0_3,     i_0_3_p_0_3);
        __m256i a_2_3_p_2_3      = _mm256_unpackhi_epi64(a_0_3_h_0_3,     i_0_3_p_0_3);
        __m256i a_4_5_p_4_5      = _mm256_unpacklo_epi64(a_4_7_h_4_7,     i_4_7_p_4_7);
        __m256i a_6_7_p_6_7      = _mm256_unpackhi_epi64(a_4_7_h_4_7,     i_4_7_p_4_7);
        __m256i a_8_9_p_8_9      = _mm256_unpacklo_epi64(a_8_11_h_8_11,   i_8_11_p_8_11);
        __m256i a_10_11_p_10_11  = _mm256_unpackhi_epi64(a_8_11_h_8_11,   i_8_11_p_8_11);
        __m256i a_12_13_p_12_13  = _mm256_unpacklo_epi64(a_12_15_h_12_15, i_12_15_p_12_15);
        __m256i a_14_15_p_14_15  = _mm256_unpackhi_epi64(a_12_15_h_12_15, i_12_15_p_12_15);
        __m256i a_16_17_p_16_17  = _mm256_unpacklo_epi64(a_16_19_h_16_19, i_16_19_p_16_19);
        __m256i a_18_19_p_18_19  = _mm256_unpackhi_epi64(a_16_19_h_16_19, i_16_19_p_16_19);
        __m256i a_20_21_p_20_21  = _mm256_unpacklo_epi64(a_20_23_h_20_23, i_20_23_p_20_23);
        __m256i a_22_23_p_22_23  = _mm256_unpackhi_epi64(a_20_23_h_20_23, i_20_23_p_20_23);
        __m256i a_24_25_p_24_25  = _mm256_unpacklo_epi64(a_24_27_h_24_27, i_24_27_p_24_27);
        __m256i a_26_27_p_26_27  = _mm256_unpackhi_epi64(a_24_27_h_24_27, i_24_27_p_24_27);
        __m256i a_28_29_p_28_29  = _mm256_unpacklo_epi64(a_28_31_h_28_31, i_28_31_p_28_31);
        __m256i a_30_31_p_30_31  = _mm256_unpackhi_epi64(a_28_31_h_28_31, i_28_31_p_28_31);

        __m256i A_0_1_P_0_1      = _mm256_unpacklo_epi64(A_0_3_H_0_3,     I_0_3_P_0_3);
        __m256i A_2_3_P_2_3      = _mm256_unpackhi_epi64(A_0_3_H_0_3,     I_0_3_P_0_3);
        __m256i A_4_5_P_4_5      = _mm256_unpacklo_epi64(A_4_7_H_4_7,     I_4_7_P_4_7);
        __m256i A_6_7_P_6_7      = _mm256_unpackhi_epi64(A_4_7_H_4_7,     I_4_7_P_4_7);
        __m256i A_8_9_P_8_9      = _mm256_unpacklo_epi64(A_8_11_H_8_11,   I_8_11_P_8_11);
        __m256i A_10_11_P_10_11  = _mm256_unpackhi_epi64(A_8_11_H_8_11,   I_8_11_P_8_11);
        __m256i A_12_13_P_12_13  = _mm256_unpacklo_epi64(A_12_15_H_12_15, I_12_15_P_12_15);
        __m256i A_14_15_P_14_15  = _mm256_unpackhi_epi64(A_12_15_H_12_15, I_12_15_P_12_15);
        __m256i A_16_17_P_16_17  = _mm256_unpacklo_epi64(A_16_19_H_16_19, I_16_19_P_16_19);
        __m256i A_18_19_P_18_19  = _mm256_unpackhi_epi64(A_16_19_H_16_19, I_16_19_P_16_19);
        __m256i A_20_21_P_20_21  = _mm256_unpacklo_epi64(A_20_23_H_20_23, I_20_23_P_20_23);
        __m256i A_22_23_P_22_23  = _mm256_unpackhi_epi64(A_20_23_H_20_23, I_20_23_P_20_23);
        __m256i A_24_25_P_24_25  = _mm256_unpacklo_epi64(A_24_27_H_24_27, I_24_27_P_24_27);
        __m256i A_26_27_P_26_27  = _mm256_unpackhi_epi64(A_24_27_H_24_27, I_24_27_P_24_27);
        __m256i A_28_29_P_28_29  = _mm256_unpacklo_epi64(A_28_31_H_28_31, I_28_31_P_28_31);
        __m256i A_30_31_P_30_31  = _mm256_unpackhi_epi64(A_28_31_H_28_31, I_28_31_P_28_31);
        
  
        //
        // STAGE 6 DU BUTTERFLY
        //
        __m256i a_0_p_0_A_0_P_0     = _mm256_unpacklo_epi128(a_0_1_p_0_1,     A_0_1_P_0_1);
        __m256i a_1_p_1_A_1_P_1     = _mm256_unpackhi_epi128(a_0_1_p_0_1,     A_0_1_P_0_1);

        __m256i a_2_p_2_A_2_P_2     = _mm256_unpacklo_epi128(a_2_3_p_2_3,     A_2_3_P_2_3);
        __m256i a_3_p_3_A_3_P_3     = _mm256_unpackhi_epi128(a_2_3_p_2_3,     A_2_3_P_2_3);

        __m256i a_4_p_4_A_4_P_4     = _mm256_unpacklo_epi128(a_4_5_p_4_5,     A_4_5_P_4_5);
        __m256i a_5_p_5_A_5_P_5     = _mm256_unpackhi_epi128(a_4_5_p_4_5,     A_4_5_P_4_5);

        __m256i a_6_p_6_A_6_P_6     = _mm256_unpacklo_epi128(a_6_7_p_6_7,     A_6_7_P_6_7);
        __m256i a_7_p_7_A_7_P_7     = _mm256_unpackhi_epi128(a_6_7_p_6_7,     A_6_7_P_6_7);

        __m256i a_8_p_8_A_8_P_8     = _mm256_unpacklo_epi128(a_8_9_p_8_9,     A_8_9_P_8_9);
        __m256i a_9_p_9_A_9_P_9     = _mm256_unpackhi_epi128(a_8_9_p_8_9,     A_8_9_P_8_9);

        __m256i a_10_p_10_A_10_P_12 = _mm256_unpacklo_epi128(a_10_11_p_10_11, A_10_11_P_10_11);
        __m256i a_11_p_11_A_11_P_12 = _mm256_unpackhi_epi128(a_10_11_p_10_11, A_10_11_P_10_11);

        __m256i a_12_p_12_A_12_P_12 = _mm256_unpacklo_epi128(a_12_13_p_12_13, A_12_13_P_12_13);
        __m256i a_13_p_13_A_13_P_13 = _mm256_unpackhi_epi128(a_12_13_p_12_13, A_12_13_P_12_13);

        __m256i a_14_p_14_A_14_P_14 = _mm256_unpacklo_epi128(a_14_15_p_14_15, A_14_15_P_14_15);
        __m256i a_15_p_15_A_15_P_15 = _mm256_unpackhi_epi128(a_14_15_p_14_15, A_14_15_P_14_15);

        __m256i a_16_p_16_A_16_P_16 = _mm256_unpacklo_epi128(a_16_17_p_16_17, A_16_17_P_16_17);
        __m256i a_17_p_17_A_17_P_17 = _mm256_unpackhi_epi128(a_16_17_p_16_17, A_16_17_P_16_17);

        __m256i a_18_p_18_A_18_P_18 = _mm256_unpacklo_epi128(a_18_19_p_18_19, A_18_19_P_18_19);
        __m256i a_19_p_19_A_19_P_19 = _mm256_unpackhi_epi128(a_18_19_p_18_19, A_18_19_P_18_19);

        __m256i a_20_p_20_A_20_P_20 = _mm256_unpacklo_epi128(a_20_21_p_20_21, A_20_21_P_20_21);
        __m256i a_21_p_21_A_21_P_21 = _mm256_unpackhi_epi128(a_20_21_p_20_21, A_20_21_P_20_21);

        __m256i a_22_p_22_A_22_P_22 = _mm256_unpacklo_epi128(a_22_23_p_22_23, A_22_23_P_22_23);
        __m256i a_23_p_23_A_23_P_23 = _mm256_unpackhi_epi128(a_22_23_p_22_23, A_22_23_P_22_23);

        __m256i a_24_p_24_A_24_P_24 = _mm256_unpacklo_epi128(a_24_25_p_24_25, A_24_25_P_24_25);
        __m256i a_25_p_25_A_25_P_25 = _mm256_unpackhi_epi128(a_24_25_p_24_25, A_24_25_P_24_25);

        __m256i a_26_p_26_A_26_P_26 = _mm256_unpacklo_epi128(a_26_27_p_26_27, A_26_27_P_26_27);
        __m256i a_27_p_27_A_27_P_27 = _mm256_unpackhi_epi128(a_26_27_p_26_27, A_26_27_P_26_27);

        __m256i a_28_p_28_A_28_P_28 = _mm256_unpacklo_epi128(a_28_29_p_28_29, A_28_29_P_28_29);
        __m256i a_29_p_29_A_29_P_29 = _mm256_unpackhi_epi128(a_28_29_p_28_29, A_28_29_P_28_29);

        __m256i a_30_p_30_A_30_P_30 = _mm256_unpacklo_epi128(a_30_31_p_30_31, A_30_31_P_30_31);
        __m256i a_31_p_31_A_31_P_31 = _mm256_unpackhi_epi128(a_30_31_p_30_31, A_30_31_P_30_31);

        STORE_SIMD_FX(p_output +  0 * constN, a_0_p_0_A_0_P_0);
        STORE_SIMD_FX(p_output +  1 * constN, a_2_p_2_A_2_P_2);
        STORE_SIMD_FX(p_output +  2 * constN, a_4_p_4_A_4_P_4);
        STORE_SIMD_FX(p_output +  3 * constN, a_6_p_6_A_6_P_6);
        STORE_SIMD_FX(p_output +  4 * constN, a_8_p_8_A_8_P_8);
        STORE_SIMD_FX(p_output +  5 * constN, a_10_p_10_A_10_P_12);
        STORE_SIMD_FX(p_output +  6 * constN, a_12_p_12_A_12_P_12);
        STORE_SIMD_FX(p_output +  7 * constN, a_14_p_14_A_14_P_14);
        STORE_SIMD_FX(p_output +  8 * constN, a_16_p_16_A_16_P_16);
        STORE_SIMD_FX(p_output +  9 * constN, a_18_p_18_A_18_P_18);
        STORE_SIMD_FX(p_output + 10 * constN, a_20_p_20_A_20_P_20);
        STORE_SIMD_FX(p_output + 11 * constN, a_22_p_22_A_22_P_22);
        STORE_SIMD_FX(p_output + 12 * constN, a_24_p_24_A_24_P_24);
        STORE_SIMD_FX(p_output + 13 * constN, a_26_p_26_A_26_P_26);
        STORE_SIMD_FX(p_output + 14 * constN, a_28_p_28_A_28_P_28);
        STORE_SIMD_FX(p_output + 15 * constN, a_30_p_30_A_30_P_30);
        
        STORE_SIMD_FX(p_output + 16 * constN, a_1_p_1_A_1_P_1);
        STORE_SIMD_FX(p_output + 17 * constN, a_3_p_3_A_3_P_3);
        STORE_SIMD_FX(p_output + 18 * constN, a_5_p_5_A_5_P_5);
        STORE_SIMD_FX(p_output + 19 * constN, a_7_p_7_A_7_P_7);
        STORE_SIMD_FX(p_output + 20 * constN, a_9_p_9_A_9_P_9);
        STORE_SIMD_FX(p_output + 21 * constN, a_11_p_11_A_11_P_12);
        STORE_SIMD_FX(p_output + 22 * constN, a_13_p_13_A_13_P_13);
        STORE_SIMD_FX(p_output + 23 * constN, a_15_p_15_A_15_P_15);
        STORE_SIMD_FX(p_output + 24 * constN, a_17_p_17_A_17_P_17);
        STORE_SIMD_FX(p_output + 25 * constN, a_19_p_19_A_19_P_19);
        STORE_SIMD_FX(p_output + 26 * constN, a_21_p_21_A_21_P_21);
        STORE_SIMD_FX(p_output + 27 * constN, a_23_p_23_A_23_P_23);
        STORE_SIMD_FX(p_output + 28 * constN, a_25_p_25_A_25_P_25);
        STORE_SIMD_FX(p_output + 29 * constN, a_27_p_27_A_27_P_27);
        STORE_SIMD_FX(p_output + 30 * constN, a_29_p_29_A_29_P_29);
        STORE_SIMD_FX(p_output + 31 * constN, a_31_p_31_A_31_P_31);

        p_output += 1;
        p_input  += 32;
    }    
//IACA_END 
}
