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

#include "CFastFixConversion.h"


CFastFixConversion::CFastFixConversion(CTrame *t, int _FACTEUR_BETA, int _vSAT_NEG_LLR, int _vSAT_POS_LLR) : CFixConversion(t) {
    FACTEUR_BETA  = _FACTEUR_BETA;
    vSAT_NEG_LLR  = _vSAT_NEG_LLR;
    vSAT_POS_LLR  = _vSAT_POS_LLR;
}

CFastFixConversion::~CFastFixConversion(){
    if( dump_q ){
        long long sum = 0;
        for(int i=0; i<256; i++){
            if( ((i-128) <= (vSAT_POS_LLR)) && ((i-128) >= (vSAT_NEG_LLR)) ){
                sum += histo[i];
            }
        }
        printf("(HISTO) START\n");
        for(int i=0; i<256; i++){
            if( ((i-128) <= (vSAT_POS_LLR+1)) && ((i-128) >= (vSAT_NEG_LLR-1)) ){
                double proba = 100.0 * ((double)histo[i]) / (double)sum;
                printf("(HISTO) %3d\t%f\n", i-128, proba);
            }
        }
        printf("(HISTO) STOP\n");
    }
}

#include <emmintrin.h>
#include <smmintrin.h>
#include <smmintrin.h>

#define OPTIMIZATION_SSE 0

void CFastFixConversion::generate(){
#if OPTIMIZATION_SSE == 0
    for(int z = 0; z<_frames; z++){
        int offset = z * _data;
	    for(int i=0; i<_data; i++){
	        int value = (FACTEUR_BETA * t_noise_data[offset+i]);  // ON TRANFORME LES FLOTTANT EN FIXED POINT
	        value = (value>vSAT_NEG_LLR) ? value : vSAT_NEG_LLR;  // ON GERE LA SATURATION DES DONNEES DANS LE
	        value = (value<vSAT_POS_LLR) ? value : vSAT_POS_LLR;  // FORMAT DEDIE AU LLR
	        t_fpoint_data[offset+i] = value;
	    }
	}
#else    
    int nb_data = _frames * _data;
    // PEUT'ON FONCTIONNER EN MODE SSE ?
    if( (nb_data & 0x03) == 0 ){

        __m128i mini  = _mm_set1_epi32(vSAT_NEG_LLR);
        __m128i maxi  = _mm_set1_epi32(vSAT_POS_LLR);
        __m128  fact  = _mm_set1_ps(FACTEUR_BETA);
        __m128*  tabI = (__m128* ) t_noise_data;
        __m128i* tabO = (__m128i*) t_fpoint_data;
        int loop = nb_data / 4;

        for(int z=0; z<loop; z++){
            __m128  iData = _mm_mul_ps( tabI[z], fact );
            __m128i resul = _mm_cvttps_epi32(iData);
            __m128i satu1 = _mm_max_epi32  ( resul, mini );
            __m128i satu2 = _mm_min_epi32  ( satu1, maxi );
            tabO[z] = satu2;
        }
    }
    // SINON ON LANCE LE CALCUL EN MODE NORMAL
    else {
        for(int z = 0; z<nb_data; z++){
            int value = (FACTEUR_BETA * t_noise_data[z]);         // ON TRANFORME LES FLOTTANT EN FIXED POINT
            value = (value>vSAT_NEG_LLR) ? value : vSAT_NEG_LLR;  // ON GERE LA SATURATION DES DONNEES DANS LE
            value = (value<vSAT_POS_LLR) ? value : vSAT_POS_LLR;  // FORMAT DEDIE AU LLR
            t_fpoint_data[z] = value;
        }
    }
#endif
    if( dump_q ){
        char* ptr = (char*)t_fpoint_data;
        unsigned int length = _frames * _data;
        while(length--) {
            histo[ ((int)(*ptr++)) + 128 ] += 1;
        }
    }
}
