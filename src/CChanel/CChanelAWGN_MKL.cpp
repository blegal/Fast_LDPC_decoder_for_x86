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

#include "CChanelAWGN_MKL.h"


#ifndef VSL_METHOD_SGAUSSIAN_BOXMULLER2
	#define VSL_METHOD_SGAUSSIAN_BOXMULLER2 1
#endif

//
// RACINE CARREE SSE OPTIMISEE A L'AIDE DE LA FONCTION
// RECIPROQUE (PRECISION DE 11 BITS SUR LA MANTISSE)
//
//inline __m256 sqrt_sse_11bits( __m256 a )
//{
//   return _mm256_mul_ps( a, _mm256_rsqrt_ps( a ) );
//}

double CChanelAWGN_MKL::inv_erf(int v){
    if (v == 3) {
        return 0.86312;
    }else if(v == 4){
        return 1.1064;
    }else if(v == 5){
        return 1.3268;
    }else if(v == 6){
        return 1.5274;
    }else if(v == 7){
        return 1.7115;
    }else if(v == 8){
        return 1.8819;
    }else if(v == 9){
        return 2.0409;
    }else if(v == 10){
        return 2.1903;
    }
    return -1;
}

double CChanelAWGN_MKL::get_R(){
    return R;
}

#define AVX_8F_LOAD(ptr)      (_mm256_load_ps(ptr))
#define AVX_8F_STORE(ptr,val) (_mm256_store_ps(ptr,val))
#define AVX_8F_SQRT(a)        (_mm256_sqrt_ps(a))
#define AVX_8F_ADD(a,b)       (_mm256_add_ps(a,b))
#define AVX_8F_SUB(a,b)       (_mm256_sub_ps(a,b))
#define AVX_8F_MUL(a,b)       (_mm256_mul_ps(a,b))
#define AVX_8F_LOG(a)         (_mm256_log_ps(a))
#define AVX_8F_DIV(a,b)       (_mm256_div_ps(a,b))
#define AVX_8F_SET1(a)        (_mm256_set1_ps(a))
#define AVX_8F_SET1i(a)       (_mm256_set1_epi32(a))
#define AVX_8F_CONV(a)        (_mm256_cvtepi32_ps(a))
#define AVX_8F_SETi(a,b,c,d,e,f,g,h)  (_mm256_set_epi32(a,b,c,d,e,f,g,h))

static int thread_id = 0;


CChanelAWGN_MKL::CChanelAWGN_MKL(CTrame *t, int _BITS_LLR, bool QPSK, bool Es_N0)
    : CChanel(t, _BITS_LLR, QPSK, Es_N0){
    int status = vslNewStream( &stream, VSL_BRNG_MT2203 + thread_id++ /*VSL_BRNG_MT2203*/, rand() );
    if( status != VSL_STATUS_OK ){
        printf("(EE) Error during vslNewStream execution\n");
        printf("(EE) thread_id = %d\n", thread_id);
        exit( 0 );
    }
    noise  = (float*)new __m128[_frames * _data / 4];
}

CChanelAWGN_MKL::~CChanelAWGN_MKL(){
    vslDeleteStream( &stream );
    delete noise;
    thread_id--;
}

void CChanelAWGN_MKL::configure(double _Eb_N0){
    
    rendement = (float) (_vars) / (float) (_data);
    if (es_n0) {
        // ES/N0 = Eb/N0 + 10*log10(R*m)
        // o√π R  = rendement
        // m     = nombre de bits par symbole de constellation (QPSK => 2)
        // Eb/N0 et ES/N0 sont en dB
        Eb_N0 = _Eb_N0 - 10.0 * log10(2 * rendement);
    } else {
        Eb_N0 = _Eb_N0;
    }    
    
    double interm = 10.0 * log10(rendement);
    interm        = -0.1*((double)Eb_N0+interm);
    SigB          = sqrt(pow(10.0,interm)/2);
    qbeta         = SigB * sqrt(2.0) * inv_erf( BITS_LLR - 1 ); // PATCH CEDRIC MARCHAND
    R             = (1.0 + qbeta);

    //
    // FACTEUR DE NORMALISATION DU CANAL PROVENENT DU DECODEUR DE CODE POLAIRE
    // A CAMILLE.
    //
    if( normalize == true ){
        norm_factor  = 2.0 / (SigB * SigB);
    }else{
        norm_factor  = 1.0;
    }
}

#define QPSK 0.707106781
#define BPSK 1.0

void CChanelAWGN_MKL::generate() {
    float pv = (qpsk) ?  QPSK :  BPSK; // ON CHOISIT LE TYPE DE CODAGE DU SIGNAL
    float mv = (qpsk) ? -QPSK : -BPSK; // BPSK OU QPSK (CODAGES LES + SIMPLES)

    int nbData = (_frames*_data);
    vsRngGaussian( VSL_METHOD_SGAUSSIAN_BOXMULLER2, stream, nbData, (float*)noise, 0.0f, SigB );        

    //
    // ON LAISSE ICC DEROULER LA BOUCLE AINSI C'EST SSE ET AVX COMPATIBLE (4 TIBO)
    //
    for (int z = 0; z < _frames*_data; z++) {
        float y = (t_coded_bits[z] == 1 ? pv : mv) + noise[z];
        t_noise_data[z] = y * norm_factor;

    }
}
#endif
