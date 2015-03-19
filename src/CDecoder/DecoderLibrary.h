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


#ifndef DECODERLIBRARY_H
#define	DECODERLIBRARY_H

#include <iostream>
using namespace std;

//
// DECODEUR MS (MIN-SUM)
//
#include "./OMS/CDecoder_OMS_fixed_SSE.h"
#include "./OMS/CDecoder_OMS_fixed_AVX.h"

#include "./NMS/CDecoder_NMS_fixed_SSE.h"
#include "./NMS/CDecoder_NMS_fixed_AVX.h"

#define UNAVAILABLE { cout << "(EE) Error, decoder unavailable"  << endl; \
                      cout << "(EE) - decoder type : " << type   << endl; \
                      cout << "(EE) - architecture : " << arch   << endl; \
                      cout << "(EE) - data format  : " << format << endl; \
                      exit( 0 ); }


CDecoder* CreateDecoder(
        string type,    //
        string arch,    //
        string format,  //
        param_decoder p_decoder,
        int vMin,       //
        int vMax,       //
        int mMin,       //
        int mMax
    ) {
    
    ////////////////////////////////////////////////////////////////////////////
    //
    // GESTION DES IMPLANTATIONS DE L'ALGORITME NULL
    //
    if (type.compare("NULL") == 0) {

        UNAVAILABLE;

        
    ////////////////////////////////////////////////////////////////////////////
    //
    // GESTION DES IMPLANTATIONS DE L'ALGORITME OFFSET MIN-SUM
    //
    } else if (type.compare("OMS") == 0) {

        if (format.compare("fixed") == 0 && arch.compare("sse") == 0) {
            CDecoder_OMS_fixed_SSE* dec = new CDecoder_OMS_fixed_SSE();
            dec->setOffset   ( p_decoder.oms_offset_fixed );
            dec->setVarRange(vMin, vMax);
            dec->setMsgRange(mMin, mMax);
            return dec;

        } else if (format.compare("fixed") == 0 && arch.compare("avx") == 0) {
            CDecoder_OMS_fixed_AVX* dec = new CDecoder_OMS_fixed_AVX();
            dec->setOffset   ( p_decoder.oms_offset_fixed );
            dec->setVarRange(vMin, vMax);
            dec->setMsgRange(mMin, mMax);
            return dec;
        } else {
                UNAVAILABLE;
        }



    ////////////////////////////////////////////////////////////////////////////
    //
    // GESTION DES IMPLANTATIONS DE L'ALGORITME OFFSET MIN-SUM
    //
    } else if (type.compare("NMS") == 0) {

        //
        // DECODEURS FLOTTANTS
        //
        if (format.compare("fixed") == 0 && arch.compare("sse") == 0) {
            CDecoder_NMS_fixed_SSE* dec = new CDecoder_NMS_fixed_SSE();
            dec->setFactor   ( p_decoder.nms_factor_fixed );
            dec->setVarRange(vMin, vMax);
            dec->setMsgRange(mMin, mMax);
            return dec;

        } else if (format.compare("fixed") == 0 && arch.compare("avx") == 0) {
            CDecoder_NMS_fixed_AVX* dec = new CDecoder_NMS_fixed_AVX();
            dec->setFactor   ( p_decoder.nms_factor_fixed );
            dec->setVarRange (vMin, vMax);
            dec->setMsgRange (mMin, mMax);
            return dec;

        } else {
                UNAVAILABLE;
        }

    } // FIN DE LA LISTE DES DECODEURS

    cout << "(EE) Requested LDPC decoder does not exist (" << arch << ":" << type << ")" << endl;
    exit(0);
    return NULL;
}

#endif	/* DECODERLIBRARY_H */

