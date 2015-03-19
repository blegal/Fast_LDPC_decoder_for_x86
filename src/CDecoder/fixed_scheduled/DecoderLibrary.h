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

#include "CDecoder_fixed.h"

#include "CDecoder_fixed_x86.h"
#include "CDecoder_x86_NULL.h"
#include "CDecoder_MS_fixed_x86.h"
#include "CDecoder_x86_OMS.h"
#include "CDecoder_x86_NMS.h"

#include "CDecoder_fixed_SSE.h"
#include "CDecoder_SSE_NULL.h"
#include "CDecoder_MS_fixed_SSE.h"
#include "CDecoder_SSE_OMS.h"
#include "CDecoder_SSE_NMS.h"

CDecoder_fixed* CreateDecoder(string type, string arch, param_decoder p_decoder, int vMin, int vMax, int mMin, int mMax, int offset, int factor) {

    //
    //
    //
    if( arch.compare("x86") == 0 ){

        if ( type.compare("MS") == 0 ) {
            CDecoder_MS_fixed_x86* dec = new CDecoder_MS_fixed_x86();
            dec->setEarlyTerm( p_decoder.early_term );
            dec->setVarRange(vMin, vMax);
            dec->setMsgRange(mMin, mMax);
            return dec;

        } else if (type.compare("OMS") == 0) {
            CDecoder_OMS_fixed_x86* dec = new CDecoder_OMS_fixed_x86();
            dec->setEarlyTerm( p_decoder.early_term );
            dec->setVarRange(vMin, vMax);
            dec->setMsgRange(mMin, mMax);
            dec->setOffset  ( offset );
            return dec;

        } else if (type.compare("NMS") == 0) {
            CDecoder_NMS_fixed_x86* dec = new CDecoder_NMS_fixed_x86();
            dec->setEarlyTerm( p_decoder.early_term );
            dec->setVarRange(vMin, vMax);
            dec->setMsgRange(mMin, mMax);
            dec->setFactor  ( factor );
            return dec;

        } else if (type.compare("NULL") == 0) {
            dec->setEarlyTerm( p_decoder.early_term );
            return new CDecoder_NULL_fixed_x86();
        }

    //
    //
    //
    }else if( arch.compare("sse") == 0 ){

        if ( type.compare("MS") == 0 ) {
            CDecoder_MS_fixed_SSE* dec = new CDecoder_MS_fixed_SSE();
            dec->setVarRange(vMin, vMax);
            dec->setMsgRange(mMin, mMax);
            return dec;

        } else if (type.compare("OMS") == 0) {
            CDecoder_OMS_fixed_SSE* dec = new CDecoder_OMS_fixed_SSE();
            dec->setVarRange(vMin, vMax);
            dec->setMsgRange(mMin, mMax);
            dec->setOffset  ( offset );
            return dec;

        } else if (type.compare("NMS") == 0) {
            CDecoder_NMS_fixed_SSE* dec = new CDecoder_NMS_fixed_SSE();
            dec->setVarRange(vMin, vMax);
            dec->setMsgRange(mMin, mMax);
            dec->setFactor  ( factor );
            return dec;

        } else if (type.compare("NULL") == 0) {
            return new CDecoder_NULL_fixed_SSE();

        }
    }
    cout << "(EE) Requested LDPC decoder does not exist (" << arch << ":" << type << ")" << endl;
    exit( 0 );
    return NULL;
}

#endif	/* DECODERLIBRARY_H */
