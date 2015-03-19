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

#ifndef CLASS_CDecoder_NMS_
#define CLASS_CDecoder_NMS_

#include "../template/CDecoder_fixed_SSE.h"

class CDecoder_NMS_fixed_SSE : public CDecoder_fixed_SSE{
protected:
    int factor_1;
    int factor_2;
    __m128i **p_vn_adr;
    
public:
    CDecoder_NMS_fixed_SSE();
    ~CDecoder_NMS_fixed_SSE();
    void setFactor(int _factor);
    void decode        (char var_nodes[], char Rprime_fix[], int nombre_iterations);
private:
    void decode_8bits  (char var_nodes[], char Rprime_fix[], int nombre_iterations);
};

#endif
