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

#ifndef CLASS_CChanelAWGN_MKL
#define CLASS_CChanelAWGN_MKL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "CChanel.h"
#include <xmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include "mkl.h"

#define small_pi  3.1415926536
#define _2pi  (2.0 * small_pi)

class CChanelAWGN_MKL : public CChanel
{
private:
    double inv_erf(int v);
    float* noise;
    VSLStreamStatePtr stream;
    
public:
    CChanelAWGN_MKL(CTrame *t, int _BITS_LLR, bool QPSK, bool Es_N0);
    ~CChanelAWGN_MKL();
    
    virtual void   configure(double _Eb_N0);
    virtual double get_R();
    virtual void   generate();
};

#endif

