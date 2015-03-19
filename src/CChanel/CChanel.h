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

#ifndef CLASS_CChanel
#define CLASS_CChanel

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../CTrame/CTrame.h"

#define small_pi  3.1415926536
#define _2pi  (2.0 * small_pi)

class CChanel
{
    
protected:
    int  _vars;
    int  _checks;
    int  _data;
    int  _frames;
    int  BITS_LLR;
    int* data_in;
    int* data_out;
    bool qpsk;
    bool es_n0;

    bool normalize; // Normalize by 2/pow(sigma, 2)
    float norm_factor;
    
    float* t_noise_data;   // taille (width)
    int*   t_coded_bits;   // taille (width)
    
    double rendement;
    double SigB;
    double Gauss;
    double Ph;
    double Qu;
    double Eb_N0;
    double qbeta;
    double R;
    
public:
    CChanel(CTrame *t, int _BITS_LLR, bool QPSK, bool Es_N0);
    virtual ~CChanel();
    virtual void configure(double _Eb_N0) = 0;  // VIRTUELLE PURE
    virtual double get_R();
    virtual double get_SigB();
    virtual void setNormalize(bool enable);
    virtual void generate() = 0;                // VIRTUELLE PURE    
};

#endif

