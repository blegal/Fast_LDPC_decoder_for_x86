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

#ifndef __CFakeDecoder__
#define __CFakeDecoder__

class CFakeDecoder{

public:
    CFakeDecoder();
    CFakeDecoder(int data);
    virtual void decode(double var_nodes[N], int Rprime_fix[N], int nombre_iterations);
};


CFakeDecoder::CFakeDecoder()
{
}

void CFakeDecoder::decode(double var_nodes[N], int Rprime_fix[N], int nombre_iterations)
{
    int i;
    for (i=0;i<N;i++){
    	Rprime_fix[i] = (var_nodes[i] > 0.0);
    }
}

#endif;
