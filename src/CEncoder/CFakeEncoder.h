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

#ifndef CLASS_CFakeEncoder
#define CLASS_CFakeEncoder

#include <stdlib.h>
#include <stdio.h>
#include "Encoder.h"

class CFakeEncoder : public Encoder
{    
//protected:
//    int  _vars;
//    int  _data;
//    int*  t_in_bits;      // taille (var)
//    int*  t_coded_bits;   // taille (data)
    
public:
    
    CFakeEncoder(CTrame *t);
    
    virtual void encode();
};

#endif
