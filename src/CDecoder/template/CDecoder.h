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

#ifndef __CDecoder__
#define __CDecoder__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

class CDecoder{
protected:
    float sigB;
    int   nb_iters;  // Use early termination criteria

public:
    CDecoder();
    virtual ~CDecoder();
    virtual void setSigmaChannel(float _sigB);
    virtual void setNumberOfIterations(int _value);
    virtual void decode(char  var_nodes[], char Rprime_fix[], int nombre_iterations) = 0;
    virtual void decode(float var_nodes[], char Rprime_fix[], int nombre_iterations) = 0;
};

#endif
