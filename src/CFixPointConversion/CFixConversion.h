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

#ifndef CLASS_CFixConversion
#define CLASS_CFixConversion

#include <stdlib.h>
#include <stdio.h>
#include "../CTrame/CTrame.h"

class CFixConversion
{
    
protected:
    int     _data;
    int     _frames;
    float*  t_noise_data;      // taille (var)
    char*   t_fpoint_data;     // taille (width)

    bool dump_q;
    long long histo[256];

public:
    CFixConversion(CTrame *t);
    virtual ~CFixConversion();
    virtual void generate();
    void ShowHistoOnDestroy( bool enable );
};

#endif

