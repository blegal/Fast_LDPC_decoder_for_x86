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

#include "CBitGenerator.h"

CBitGenerator::CBitGenerator(CTrame *t, bool zero_only){
    _vars        = t->nb_vars();
    t_in_bits    = t->get_t_in_bits();
    _zero_mode   = zero_only;

    for(int i=0; i<_vars; i++){
        t_in_bits[i] = 0;
    }
}

void CBitGenerator::generate(){
    if( _zero_mode == false ){
        for(int i=0; i<_vars; i++){
            t_in_bits[i] = rand()%2;
        }
    }
}
