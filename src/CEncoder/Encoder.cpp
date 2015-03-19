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

#include "Encoder.h"

Encoder::Encoder(CTrame *t)
{
    _size_in   = t->nb_vars();
    _size_out  = t->nb_data();
    data_in    = t->get_t_in_bits();
    data_out   = t->get_t_coded_bits();
    _frames    = t->nb_frames();
}

Encoder::~Encoder()
{
}

int Encoder::size_in()
{
    return _size_in;
}


int Encoder::size_out()
{
    return _size_out;
}


void Encoder::encode()
{
    for(int y=0; y<_size_out; y++)
    {
        data_out[y] = 0;
    }
}


void Encoder::sum_bits(){
    int sum = 0;
    for(int y=0; y<_size_out; y++)
    {
        sum += data_out[y];
    }
    printf("SOMME BITS = %d\n", sum);
}


void Encoder::sum_pos(){
    int sum = 0;
    for(int y=0; y<_size_out; y++)
    {
        sum += (data_out[y]!=0)?(y+1):0;
    }
    printf("SOMME POSI = %d\n", sum);
}
