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

#include "CChanel.h"

double CChanel::get_R(){
    return R;
}

double CChanel::get_SigB(){
    return SigB;
}

CChanel::~CChanel(){
}


CChanel::CChanel(CTrame *t, int _BITS_LLR, bool QPSK, bool ES_N0){
    qbeta        = 0.0;
    R            = 0.0;
    _vars        = t->nb_vars();
    _data        = t->nb_data();
    _checks      = t->nb_checks();
    t_coded_bits = t->get_t_coded_bits();
    t_noise_data = t->get_t_noise_data();
    _frames      = t->nb_frames();
    BITS_LLR     = _BITS_LLR;
    qpsk         = QPSK;
    es_n0        = ES_N0;
    normalize    = false;
    norm_factor  = 0.0f;
}

void CChanel::setNormalize(bool enable){
    normalize = enable;
}
