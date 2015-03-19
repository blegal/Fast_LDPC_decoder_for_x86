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

#ifdef __AVX2__
#include "CDecoder_fixed_AVX.h"

CDecoder_fixed_AVX::CDecoder_fixed_AVX()
{
    var_nodes = new __m256i[NOEUD];
    var_mesgs = new __m256i[MESSAGE];
}

CDecoder_fixed_AVX::~CDecoder_fixed_AVX()
{
    delete var_nodes;
    delete var_mesgs;
}

void CDecoder_fixed_AVX::decode(float var_nodes[], char Rprime_fix[], int nombre_iterations)
{
    // ON NE FAIT RIEN !
    // CETTE METHODE ASSURE JUSTE LA COMPATIBILITE ENTRE LES CLASSES MANIPULANT
    // DES DONNEES FLOTTANTES ET CELLES MANIPULANT DES DONNEES EN VIRGULE FIXE.
}
#endif
