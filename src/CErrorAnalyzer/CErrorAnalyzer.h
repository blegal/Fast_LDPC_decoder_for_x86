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

#ifndef CLASS_CError_Analyzer
#define CLASS_CError_Analyzer

#include "../CTrame/CTrame.h"

class CErrorAnalyzer
{
private:
	long int nb_bit_errors;
	long int nb_frame_errors;
	long int nb_analyzed_frames;
    int *buf_en_bits;

protected:
    int   _vars;
    int   _data;
    int   _frames;
    int*  t_in_bits;
    char* t_decode_data;

    int _max_fe;
    bool _auto_fe_mode;
    bool _worst_case;

    int* faulty_bits;
    
public:
    CErrorAnalyzer(CTrame *t);
    CErrorAnalyzer(CTrame *t, int max_fe);
    CErrorAnalyzer(CTrame *t, int max_fe, bool auto_fe_mode, bool worst_case_fer);
    virtual ~CErrorAnalyzer();
    virtual void generate();
    virtual void store_enc_bits();
    virtual void generate(int nErrors);

    long int nb_processed_frames();
    long int nb_fe();
    long int nb_be();

    long int fe_limit();
    bool fe_limit_achieved();

    double fer_value();
    double ber_value();

    long int nb_processed_frames(int add);
    long int nb_fe(int add);
    long int nb_be(int add);

    int nb_data();
    int nb_vars();
    int nb_checks();
};

#endif
