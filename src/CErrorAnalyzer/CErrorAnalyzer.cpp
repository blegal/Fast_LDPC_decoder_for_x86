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

#include "CErrorAnalyzer.h"

CErrorAnalyzer::CErrorAnalyzer(CTrame *t){
    _data         = t->nb_data();
    _vars         = t->nb_vars();
    _frames       = t->nb_frames();
    t_decode_data = t->get_t_decode_data();
    t_in_bits     = t->get_t_in_bits();
    nb_bit_errors      = 0;
    nb_frame_errors    = 0;
    nb_analyzed_frames = 0;
    _max_fe            = 200;
    _auto_fe_mode      = true;
    buf_en_bits        = new int[_data * _frames];
    _worst_case        = false;
    faulty_bits        = new int[_data];
    for(int i=0; i<_data; i++){
        faulty_bits[i] = 0;
    }
}

CErrorAnalyzer::CErrorAnalyzer(CTrame *t, int max_fe){
    _data              = t->nb_data();
    _vars              = t->nb_vars();
    _frames            = t->nb_frames();
    t_decode_data      = t->get_t_decode_data();
    t_in_bits          = t->get_t_in_bits();
    nb_bit_errors      = 0;
    nb_frame_errors    = 0;
    nb_analyzed_frames = 0;
    _max_fe            = max_fe;
    _auto_fe_mode      = false;
    buf_en_bits        = new int[_data * _frames];
    _worst_case        = false;
    faulty_bits        = new int[_data];
    for(int i=0; i<_data; i++){
        faulty_bits[i] = 0;
    }
}

CErrorAnalyzer::CErrorAnalyzer(CTrame *t, int max_fe, bool auto_fe_mode, bool worst_case_fer){
    _data              = t->nb_data();
    _vars              = t->nb_vars();
    _frames            = t->nb_frames();
    t_decode_data      = t->get_t_decode_data();
    t_in_bits          = t->get_t_in_bits();
    nb_bit_errors      = 0;
    nb_frame_errors    = 0;
    nb_analyzed_frames = 0;
    _max_fe            = max_fe;
    _auto_fe_mode      = auto_fe_mode;
    buf_en_bits        = new int[_data * _frames];
    _worst_case        = worst_case_fer;

    faulty_bits        = new int[_data];
    for(int i=0; i<_data; i++){
        faulty_bits[i] = 0;
    }
}

#define SHOW_FAULTY_BITS    false
CErrorAnalyzer::~CErrorAnalyzer(){
    delete buf_en_bits;

    if( SHOW_FAULTY_BITS ){
        for(int i=0; i<_data; i++){
            if( faulty_bits[i] != 0 ){
                printf("Faulty bit : VN(%6d) = %4d\n", i, faulty_bits[i]);
            }
        }
    }
    delete faulty_bits;
}


long int CErrorAnalyzer::fe_limit()
{
    #define SEUIL_MINI  100
    if( _auto_fe_mode == false ){
        return _max_fe;
    }else{
        double tBER = ber_value();
        if( tBER < 1.0e-9){
            return (_max_fe/16) >= SEUIL_MINI ? (_max_fe/16) : SEUIL_MINI;
        }else if( tBER < 1.0e-8){
            return (_max_fe/8)  >= SEUIL_MINI ? (_max_fe/8)  : SEUIL_MINI;
        }else if( tBER < 1.0e-7){
            return (_max_fe/4)  >= SEUIL_MINI ? (_max_fe/4)  : SEUIL_MINI;
        }else if( tBER < 1.0e-6){
            return (_max_fe/2)  >= SEUIL_MINI ? (_max_fe/2)  : SEUIL_MINI;
        }else{
            return (_max_fe);
        }
    }
}


bool CErrorAnalyzer::fe_limit_achieved()
{
    return (nb_fe() >= fe_limit());
}

void CErrorAnalyzer::generate(){
    for(int z = 0; z<_frames; z++){
        int offset1 = z * _data;
        int offset2 = z * _vars;
        int nErrors = 0;

        if( _worst_case == false ){
            for (int i = 0; i < _vars; i++) {
                nErrors += (t_decode_data[offset1 + i] != buf_en_bits[offset2 + i]);
            }
        }else{
            for (int i = 0; i < _data; i++) {
                nErrors += (t_decode_data[offset1 + i] != 0);
            }
        }

        //
        // ON FAIT UN HISTOGRAMME DES FAULTY BITS
        //
        if ( nErrors > 0 ) {
            for (int i = 0; i < _data; i++) {
                if (t_decode_data[offset1 + i] != 0) {
                    faulty_bits[i] += 1;
                }
            }
        }
        
        nb_bit_errors      += nErrors;
        nb_frame_errors    += (nErrors != 0);
        nb_analyzed_frames += 1;
    }
}

void CErrorAnalyzer::generate(int nErrors){
    nb_bit_errors      += nErrors;
    nb_frame_errors    += (nErrors != 0);
    nb_analyzed_frames += 1;
}

void CErrorAnalyzer::store_enc_bits(){
    for(int z = 0; z<_frames; z++){
        int offset = z * _vars;
        for(int i=0; i<_vars; i++){
            buf_en_bits[offset+i] = t_in_bits[offset+i];
        }
    }
}

long int CErrorAnalyzer::nb_processed_frames(){
    return nb_analyzed_frames;
}

long int CErrorAnalyzer::nb_fe(){
    return nb_frame_errors;
}

long int CErrorAnalyzer::nb_be(){
    return nb_bit_errors;
}

double CErrorAnalyzer::fer_value(){
    double tFER = (((double)nb_fe())/(nb_processed_frames()));
    return tFER;
}

double CErrorAnalyzer::ber_value(){
    double tBER = (((double)nb_be())/(nb_processed_frames())/(_vars));
    if( _worst_case == true ){
        tBER = (((double)nb_be())/(nb_processed_frames())/(_data));
    }
    return tBER;
}

long int CErrorAnalyzer::nb_be(int data){
    nb_bit_errors = data;
    return nb_bit_errors;
}

long int CErrorAnalyzer::nb_processed_frames(int data){
    nb_analyzed_frames = data;
    return nb_analyzed_frames;
}

long int CErrorAnalyzer::nb_fe(int data){
    nb_frame_errors = data;
    return nb_frame_errors;
}

int CErrorAnalyzer::nb_data(){
    return _data;
}

int CErrorAnalyzer::nb_vars(){
    return _vars;
}

int CErrorAnalyzer::nb_checks(){
    return (_data - _vars);
}
