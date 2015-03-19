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

#include "CTerminal.h"

void CTerminal::ShowTime(unsigned long secondes)
{
    int ss = secondes % 60;
    int mn = (secondes / 60) % 60;
    int hh = (secondes / 3600);
    printf("%2.2dh%2.2d'%2.2d", hh, mn, ss);
}

CTerminal::CTerminal(CErrorAnalyzer *_counter, CTimer *_timer, double _eb_n0){
    counter   = _counter;
    timer     = _timer;
    Eb_N0     = _eb_n0;
}


void CTerminal::temp_report(){
    if(counter->nb_be() != 0 ){
        double tBER         = counter->ber_value();
        double tFER         = counter->fer_value();
        unsigned long temps = (timer->get_time_sec() < 1) ? 1 : timer->get_time_sec();
        unsigned long fpmn  = (60 * counter->nb_processed_frames()) / temps;
        double        bps   = ((double)fpmn * (double)counter->nb_data()) / 60.0 / 1000.0 / 1000.0;
        printf("(RT) FRA: %8ld | FE: %3d | FER: %2.2e | BE : %5d | BER: %2.2e | [BE/FE] : %4f | FPM: %3ld | BPS: %2.2f | ETA: ",
                counter->nb_processed_frames(),
                (int)counter->nb_fe(), tFER,
                (int)counter->nb_be(), tBER,
                (float)counter->nb_be() / (float)counter->nb_fe(),
                fpmn, bps);
//        printf("(RT) FRA: %8ld | FE: %3d | BER: %2.2e | FPM: %3ld | BPS: %2.2f | ETA: ", counter->nb_processed_frames(), (int)counter->nb_fe(), tBER, fpmn, bps);
        ShowTime( temps );
        int eta  = (int)(((double)temps / (double)counter->nb_fe()) * (counter->fe_limit()));
        printf(" | ETR: ");
        ShowTime( eta );
        printf("\r");
    }else{
        double tBER         = ( 1.0) / ((double)counter->nb_processed_frames()) / counter->nb_vars();
        double tFER         = ( 1.0) / ((double)counter->nb_processed_frames());
        unsigned long temps = (timer->get_time_sec() < 1) ? 1 : timer->get_time_sec();
        unsigned long fpmn  = (60 * counter->nb_processed_frames()) / temps;
        double        bps   = ((double)fpmn * (double)counter->nb_data()) / 60.0 / 1000.0 / 1000.0;
        printf("(RT) FRA: %8ld | FE: %3d | FER: %2.2e | BE : %5d | BER: %2.2e | [BE/FE] : %4f | FPM: %3ld | BPS: %2.2f | ETA: ",
                counter->nb_processed_frames(),
                (int)counter->nb_fe(), tFER,
                (int)counter->nb_be(), tBER,
                (float)counter->nb_be() / 1.0,
                fpmn, bps);
//        printf("(RT) FRA: %8ld | FE: %3d | BER: <%2.2e | FPM: %3ld | BPS: %2.2f | ETA: ", counter->nb_processed_frames(), (int)counter->nb_fe(), tBER, fpmn, bps);
        ShowTime( temps );
        printf(" | ETR: INF.");
        printf("\r");
    }
    fflush(stdout);
}


void CTerminal::final_report(){
    double tBER = counter->ber_value();
    double tFER = counter->fer_value();
    unsigned long temps = timer->get_time_sec() + 1;
    unsigned long fpmn  = (60 * counter->nb_processed_frames()) / temps;
    double        bps   = ((double)fpmn * (double)counter->nb_data()) / 60.0 / 1000.0 / 1000.0;
    double be_par_fe = (double) counter->nb_be() / (double) counter->nb_fe();
    printf("SNR = %.2f | BER =  %2.3e | FER =  %2.3e | BPS =  %2.2f | MATRICES = %10ld| FE = %d | BE = %d | BE/FE = %.1f | RUNTIME = ", Eb_N0, tBER, tFER, bps, counter->nb_processed_frames(), (int) counter->nb_fe(), (int) counter->nb_be(), be_par_fe);
//    printf("SNR = %.2f | BER =  %2.3e | FER =  %2.3e | BPS =  %2.2f | MATRICES = %10ld| FE = %d | BE = %d | RUNTIME = ", Eb_N0, tBER, tFER, bps, counter->nb_processed_frames(), (int)counter->nb_fe(), (int)counter->nb_be());
    ShowTime( temps );
    printf("\n");
    fflush(stdout);
}
