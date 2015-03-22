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
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <limits.h>

using namespace std;

#include <omp.h>

#define pi  3.1415926536
#include "./CTimer/CTimer.h"
#include "./CTrame/CTrame.h"
#include "CChanel/ChanelLibrary.h"
#include "./CEncoder/EncoderLibrary.h"
#include "./CErrorAnalyzer/CErrorAnalyzer.h"
#include "./CTerminal/CTerminal.h"
#include "./Constantes/constantes_sse.h"
#include "./CFixPointConversion/CFastFixConversion.h"
#include "./CChanel/ChanelLibrary.h"
#include "./CTools/CTools.h"

struct param_decoder {
    float oms_offset_float;
    int   oms_offset_fixed;

    float nms_factor_float;
    int   nms_factor_fixed;

    int  nb_iters;
};

struct param_simulation {
    double snr_min;
    double snr_max;
    double snr_pas;

    int    fe_limit;

    int    channel_type;
    bool   qpsk_channel;
    bool   Es_N0;
    bool   norm_channel;

    bool   ber_limit;
    double ber_limit_value;

    bool   fer_limit;
    double fer_limit_value;

    bool   llr_optimization;
    bool   real_encoder;
    int    STOP_TIMER_SECOND;
    int    nb_frames;
    bool   worst_case_fer;
    bool   show_llr_histo;
};

#include "./CDecoder/DecoderLibrary.h"

#define NOEUD       _N
#define PARITE      _K
#define MESSAGE     _M
#define INFORMATION NmoinsK

//
// VARIABLES UTILISEE AFIN DE CALCULER LES VALEURS DES SATURATIONS
// DES MESSAGES ET DES VARIABLES
//
int BITS_LLR     = NB_BITS_MESSAGES;
int BITS_VAR     = NB_BITS_VARIABLES;
int BITS_MSG     = NB_BITS_MESSAGES;

int vSAT_NEG_MSG = SAT_NEG_MSG;
int vSAT_POS_MSG = SAT_POS_MSG;

int vSAT_NEG_VAR = SAT_NEG_VAR;
int vSAT_POS_VAR = SAT_POS_VAR;

int vSAT_NEG_LLR  = SAT_NEG_MSG;
int vSAT_POS_LLR  = SAT_POS_MSG;
int vFRAQ_LLR     = NB_BITS_VARIABLES / 2;

int FACTEUR_BETA = (0x0001<<(NB_BITS_MESSAGES/2));

#define MAX_THREADS 4

////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
    srand(0);
    printf("(II) LDPC DECODER - Scheduled Fixed-point decoder (SSE version - 16 frames)\n");
    printf("(II) MANIPULATION DE DONNEES (Fixed Point, %ld bits)\n", 8 * sizeof (int));
    printf("(II) GENEREE : %s - %s\n", __DATE__, __TIME__);

    param_simulation p_simulation;
    p_simulation.snr_min      = 0.5;
    p_simulation.snr_max      = 3.0;
    p_simulation.snr_pas      = 0.1;
    p_simulation.fe_limit     = 100;
    p_simulation.channel_type = 3;
    p_simulation.norm_channel = false;

    p_simulation.llr_optimization = false;
    p_simulation.real_encoder     = false;
    p_simulation.qpsk_channel     = false;
    p_simulation.Es_N0            = false;
    p_simulation.worst_case_fer   = false;
    p_simulation.show_llr_histo   = false;

    param_decoder p_decoder;
    p_decoder.nb_iters   = 30;

    p_decoder.nms_factor_fixed = 29;
    p_decoder.nms_factor_float = 0.75;

    p_decoder.oms_offset_fixed = 1;
    p_decoder.oms_offset_float = 0.15;


    int    NOMBRE_ITERATIONS = 30;
    int    STOP_TIMER_SECOND = -1;
    int    nb_frames         = 16;

    string format = "float";
    string type = "OMS";
    string arch = "sse";

    //
    // ON SELECTIONNE LE FORMAT DE DONNEES DU DECODEUR
    //
    if (strcmp(argv[1], "-fixed") == 0) {
        format = "fixed";
    } else if (strcmp(argv[1], "-float") == 0) {
        format = "float";
    } else {
        printf("(EE) First argument must be the decoder data format [-float, -fixed]\n");
        printf("(EE) -float : floatting   point format decoder\n");
        printf("(EE) -fixed : fixed-point point format decoder\n");
        exit(0);
    }

    //
    // ON CONFIGURE LE NOMBRE DE THREAD A UTILISER PAR DEFAUT
    //
    int NUM_ACTIVE_THREADS = 1;
    omp_set_num_threads(NUM_ACTIVE_THREADS);

    //
    // ON VA PARSER LES ARGUMENTS DE LIGNE DE COMMANDE
    //
    for (int p = 2; p < argc; p++) {

        //
        // REGLAGE DES PARAMETRES DE SIMULATION
        //
        if (strcmp(argv[p], "-min") == 0) {
            p_simulation.snr_min = atof(argv[p + 1]);
            p += 1;

        } else if (strcmp(argv[p], "-max") == 0) {
            p_simulation.snr_max = atof(argv[p + 1]);
            p += 1;

        } else if (strcmp(argv[p], "-pas") == 0) {
            p_simulation.snr_pas = atof(argv[p + 1]);
            p += 1;

        } else if (strcmp(argv[p], "-fer") == 0) {
            p_simulation.fe_limit = atoi(argv[p + 1]);
            p += 1;

        } else if (strcmp(argv[p], "-wc_fer") == 0) {
            p_simulation.worst_case_fer = true;

        } else if (strcmp(argv[p], "-histo") == 0) {
            p_simulation.show_llr_histo = true;

        //
        // LIMITATION DU TEMPS DE SIMULTATION
        //
        } else if (strcmp(argv[p], "-timer") == 0) {
            STOP_TIMER_SECOND = atoi(argv[p + 1]);
            p += 1;

        } else if (strcmp(argv[p], "-qef") == 0) {
            p_simulation.ber_limit       = true;
            p_simulation.ber_limit_value = (atof(argv[p + 1]));
            p += 1;


        } else if (strcmp(argv[p], "-tfer") == 0) {
            p_simulation.fer_limit       = true;
            p_simulation.fer_limit_value = (atof(argv[p + 1]));
            p += 1;


            //
            // REGLAGE DES DU MODELE DU CANAL
            //
        } else if (strcmp(argv[p], "-bpsk") == 0) {
            p_simulation.qpsk_channel = false;

        } else if (strcmp(argv[p], "-qpsk") == 0) {
            p_simulation.qpsk_channel = true;

        } else if (strcmp(argv[p], "-Eb/N0") == 0) {
            p_simulation.Es_N0 = false;

        } else if (strcmp(argv[p], "-Es/N0") == 0) {
            p_simulation.Es_N0 = true;

        } else if (strcmp(argv[p], "-norm-channel") == 0) {
            p_simulation.norm_channel = true;

        } else if (strcmp(argv[p], "-encoder") == 0) {
            p_simulation.real_encoder = true;


            //
            // REGLAGE DES DU MODELE DU CANAL
            //
        } else if (strcmp(argv[p], "-thread") == 0) {
            int nThreads = atoi(argv[p + 1]);
            if (nThreads > 4) {
                printf("(WW) Number of thread can be higher than 4 => Using 4 threads.");
                NUM_ACTIVE_THREADS = 4;
            } else if (nThreads < 1) {
                printf("(WW) Number of thread can be lower than 1 => Using 1 thread.");
                NUM_ACTIVE_THREADS = 1;
            } else {
                NUM_ACTIVE_THREADS = nThreads;
            }
            omp_set_num_threads(NUM_ACTIVE_THREADS);
            p += 1;

        //
        // INITIALISATION ALEATOIRE DU GENERATEUR ALEATOIRE
        //
        } else if (strcmp(argv[p], "-random") == 0) {
            printf("(II) Random Generator REAL initialization\n");
            srand(time(NULL));


            //
            // REGLAGE DES DU MODELE DU CANAL
            //
        } else if (strcmp(argv[p], "-sse") == 0) {
            arch = "sse";
            if (format.compare("float") == 0)
                nb_frames = 4;
            else
                nb_frames = 16;

        } else if (strcmp(argv[p], "-avx") == 0) {
            arch = "avx";
            if (format.compare("float") == 0)
                nb_frames = 8;
            else
                nb_frames = 32;

        } else if (strcmp(argv[p], "-x86") == 0) {
            arch = "x86";
            nb_frames = 1;

        } else if (strcmp(argv[p], "-NMS") == 0) {
            type = "NMS";
            p_decoder.nms_factor_float = atof(argv[p + 1]);
            p_decoder.nms_factor_fixed = atoi(argv[p + 1]);
            p += 1;

        } else if (strcmp(argv[p], "-OMS") == 0) {
            type = "OMS";
            p_decoder.oms_offset_float = atof(argv[p + 1]);
            p_decoder.oms_offset_fixed = atoi(argv[p + 1]);
            p += 1;

        } else if (strcmp(argv[p], "-iter") == 0) {
            NOMBRE_ITERATIONS      = atoi(argv[p + 1]);
            p_decoder.nb_iters = atoi(argv[p + 1]);
            p += 1;

        } else if (strcmp(argv[p], "-var") == 0) {
            vSAT_NEG_VAR = (-(0x0001 << (atoi(argv[p + 1]) - 1)) + 1);
            vSAT_POS_VAR = ( (0x0001 << (atoi(argv[p + 1]) - 1)) - 1);
            BITS_VAR     = atoi(argv[p + 1]);
            p           += 1;

        } else if (strcmp(argv[p], "-msg") == 0) {
            vSAT_NEG_MSG = (-(0x0001 << (atoi(argv[p + 1]) - 1)) + 1);
            vSAT_POS_MSG = ( (0x0001 << (atoi(argv[p + 1]) - 1)) - 1);
            BITS_MSG     = atoi(argv[p + 1]);
            p           += 1;

        } else if (strcmp(argv[p], "-llr") == 0) {
            vSAT_NEG_LLR = (-(0x0001 << (atoi(argv[p + 1]) - 1)) + 1);
            vSAT_POS_LLR = ((0x0001 << (atoi(argv[p + 1]) - 1)) - 1);
            BITS_LLR     = atoi(argv[p + 1]);
            vFRAQ_LLR    = BITS_LLR / 2;
            FACTEUR_BETA = (0x0001 << (vFRAQ_LLR));
            p           += 1;

        } else if (strcmp(argv[p], "-fraq") == 0) {
            vFRAQ_LLR    = atoi(argv[p + 1]);
            FACTEUR_BETA = (0x0001 << (vFRAQ_LLR));
            p           += 1;

        } else {
            printf("(EE) Unknown argument (%d) => [%s]\n", p, argv[p]);
            exit(0);
        }
    }

    double rendement = (float) (INFORMATION) / (float) (NOEUD);
    printf("(II) NUMBER OF // THREAD  : %d\n", NUM_ACTIVE_THREADS);
    printf("(II) Code LDPC (N, K)     : (%d,%d)\n", NOEUD, PARITE);
    printf("(II) Rendement du code    : %.3f\n", rendement);
    if (type == "MS") {
        printf("(II) LDPC HEURISTIC (CN)  : MIN-SUM\n");
    } else if (type == "OMS") {
        printf("(II) LDPC HEURISTIC (CN)  : OFFSET MIN-SUM\n");
    } else if (type == "NMS") {
        printf("(II) LDPC HEURISTIC (CN)  : NORMALIZED-MIN-SUM\n");
    }
    printf("(II) # ITERATIONs du CODE : %d\n", NOMBRE_ITERATIONS);
    printf("(II) FER LIMIT FOR SIMU   : %d\n", p_simulation.fe_limit);
    printf("(II) SIMULATION  RANGE    : [%.2f, %.2f], STEP = %.2f\n", p_simulation.snr_min, p_simulation.snr_max, p_simulation.snr_pas);

    if (format.compare("fixed") == 0) {
        printf("(II) LLR DATA    Q(%d,%d)   : %d bits [%d, %d]\n", (BITS_LLR - vFRAQ_LLR), (vFRAQ_LLR), BITS_LLR, vSAT_NEG_LLR, vSAT_POS_LLR);
        printf("(II) MESSAGE     Q(%d,%d)   : %d bits [%d, %d]\n", (BITS_MSG - vFRAQ_LLR), (vFRAQ_LLR), BITS_MSG, vSAT_NEG_MSG, vSAT_POS_MSG);
        printf("(II) VARIABLE    Q(%d,%d)   : %d bits [%d, %d]\n", (BITS_VAR - vFRAQ_LLR), (vFRAQ_LLR), BITS_VAR, vSAT_NEG_VAR, vSAT_POS_VAR);

        if (type.compare("OMS") == 0)
            printf("(II) VARIABLE BETA FIX    : %d\n", p_decoder.oms_offset_fixed);
        if (type.compare("NMS") == 0)
            printf("(II) NORMALIZE FACTOR     : %d\n", p_decoder.nms_factor_fixed);

        printf("(II) FACTEUR BETA (LLR)   : %d\n", FACTEUR_BETA);
    } else {
        if (type.compare("OMS") == 0)
            printf("(II) OFFSET FACTOR        : %f\n", p_decoder.oms_offset_float);
        if (type.compare("NMS") == 0)
            printf("(II) NORMALIZE FACTOR     : %f\n", p_decoder.nms_factor_float);
    }


    CTimer simu_timer(true);

    //
    // ALLOCATION DYNAMIQUE DES DONNESS NECESSAIRES A LA SIMULATION DU SYSTEME
    //
    CTrame* simu_data[MAX_THREADS];
    for(int i=0; i<4; i++){
        simu_data[i] = new CTrame(NOEUD, PARITE, nb_frames);
    }

    CDecoder* decoder[MAX_THREADS];
    for(int i=0; i<4; i++){
        decoder[i] = CreateDecoder(type, arch, format, p_decoder, vSAT_NEG_VAR, vSAT_POS_VAR, vSAT_NEG_MSG, vSAT_POS_MSG/*, msOffset, msFactor, OFFSET_FACTOR, NORMALIZED_FACTOR*/);
    }

    Encoder *encoder[MAX_THREADS];
    for(int i=0; i<4; i++){
        encoder[i] = EncoderLibrary(p_simulation.real_encoder, simu_data[i]);
    }

    CChanel* noise[MAX_THREADS];
    for(int i=0; i<4; i++){
        noise[i] = CreateChannel(arch, p_simulation.channel_type, simu_data[i], p_simulation.qpsk_channel, p_simulation.Es_N0);
        noise[i]->setNormalize( p_simulation.norm_channel );
    }

    //
    // ON CREE L'OBJET EN CHARGE DE LA CONVERSION EN VIRGULE FIXE DE L'INFORMATION DU CANAL
    //
    CFixConversion* conv_fp   [MAX_THREADS];
    CErrorAnalyzer* errCounter[MAX_THREADS];

    double Eb_N0 = p_simulation.snr_min;
    while (Eb_N0 <= p_simulation.snr_max) {

        //
        // ON CREE LE CANAL DE COMMUNICATION (BRUIT GAUSSIEN)
        //

        for(int i=0; i<4; i++){
            noise[i]->configure(Eb_N0);
        }

        for(int i=0; i<4; i++){
            decoder[i]->setSigmaChannel(noise[i]->get_SigB());
        }

        for (int i = 0; i < 4; i++) {
            conv_fp[i] = new CFastFixConversion(simu_data[i], FACTEUR_BETA, vSAT_NEG_LLR, vSAT_POS_LLR);
            conv_fp[i]->ShowHistoOnDestroy(p_simulation.show_llr_histo);
        }

        bool auto_fe_mode = false;
        CErrorAnalyzer  errCounters  (simu_data[0], p_simulation.fe_limit, auto_fe_mode, p_simulation.worst_case_fer);
        for(int i=0; i<4; i++){
            errCounter[i] = new CErrorAnalyzer(simu_data[i], p_simulation.fe_limit, auto_fe_mode, p_simulation.worst_case_fer);
        }

        // ON GENERE LA PREMIERE TRAME BRUITEE
        for(int i=0; i<4; i++){
            encoder[i]->encode();
        }

        for(int i=0; i<4; i++){
            noise[i]->generate();
        }

        for(int i=0; i<4; i++){
            conv_fp[i]->generate();
        }

        for(int i=0; i<4; i++){
            errCounter[i]->store_enc_bits();
        }

        //
        // ON CREE UN OBJET POUR LA MESURE DU TEMPS DE SIMULATION (REMISE A ZERO POUR CHAQUE Eb/N0)
        //
        CTimer temps_ecoule(true);

        //
        // ON CREE L'OBJET EN CHARGE DES INFORMATIONS DANS LE TERMINAL UTILISATEUR
        //
        CTerminal terminal(&errCounters, &temps_ecoule, Eb_N0);

        CTimer timer[MAX_THREADS];
        long int etime[MAX_THREADS] = {0, 0, 0, 0};

        CTimer biais_measure(false);
        biais_measure.start();
        biais_measure.stop();
        unsigned int biais = biais_measure.get_time_us();

        while (1) {
            const int maxLoopF = 32768;
            int loopf  = (8 * NUM_ACTIVE_THREADS) * (64800 / NOEUD);
            loopf      = loopf > maxLoopF ? maxLoopF: loopf;

            int d1[maxLoopF], d2[maxLoopF], d3[maxLoopF], d4[maxLoopF];
            int f1[maxLoopF], f2[maxLoopF], f3[maxLoopF], f4[maxLoopF];

            #pragma omp parallel sections
            {

                #pragma omp section
                {
                    for (int q = 0; q < loopf; q++) {
                        float *f_llr = simu_data[0]->get_t_noise_data(); // [NOEUD];
                        char *i_llr   = simu_data[0]->get_t_fpoint_data();  // [NOEUD];
                        char *o_llr   = simu_data[0]->get_t_decode_data();  // [NOEUD];

                        timer[0].start();
                        decoder[0]->decode(f_llr, o_llr, NOMBRE_ITERATIONS);
                        decoder[0]->decode(i_llr, o_llr, NOMBRE_ITERATIONS);
                        timer[0].stop();
                        etime[0] += (timer[0].get_time_us() - biais);
//                        printf("show = %d\n", timer[0].get_time_us());

                        encoder[0]->encode();
                        noise[0]->generate();  // ON GENERE LE BRUIT DU CANAL
                        conv_fp[0]->generate(); // ON CONVERTIT LES DONNEES EN VIRGULE FIXE
                        int q1 = errCounter[0]->nb_be();
                        int fr = errCounter[0]->nb_fe();
                        errCounter[0]->generate();
                        d1[q] = errCounter[0]->nb_be() - q1;
                        f1[q] = errCounter[0]->nb_fe() - fr;
                        errCounter[0]->store_enc_bits();
                    }
                }

                #pragma omp section
                {
                    for (int q = 0; q < loopf; q++) {
                        float *f_llr = simu_data[1]->get_t_noise_data(); // [NOEUD];
                        char *i_llr  = simu_data[1]->get_t_fpoint_data();  // [NOEUD];
                        char *o_llr  = simu_data[1]->get_t_decode_data();  // [NOEUD];

                        timer[1].start();
                        decoder[1]->decode(f_llr, o_llr, NOMBRE_ITERATIONS);
                        decoder[1]->decode(i_llr, o_llr, NOMBRE_ITERATIONS);
                        timer[1].stop();
                        etime[1] += timer[1].get_time_us() - biais;

                        encoder[1]->encode();
                        noise[1]->generate();  // ON GENERE LE BRUIT DU CANAL
                        conv_fp[1]->generate(); // ON CONVERTIT LES DONNEES EN VIRGULE FIXE
                        int fr = errCounter[1]->nb_fe();
                        int q2 = errCounter[1]->nb_be();
                        errCounter[1]->generate();
                        d2[q] = errCounter[1]->nb_be() - q2;
                        f2[q] = errCounter[1]->nb_fe() - fr;
                        errCounter[1]->store_enc_bits();
                    }
                }

                #pragma omp section
                {
                    for (int q = 0; q < loopf; q++) {
                        float *f_llr = simu_data[2]->get_t_noise_data();   // [NOEUD];
                        char *i_llr  = simu_data[2]->get_t_fpoint_data();  // [NOEUD];
                        char *o_llr  = simu_data[2]->get_t_decode_data();  // [NOEUD];

                        timer[2].start();
                        decoder[2]->decode(f_llr, o_llr, NOMBRE_ITERATIONS);
                        decoder[2]->decode(i_llr, o_llr, NOMBRE_ITERATIONS);
                        timer[2].stop();
                        etime[2] += timer[2].get_time_us() - biais;

                        encoder[2]->encode();
                        noise[2]->generate();   // ON GENERE LE BRUIT DU CANAL
                        conv_fp[2]->generate(); // ON CONVERTIT LES DONNEES EN VIRGULE FIXE
                        int q3 = errCounter[2]->nb_be();
                        int fr = errCounter[2]->nb_fe();
                        errCounter[2]->generate();
                        d3[q] = errCounter[2]->nb_be() - q3;
                        f3[q] = errCounter[2]->nb_fe() - fr;
                        errCounter[2]->store_enc_bits();
                    }
                }

                #pragma omp section
                {
                    for (int q = 0; q < loopf; q++) {
                        float *f_llr = simu_data[3]->get_t_noise_data(); // [NOEUD];
                        char *i_llr = simu_data[3]->get_t_fpoint_data();  // [NOEUD];
                        char *o_llr = simu_data[3]->get_t_decode_data();  // [NOEUD];

                        timer[3].start();
                        decoder[3]->decode(f_llr, o_llr, NOMBRE_ITERATIONS);
                        decoder[3]->decode(i_llr, o_llr, NOMBRE_ITERATIONS);
                        timer[3].stop();
                        etime[3] += timer[3].get_time_us() - biais;

                        encoder[3]->encode();
                        noise[3]->generate();  // ON GENERE LE BRUIT DU CANAL
                        conv_fp[3]->generate(); // ON CONVERTIT LES DONNEES EN VIRGULE FIXE
                        int q4 = errCounter[3]->nb_be();
                        int fr = errCounter[3]->nb_fe();
                        errCounter[3]->generate();
                        d4[q] = errCounter[3]->nb_be() - q4;
                        f4[q] = errCounter[3]->nb_fe() - fr;
                        errCounter[3]->store_enc_bits();
                    }
                }
            }
            //
            // ON COMPTE LE NOMBRE D'ERREURS DANS LA TRAME DECODE
            //
            for (int q = 0; q < loopf; q++) {

                int diff = ((f1[q] - 1) > 0) ? (f1[q] - 1) : 0;
                errCounters.generate(d1[q] - diff);
                for (int z = 1; z < nb_frames; z++) errCounters.generate(f1[q] > z ? 1 : 0);

                diff = ((f2[q] - 1) > 0) ? (f2[q] - 1) : 0;
                errCounters.generate(d2[q] - diff);
                for (int z = 1; z < nb_frames; z++) errCounters.generate(f2[q] > z ? 1 : 0);

                diff = ((f3[q] - 1) > 0) ? (f3[q] - 1) : 0;
                errCounters.generate(d3[q] - diff);
                for (int z = 1; z < nb_frames; z++) errCounters.generate(f3[q] > z ? 1 : 0);

                diff = ((f4[q] - 1) > 0) ? (f4[q] - 1) : 0;
                errCounters.generate(d4[q] - diff);
                for (int z = 1; z < nb_frames; z++) errCounters.generate(f4[q] > z ? 1 : 0);
            }
            //
            // ON compare le Frame Error avec la limite imposee par l'utilisateur. Si on depasse
            // alors on affiche les resultats sur Eb/N0 courant.
            //
            if (errCounters.fe_limit_achieved() == true) {
                break;
            }

            //
            // ON REGARDE SI L'UTILISATEUR A LIMITE LE TEMPS DE SIMULATION...
            //
            if ((simu_timer.get_time_sec() >= STOP_TIMER_SECOND) && (STOP_TIMER_SECOND != -1)) {
                break;
            }

            //
            // AFFICHAGE A L'ECRAN DE L'EVOLUTION DE LA SIMULATION SI NECESSAIRE
            //
            terminal.temp_report();
        }


        terminal.final_report();

        if (STOP_TIMER_SECOND != -1) {
            printf("(PERF) H. LAYERED %d fixed, %dx%d LDPC code, %d its, %d threads\n", nb_frames, NOEUD, PARITE, NOMBRE_ITERATIONS, NUM_ACTIVE_THREADS);
            float sum = 0.0;
            for (int z = 0; z < NUM_ACTIVE_THREADS; z++) {
                float nf = (errCounters.nb_processed_frames() / 4); // 4 car 4 threads...
                float nb = ((nf) * (1000000.0 / etime[z]) * NOEUD) / 1000.0 / 1000.0;
                printf("(PERF) Kernel Execution time = %ld us for %.0f frames => %1.3f Mbps\n", etime[z], nf, nb);
                sum += nb;
            }
            float latenc1 = etime[0] * nb_frames / (errCounters.nb_processed_frames()/4); // en us
            printf("(PERF) SNR = %.2f, ITERS = %d, LATENCY    = %1.3f us\n", Eb_N0, NOMBRE_ITERATIONS, latenc1);
            printf("(PERF) SNR = %.2f, ITERS = %d, THROUGHPUT = %1.3f Mbps\n", Eb_N0, NOMBRE_ITERATIONS, sum);
            printf("(PERF) Total Kernel throughput = %1.3f Mbps\n", sum);
        }

        Eb_N0 = Eb_N0 + p_simulation.snr_pas;

        if ((simu_timer.get_time_sec() >= STOP_TIMER_SECOND) && (STOP_TIMER_SECOND != -1)) {
            printf("(II) THE SIMULATION HAS STOP DUE TO THE (USER) TIME CONTRAINT.\n");
            break;
        }

        if (p_simulation.ber_limit == true) {
            if (errCounters.ber_value() < p_simulation.ber_limit_value) {
                printf("(II) THE SIMULATION HAS STOP DUE TO THE (USER) QUASI-ERROR FREE CONTRAINT (on BER).\n");
                break;
            }
        }

        if (p_simulation.fer_limit == true) {
            if (errCounters.fer_value() < p_simulation.fer_limit_value) {
                printf("(II) THE SIMULATION HAS STOP DUE TO THE (USER) QUASI-ERROR FREE CONTRAINT (on FER).\n");
                break;
            }
        }

    }

    ////////////////////////////////////////////////////////////////////////////////
    //
    //
    // SECOND EVALUATION OF THE THROUGHPUT WITHOUT ENCODED FRAME REGENERATION
    //
    //
    if( STOP_TIMER_SECOND != -1 )
    {
        int exec = 0;
        const int t_eval = STOP_TIMER_SECOND;


        //
        // ONE THREAD MODE
        //
        if (NUM_ACTIVE_THREADS == 1) {
            CTimer t_Timer1(true);
            while (t_Timer1.get_time_sec() < t_eval) {
                for (int qq = 0; qq < 256; qq++) {
                    // to limit timer runtime impact on performances (for very small LDPC codes)
                    // Indeed, depending on OS and CTimer implementations, time read can be long...
                    decoder[0]->decode(simu_data[0]->get_t_fpoint_data(), simu_data[0]->get_t_decode_data(), NOMBRE_ITERATIONS);
                    exec += 1;
                }
            }
            t_Timer1.stop();
            float debit = _N * ((exec * nb_frames) / ((float) t_Timer1.get_time_sec()));
            debit /= 1000000.0f;
            printf("(PERF1) Total Kernel throughput = %1.3f Mbps\n", debit);
        }

        //
        // TWO THREAD MODE
        //
        if (NUM_ACTIVE_THREADS == 2) {
            exec = 0;
            omp_set_num_threads(2);
            CTimer t_Timer2(true);

            while (t_Timer2.get_time_sec() < t_eval) {
                const int looper = 256;
                #pragma omp parallel sections
                {
                    #pragma omp section
                    {
                        for (int qq = 0; qq < looper; qq++)
                            decoder[0]->decode(simu_data[0]->get_t_fpoint_data(), simu_data[1]->get_t_decode_data(), NOMBRE_ITERATIONS);
                    }
                    #pragma omp section
                    {
                        for (int qq = 0; qq < looper; qq++)
                            decoder[1]->decode(simu_data[1]->get_t_fpoint_data(), simu_data[2]->get_t_decode_data(), NOMBRE_ITERATIONS);
                    }
                }
                exec += 2 * looper;
            }
            t_Timer2.stop();

            // for each decoder run, we decoded nb_frames codewords (depending on the SIMD width)
            float debit = _N * ((exec * nb_frames) / ((float) t_Timer2.get_time_sec()));
            debit /= 1000000.0f;
            printf("(PERF2) Total Kernel throughput = %1.3f Mbps\n", debit);
        }


        //
        // FOUR THREAD MODE
        //
        if (NUM_ACTIVE_THREADS == 4) {
            exec = 0;
            omp_set_num_threads(4);
            CTimer t_Timer3(true);

            while (t_Timer3.get_time_sec() < t_eval) {
                const int looper = 256;
                #pragma omp parallel sections
                {
                    #pragma omp section
                    {
                        for (int qq = 0; qq < looper; qq++)
                            decoder[0]->decode(simu_data[0]->get_t_fpoint_data(), simu_data[1]->get_t_decode_data(), NOMBRE_ITERATIONS);
                    }
                    #pragma omp section
                    {
                        for (int qq = 0; qq < looper; qq++)
                            decoder[1]->decode(simu_data[1]->get_t_fpoint_data(), simu_data[2]->get_t_decode_data(), NOMBRE_ITERATIONS);
                    }
                    #pragma omp section
                        {
                        for (int qq = 0; qq < looper; qq++)
                            decoder[2]->decode(simu_data[2]->get_t_fpoint_data(), simu_data[3]->get_t_decode_data(), NOMBRE_ITERATIONS);
                    }
                    #pragma omp section
                    {
                        for (int qq = 0; qq < looper; qq++)
                            decoder[3]->decode(simu_data[3]->get_t_fpoint_data(), simu_data[4]->get_t_decode_data(), NOMBRE_ITERATIONS);
                    }
                }
                exec += 4 * looper;
            }
            t_Timer3.stop();

            float debit = _N * ((exec * nb_frames) / ((float) t_Timer3.get_time_sec()));
            debit /= 1000000.0f;
            printf("(PERF4) Total Kernel throughput = %1.3f Mbps\n", debit);
        }
        exit(0);
    }
    //
    //
    //
    ////////////////////////////////////////////////////////////////////////////////

    //
    // ON FAIT LE MENAGE PARMIS TOUS LES OBJETS CREES DYNAMIQUEMENT...
    //
    for(int i=0; i<4; i++){
        delete simu_data[i];
        delete noise[i];
        delete decoder[i];
        delete encoder[i];
        delete errCounter[i];
        delete conv_fp[i];
    }

    return 1;
}
