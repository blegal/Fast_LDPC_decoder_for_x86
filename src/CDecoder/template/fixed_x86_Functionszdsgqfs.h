/*
 *  ldcp_decoder.h
 *  ldpc3
 *
 *  Created by legal on 02/04/11.
 *  Copyright 2011 ENSEIRB. All rights reserved.
 *
 */

/*----------------------------------------------------------------------------*/

//#include "../shared/genere_par_java_simu.h"

#ifndef __x86_Functions__
#define __x86_Functions__

#define llr_from_input(v)  ((2.0 * v)/(sigB * sigB))

static inline int Signe_de_contrib(int entree){
    return ((entree<0)?0:1);
}

static inline int f_abs_fix(int entree){
    int sortie;
    if (entree < 0)
        sortie = -entree;
    else
        sortie = entree;
    return(sortie);
}

static inline int i32b_fix_CondInvSign(int signe, int value){
    if (signe == 1){
        return value;
    }
    return -value;
}

static inline void i32b_fix_min_update(int input, int *i32_min, int *i32_min2)
{
    if(input < *i32_min){
        *i32_min2 = *i32_min;
        *i32_min  = input;
    }else if(input < *i32_min2){
        *i32_min2 = input;
        *i32_min  = *i32_min;
    }
}

//extern int vSAT_NEG_VAR;
//extern int vSAT_POS_VAR;
//static inline int i_contrib_sym_Saturate(int _c){
//    if     (_c < vSAT_NEG_VAR) return vSAT_NEG_VAR;
//    else if(_c > vSAT_POS_VAR) return vSAT_POS_VAR;
//    return _c;
//}


//extern int vSAT_NEG_MSG;
//extern int vSAT_POS_MSG;
//static inline int i_mesg_Saturate(int _c){
//    int  c = _c;
//    if     (_c < vSAT_NEG_MSG) c = vSAT_NEG_MSG;
//    else if(_c > vSAT_POS_MSG) c = vSAT_POS_MSG;
//    return c;
//}


//static inline int i_contrib_sym_fix_sub(int a, int b){
//    int c = a - b;
//    return i_contrib_sym_Saturate( c );
//}

#endif
