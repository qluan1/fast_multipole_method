#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#include <time.h>
#include <string.h>
#include "../include/mp_fmm.h"

/* 
    Utilize MPC library https://www.multiprecision.org/mpc/home.html.
*/

void multipoleExpansionReciprocalMP(mpc_t *coeff, mpc_t center,
unsigned long int deg,
unsigned long int numSources, unsigned long int *sourcesIndices,
mpc_t *sources, mpc_t *vec, 
mpc_t *tempVec, unsigned long int tempLen){
    // assert that tempLen >= numSources + 1
    unsigned long int tempIdx = tempLen - 1;
    for (unsigned long int i = 0; i < numSources; i++){
        mpc_set_si(tempVec[i], 1, MPC_RNDNN);     
    }

    for (int k = 1; k <= deg; k++){
        for (int i = 0; i < numSources; i++){
            mpc_mul(tempVec[tempIdx], vec[sourcesIndices[i]], tempVec[i], MPC_RNDNN);
            mpc_add(coeff[k], coeff[k], tempVec[tempIdx], MPC_RNDNN);
            mpc_sub(tempVec[tempIdx], sources[sourcesIndices[i]], center, MPC_RNDNN);
            mpc_mul(tempVec[i], tempVec[i], tempVec[tempIdx], MPC_RNDNN);
        }
    }
    return;
}

void multipoleExpansionLogMP(mpc_t *coeff, mpc_t center, 
unsigned int deg,  
unsigned long int numSources, unsigned long int *sourcesIndices,
mpc_t * sources, mpc_t * vec,
mpc_t * tempVec, unsigned long int tempLen)
{
    // assert that tempLen >= numSources + 1
    unsigned long int tempIdx = tempLen - 1;
    for (int i = 0; i < numSources; i++){
        mpc_add(coeff[0], coeff[0], vec[sourcesIndices[i]], MPC_RNDNN);
        mpc_sub(tempVec[tempIdx], sources[sourcesIndices[i]], center, MPC_RNDNN);
        for (int k = 1; k <= deg; k++){
            mpc_set_si(tempVec[0], -1, MPC_RNDNN);
            mpc_div_ui(tempVec[0], tempVec[0], k, MPC_RNDNN);
            mpc_mul(tempVec[0], tempVec[0], vec[sourcesIndices[i]], MPC_RNDNN);
            mpc_mul(tempVec[0], tempVec[0], tempVec[tempIdx], MPC_RNDNN);
            mpc_add(coeff[k], coeff[k], tempVec[0], MPC_RNDNN);
            mpc_sub(tempVec[0], sources[sourcesIndices[i]], center, MPC_RNDNN);
            mpc_mul(tempVec[tempIdx], tempVec[tempIdx], tempVec[0], MPC_RNDNN);
            //coeff[k] += (-1.0) * vec[sourcesIndices[i]] * temp / (double) k;
            //temp *= sources[sourcesIndices[i]] - center;
        }
    }
    return;
}

void multipoleExpansionMP
(mpc_t *coeff, mpc_t center, 
unsigned int deg,  
unsigned long int numSources, unsigned long int *sourcesIndices,
mpc_t * sources, mpc_t * vec,
mpc_t * tempVec, unsigned long int tempLen, int m)
{
    // m = 0 reciprocal kernel
    // m = 1 log kernel
    if (m == 0){
        multipoleExpansionReciprocalMP(coeff, center, deg, 
                                            numSources, sourcesIndices, 
                                            sources, vec, tempVec, tempLen);
    } else {
        multipoleExpansionLogMP(coeff, center, deg, 
                                            numSources, sourcesIndices, 
                                            sources, vec, tempVec, tempLen);
    }
    return;
}

void shiftMultipoleExpansionReciprocalMP(mpc_t *coeff, 
mpc_t nCenter, mpc_t oCenter, 
unsigned int deg, mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose)
{
    // assert that tempLen >= multipole expansion deg + 2
    // assert that zLen >= multipole expansion deg + 1
    unsigned long int tempIdx = tempLen - 1;
    mpc_set_si(zVec[0], 1, MPC_RNDNN);
    mpc_sub(tempVec[tempIdx], oCenter, nCenter, MPC_RNDNN);
    for (int i = 1; i <= deg; ++i){
        mpc_mul(zVec[i], zVec[i-1], tempVec[tempIdx], MPC_RNDNN);
    }

    mpc_set(tempVec[0], oCoeff[0], MPC_RNDNN);
    for (int l = 1; l <= deg; l ++){
        mpc_set_si(tempVec[l], 0, MPC_RNDNN);
        for (int k = 1; k <= l; k++){
            int stz = mpz_set_str(choose, CHOOSE_STRING[l-1][k-1], 10);
            if (stz == -1){
                printf("error in parsing M choose N\n");
            }
            mpc_mul(tempVec[tempIdx], oCoeff[k], zVec[l-k], MPC_RNDNN);
            mpfr_mul_z(mpc_realref(tempVec[tempIdx]), mpc_realref(tempVec[tempIdx]), choose, MPFR_RNDN);
            mpfr_mul_z(mpc_imagref(tempVec[tempIdx]), mpc_imagref(tempVec[tempIdx]), choose, MPFR_RNDN);
            mpc_add(tempVec[l], tempVec[l], tempVec[tempIdx], MPC_RNDNN);
        }
    }

    for (int i = 0; i <= deg; ++i){
        mpc_add(coeff[i], coeff[i], tempVec[i], MPC_RNDNN);
    }
    return;
}


void shiftMultipoleExpansionLogMP(mpc_t *coeff, 
mpc_t nCenter, mpc_t oCenter, 
unsigned int deg, mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose){
    unsigned long int tempIdx = tempLen - 1;
    mpc_set_si(zVec[0], 1, MPC_RNDNN);
    mpc_sub(tempVec[tempIdx], oCenter, nCenter, MPC_RNDNN);
    for (int i = 1; i <= deg; ++i){
        mpc_mul(zVec[i], zVec[i-1], tempVec[tempIdx], MPC_RNDNN);
    }

    mpc_set(tempVec[0], oCoeff[0], MPC_RNDNN);
    for (int l = 1; l <= deg; l++){
        mpc_mul(tempVec[l], tempVec[0], zVec[l], MPC_RNDNN);
        mpfr_mul_si(mpc_realref(tempVec[l]), mpc_realref(tempVec[l]), -1, MPFR_RNDN);
        mpfr_mul_si(mpc_imagref(tempVec[l]), mpc_imagref(tempVec[l]), -1, MPFR_RNDN);
        mpfr_div_ui(mpc_realref(tempVec[l]), mpc_realref(tempVec[l]), l, MPFR_RNDN);
        mpfr_div_ui(mpc_imagref(tempVec[l]), mpc_imagref(tempVec[l]), l, MPFR_RNDN);
        for (int k = 1; k <= l; k++){
            int stz = mpz_set_str(choose, CHOOSE_STRING[l-1][k-1], 10);
            if (stz == -1){
                printf("error in parsing M choose N\n");
            }
            mpc_mul(tempVec[tempIdx], oCoeff[k], zVec[l-k], MPC_RNDNN);
            mpfr_mul_z(mpc_realref(tempVec[tempIdx]), mpc_realref(tempVec[tempIdx]), choose, MPFR_RNDN);
            mpfr_mul_z(mpc_imagref(tempVec[tempIdx]), mpc_imagref(tempVec[tempIdx]), choose, MPFR_RNDN);
            mpc_add(tempVec[l], tempVec[l], tempVec[tempIdx], MPC_RNDNN);            
        }        
    }

    for (int i = 0; i <= deg; ++i){
        mpc_add(coeff[i], coeff[i], tempVec[i], MPC_RNDNN);
    }
    return;
}

void shiftMultipoleExpansionMP(mpc_t *coeff, 
mpc_t nCenter, mpc_t oCenter, 
unsigned int deg, mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose, int mode){
    if (mode == 0){
        shiftMultipoleExpansionReciprocalMP(coeff, 
        nCenter, oCenter, 
        deg, oCoeff, 
        tempVec, tempLen,
        zVec, zLen,
        choose);
    } else {
        shiftMultipoleExpansionLogMP(coeff, 
        nCenter, oCenter, 
        deg, oCoeff, 
        tempVec, tempLen,
        zVec, zLen,
        choose);
    }
    return;
}

void multipoleToLocalReciprocalMP(mpc_t *coeff,
mpc_t nCenter, mpc_t oCenter,
unsigned long int md, unsigned long int ld,
mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose){

//  tempLen >= md + 2 or ld + 2
//  zLen >= md + 1 or ld + 1
    unsigned int d;
    d = (ld > md) ? ld : md;
    unsigned long int tempIdx = tempLen - 1;

    mpc_sub(tempVec[tempIdx], oCenter, nCenter, MPC_RNDNN);
    mpc_set_si(zVec[0], 1, MPC_RNDNN);
    for (unsigned long int i = 1; i <= d; i++){
        mpc_mul(zVec[i], zVec[i-1], tempVec[tempIdx], MPC_RNDNN);
    }

    mpc_set_si(tempVec[0], 0, MPC_RNDNN);

    int c = -1;
    for (int k = 1; k <= md; k++){
        mpc_mul_si(tempVec[tempIdx], oCoeff[k], c, MPC_RNDNN);
        mpc_div(tempVec[tempIdx], tempVec[tempIdx], zVec[k], MPC_RNDNN);
        mpc_add(tempVec[0], tempVec[0], tempVec[tempIdx], MPC_RNDNN);
        c *= -1;
    }

    for (int l = 1; l <= ld; l++){
        mpc_set_si(tempVec[l], 0, MPC_RNDNN);
        c = -1;
        for (int k = 1; k <= md; k++){
            mpc_mul_si(tempVec[tempIdx], oCoeff[k], c, MPC_RNDNN);
            int stz = mpz_set_str(choose, CHOOSE_STRING[l+k-1][k-1], 10);
            if (stz == -1){
                printf("error in parsing M choose N\n");
            }
            mpfr_mul_z(mpc_realref(tempVec[tempIdx]), mpc_realref(tempVec[tempIdx]), choose, MPFR_RNDN);
            mpfr_mul_z(mpc_imagref(tempVec[tempIdx]), mpc_imagref(tempVec[tempIdx]), choose, MPFR_RNDN);
            mpc_div(tempVec[tempIdx], tempVec[tempIdx], zVec[k], MPC_RNDNN);
            mpc_div(tempVec[tempIdx], tempVec[tempIdx], zVec[l], MPC_RNDNN);
            mpc_add(tempVec[l], tempVec[l], tempVec[tempIdx], MPC_RNDNN);
            c *= -1;
        }
    }

    for (int i = 0; i <= ld; i++){
        mpc_add(coeff[i], coeff[i], tempVec[i], MPC_RNDNN);
    }
    return;
}


void multipoleToLocalLogMP(mpc_t *coeff,
mpc_t nCenter, mpc_t oCenter,
unsigned long int md, unsigned long int ld,
mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose){

    unsigned int d;
    d = (ld > md) ? ld : md;
    unsigned long int tempIdx = tempLen - 1;

    mpc_sub(tempVec[tempIdx], oCenter, nCenter, MPC_RNDNN);
    mpc_set_si(zVec[0], 1, MPC_RNDNN);
    for (unsigned long int i = 1; i <= d; i++){
        mpc_mul(zVec[i], zVec[i-1], tempVec[tempIdx], MPC_RNDNN);
    }

    mpc_mul_si(tempVec[tempIdx], tempVec[tempIdx], -1, MPC_RNDNN);
    mpc_log(tempVec[tempIdx], tempVec[tempIdx], MPC_RNDNN);
    mpc_mul(tempVec[0], oCoeff[0], tempVec[tempIdx], MPC_RNDNN);

    int c = -1;
    for (int k = 1; k <= md; k++){
        mpc_mul_si(tempVec[tempIdx], oCoeff[k], c, MPC_RNDNN);
        mpc_div(tempVec[tempIdx], tempVec[tempIdx], zVec[k], MPC_RNDNN);
        mpc_add(tempVec[0], tempVec[0], tempVec[tempIdx], MPC_RNDNN);
        c *= -1;
    }

    for (int l = 1; l <= ld; l++){
        mpc_mul_si(tempVec[l], oCoeff[0], -1, MPC_RNDNN);
        mpc_div(tempVec[l], tempVec[l], zVec[l], MPC_RNDNN);
        mpc_div_ui(tempVec[l], tempVec[l], l, MPC_RNDNN);
        c = -1;
        for (int k = 1; k <= md; k++){
            mpc_mul_si(tempVec[tempIdx], oCoeff[k], c, MPC_RNDNN);
            int stz = mpz_set_str(choose, CHOOSE_STRING[l+k-1][k-1], 10);
            if (stz == -1){
                printf("error in parsing M choose N\n");
            }
            mpfr_mul_z(mpc_realref(tempVec[tempIdx]), mpc_realref(tempVec[tempIdx]), choose, MPFR_RNDN);
            mpfr_mul_z(mpc_imagref(tempVec[tempIdx]), mpc_imagref(tempVec[tempIdx]), choose, MPFR_RNDN);
            mpc_div(tempVec[tempIdx], tempVec[tempIdx], zVec[k], MPC_RNDNN);
            mpc_div(tempVec[tempIdx], tempVec[tempIdx], zVec[l], MPC_RNDNN);
            mpc_add(tempVec[l], tempVec[l], tempVec[tempIdx], MPC_RNDNN);
            c *= -1;
        }
    }

    for (int i = 0; i <= ld; i++){
        mpc_add(coeff[i], coeff[i], tempVec[i], MPC_RNDNN);
    }
    return;
}


void multipoleToLocalMP(mpc_t *coeff,
mpc_t nCenter, mpc_t oCenter,
unsigned long int md, unsigned long int ld,
mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose, int mode){
    if (mode == 0){
        multipoleToLocalReciprocalMP(coeff, 
        nCenter, oCenter, 
        md, ld,
        oCoeff,
        tempVec, tempLen,
        zVec, zLen,
        choose);
    } else {
        multipoleToLocalLogMP(coeff, 
        nCenter, oCenter, 
        md, ld,
        oCoeff,
        tempVec, tempLen,
        zVec, zLen,
        choose);
    }
    return;
}

void shiftLocalExpansionMP(mpc_t *coeff,
mpc_t nCenter, mpc_t oCenter,
unsigned int d, mpc_t *oCoeff,
mpc_t *tempVec, unsigned long int tempLen, mpc_t shift){
    // tempLen >= ld + 2  (local degree + 2)
    unsigned long int tempIdx = tempLen - 1;
    mpc_sub(shift, oCenter, nCenter, MPC_RNDNN);
    
    for (int i = 0; i <= d; i++){
        mpc_set(tempVec[i], oCoeff[i], MPC_RNDNN);
    }

    for (int j = 0; j <= d; j++){
        for (int k = d - 1 - j; k < d; k++){
            mpc_mul(tempVec[tempIdx], shift, tempVec[k+1], MPC_RNDNN);
            mpc_sub(tempVec[k], tempVec[k], tempVec[tempIdx], MPC_RNDNN);
        }
    }

    for (int i = 0; i <= d; i++){
        mpc_add(coeff[i], coeff[i], tempVec[i], MPC_RNDNN);
    }
    return;
}

void evalLocalExpansionMP(mpc_t res, mpc_t z,
mpc_t center, unsigned int d, mpc_t *coeff, mpc_t shift){
    mpc_set(res, coeff[d], MPC_RNDNN);
    mpc_sub(shift, z, center, MPC_RNDNN);
    for (int k = d-1; k >= 0; k--) {
        mpc_mul(res, res, shift, MPC_RNDNN);
        mpc_add(res, res, coeff[k], MPC_RNDNN);
    }
    return;
}

void evalMultipoleExpansionLogMP(mpc_t res, mpc_t z,
mpc_t center, unsigned int d, mpc_t *coeff){
    // NOTE: center is modified in the process 
    mpc_set(res, coeff[d], MPC_RNDNN);
    mpc_sub(center, z, center, MPC_RNDNN);
    mpc_ui_div(center, 1, center, MPC_RNDNN);
    for (int k = d-1; k > 0; k--){
        mpc_mul(res, res, center, MPC_RNDNN);
        mpc_add(res, res, coeff[k], MPC_RNDNN);
    }
    mpc_mul(res, res, center, MPC_RNDNN);
    mpc_log(center, center, MPC_RNDNN);
    mpc_mul(center, center, coeff[0], MPC_RNDNN);
    mpc_sub(res, res, center, MPC_RNDNN);
    return;
}

void evalMultipoleExpansionReciprocalMP(mpc_t res, mpc_t z,
mpc_t center, unsigned int d, mpc_t *coeff){
    //NOTE: center is modified in the process 
    mpc_set(res, coeff[d], MPC_RNDNN);
    mpc_sub(center, z, center, MPC_RNDNN);
    mpc_ui_div(center, 1, center, MPC_RNDNN);
    for (int k = d-1; k >= 0; k--){
        mpc_mul(res, res, center, MPC_RNDNN);
        mpc_add(res, res, coeff[k], MPC_RNDNN);
    }
    return;
}

void evalMultipoleExpansionMP(mpc_t res, mpc_t z,
mpc_t center, unsigned int d, 
mpc_t *coeff, int mode){
    if (mode == 0){
        evalMultipoleExpansionReciprocalMP(res, z, center, d, coeff);
    } else {
        evalMultipoleExpansionLogMP(res, z, center, d, coeff);
    }
    return;
}


