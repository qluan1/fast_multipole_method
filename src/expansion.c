#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "../include/mp_fmm.h"


/* 
Expansion based on double precision computation
*/

void multipoleExpansionReciprocalDP(double complex *coeff, double complex center, 
unsigned int deg, 
unsigned long int numSources, unsigned long int *sourcesIndices, 
double complex * sources, double complex *vec, double complex *tempS){

    for (int i = 0; i < numSources; i++){
        tempS[i] = 1.0 + 0.0*I;
    }

    for (int k = 1; k <= deg; k++){
        for (int i = 0; i < numSources; i++){
            coeff[k] += vec[sourcesIndices[i]]*tempS[i];
            tempS[i] *= sources[sourcesIndices[i]] - center;
        }
    }
    return;
}


void multipoleExpansionLogDP(double complex *coeff, double complex center, 
unsigned int deg,  
unsigned long int numSources, unsigned long int *sourcesIndices,
double complex * sources, double complex * vec){

    double complex temp;
    for (int i = 0; i < numSources; i++){
        coeff[0] += vec[sourcesIndices[i]];
        temp = sources[sourcesIndices[i]] - center;
        for (int k = 1; k <= deg; k++){
            coeff[k] += (-1.0) * vec[sourcesIndices[i]] * temp / (double) k;
            temp *= sources[sourcesIndices[i]] - center;
        }
    }
    return;
}



void multipoleExpansionDP(double complex *coeff, double complex center, unsigned int deg,  
unsigned long int numSources, unsigned long int *sourcesIndices,
double complex * sources, double complex * vec, int m, double complex *tempS){
    if (m == 0){
        multipoleExpansionReciprocalDP(coeff, center, deg, 
                                            numSources, sourcesIndices, 
                                            sources, vec, tempS);
    } else {
        multipoleExpansionLogDP(coeff, center, deg, 
                                            numSources, sourcesIndices, 
                                            sources, vec);
    }
    return;
}


void shiftMultipoleExpansionReciprocalDP(double complex* coeff, double complex nCenter, double complex oCenter, 
unsigned int deg, double complex *oCoeff, double complex *temp, double complex *z0){
    char * endptr;
    z0[0] = 1;
    double complex shift = oCenter - nCenter;
    for (int i = 1; i <= deg; ++i){
        z0[i] = z0[i-1]*shift;
    }

    temp[0] = oCoeff[0];
    for (int l = 1; l <= deg; l++){
        temp[l] = 0;
        for (int k=1; k <= l; k++){
            temp[l] += oCoeff[k]*z0[l-k]* strtod(CHOOSE_STRING[l-1][k-1], &endptr);
        } 
    }
    for (int i = 0; i <= deg; ++i){
        coeff[i] += temp[i];
    }
    return;
}

void shiftMultipoleExpansionLogDP(double complex* coeff, double complex nCenter, double complex oCenter, 
unsigned int deg, double complex *oCoeff, double complex *temp, double complex *z0){
    char *endptr;
    z0[0] = 1;
    for (int i = 1; i <= deg; i++){
        z0[i] = z0[i-1] * (oCenter - nCenter);
    }
    temp[0] = oCoeff[0];
    for (int l = 1; l <= deg; l++){
        temp[l] = -1.0 * temp[0] * z0[l] / l;
        for (int k = 1; k <= l; k++){
            temp[l] += oCoeff[k] * z0[l-k] * strtod(CHOOSE_STRING[l-1][k-1], &endptr);
        }
    }
    for (int i = 0; i <= deg; ++i){
        coeff[i] += temp[i];
    }
    return;
}

void shiftMultipoleExpansionDP(double complex* coeff, double complex nCenter, 
double complex oCenter, 
unsigned int deg, double complex *oCoeff, 
double complex *temp, double complex *z0, int m){
    if (m == 0){ // shift for expansion with sums of reciprocals
        shiftMultipoleExpansionReciprocalDP( coeff,  nCenter,  oCenter, deg, oCoeff, temp, z0);
    } else { // shift for expansion with sums of logs
        shiftMultipoleExpansionLogDP( coeff,  nCenter,  oCenter, deg, oCoeff, temp, z0);
    }
    return;
}

void multipoleToLocalReciprocalDP(double complex* coeff, double complex nCenter, double complex oCenter,
unsigned int md, unsigned int ld, double complex *oCoeff, 
double complex *temp, double complex *z){
/*
unsigned int md: multipole expansion degree
unsigned int ld: local expansion degree
*/
    char * endptr;
    unsigned int d;
    if (ld >= md){
        d = ld;
    } else {
        d = md;
    }
    z[0] = 1.0 + 0.0*I;
    double complex z0 = oCenter - nCenter;
    for (int i = 1; i <= d; i++){
        z[i] = z0 * z[i-1];
    }

    temp[0] = 0;
    double c = -1.0;
    for (int k = 1; k <= md; k++){
        temp[0] += c*oCoeff[k]/z[k];
        c *= -1.0;
    }

    for (int l = 1; l <= ld; l++){
        temp[l] = 0;
        c = -1.0;
        for (int k = 1; k <= md; k++){
            temp[l] += c * oCoeff[k] * ((strtod(CHOOSE_STRING[l+k-1][k-1], &endptr)/z[k])/z[l]);
            c *= -1.0;
        }
    }

    for (int i = 0; i <= ld; i++){
        coeff[i] += temp[i];
    }
    return;
}

void multipoleToLocalLogDP(double complex* coeff, double complex nCenter, double complex oCenter,
unsigned int md, unsigned int ld, double complex *oCoeff, 
double complex *temp, double complex *z){
    char * endptr;
    unsigned int d;
    if (ld >= md){
        d = ld;
    } else {
        d = md;
    }
    z[0] = 1.0 + 0.0*I;
    double complex z0 = oCenter - nCenter;
    for (int i = 1; i <= d; i++){
        z[i] = z0 * z[i-1];
    }

    temp[0] = oCoeff[0]*clog(-1.0*z0);
    double c = -1.0;
    for (int k = 1; k <= md; k++){
        temp[0] += c*oCoeff[k]/z[k];
        c *= -1.0;
    }

    for (int l = 1; l <= ld; l++){
        temp[l] = -1.0 * oCoeff[0] / (l *z[l]);
        c = -1.0;
        for (int k = 1; k <= md; k++){
            temp[l] += c * oCoeff[k] * (strtod(CHOOSE_STRING[l+k-1][k-1], &endptr)/z[k])/z[l];
            c *= -1;
        }
    }

    for (int i = 0; i <= ld; i++){
        coeff[i] += temp[i];
    }
    return;
}


void multipoleToLocalDP(double complex* coeff, double complex nCenter, 
double complex oCenter,
unsigned int md, unsigned int ld, double complex *oCoeff, 
double complex *temp, double complex *z, int m){
    if (m == 0){// reciprocal
        multipoleToLocalReciprocalDP(coeff, nCenter, oCenter,
        md, ld, oCoeff, 
        temp, z);        
    } else {// log
        multipoleToLocalLogDP(coeff, nCenter, oCenter,
        md, ld, oCoeff, 
        temp, z); 
    }
    return;
}


void shiftLocalExpansionDP(double complex *coeff, double complex nCenter, double complex oCenter,
unsigned int d, double complex *oCoeff, double complex *temp){
    double complex shift = oCenter - nCenter;
    for (int i = 0; i <= d; i++){
        temp[i] = oCoeff[i];
    }
    for (int j = 0; j <= d; j ++){
        for (int k = d - 1 - j; k < d; k++){
            temp[k] = temp[k] - shift*temp[k+1];
        }
    }

    for (int i = 0; i <= d; i++){
        coeff[i] += temp[i];
    }
    return;
}

double complex evalLocalExpansionDP(double complex z, double complex center,
unsigned int d, double complex *coeff){
    double complex res = coeff[d];
    double complex z0 = z - center;
    for (int k = d-1; k >= 0; k--){
        res = res*z0 + coeff[k];
    }
    return res;
}


double complex evalMultipoleExpansionLogDP(double complex z, double complex center,
unsigned int d, double complex *coeff){
    // double complex res = 0;
    // double complex shift = z - center;
    // double complex temp = z - center;
    // res += coeff[0]*clog(shift);
    // for (int k = 1; k <= d; k++){
    //     res += coeff[k]/temp;
    //     temp *= shift;
    // }
    // return res;
    double complex res = coeff[d];
    double complex shift = 1/(z - center);
    for(int k = d-1; k > 0; k--){
        res = res*shift + coeff[k];
    }
    res = res*shift - coeff[0]*clog(shift);
    return res;
}

double complex evalMultipoleExpansionReciprocalDP(double complex z, double complex center,
unsigned int d, double complex *coeff){
    // double complex res = 0.0 + 0.0*I;
    // double complex shift = z - center;
    // double complex temp = 1.0 + 0.0*I;
    // for (int k = 0; k <= d; k++){
    //     res += coeff[k]/temp;
    //     temp *= shift;
    // }

    double complex res = coeff[d];
    double complex shift = 1/(z-center);
    for (int k = d-1; k >= 0; k--){
        res = res*shift + coeff[k];
    }
    return res;    
}

double complex evalMultipoleExpansionDP(double complex z, double complex center,
unsigned int d, double complex *coeff, int m){
    if (m == 0) {
        return evalMultipoleExpansionReciprocalDP(z, center, d, coeff);
    } else {
        return evalMultipoleExpansionLogDP(z, center, d, coeff);
    }
}

/* 
Expansion computation based on multiple precision
*/


