#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "include/mp_fmm.h"



int fmmDP(char *sourceFileName, char* targetFileName, 
char* vectorFileName, char* resultFileName,
unsigned int med, unsigned int led,
unsigned long int mx, int mode
){
    // load sources
    double complex * sources = NULL;
    unsigned long int s = 0;
    unsigned long int ns;
    FILE *stream;
    stream = fopen(sourceFileName, "r");
    ns = complexDPReadall(stream, &sources, &s);
    if (ns == 0){
        printf("Not able to load source points as double complex.\n");
        fclose(stream);
        free(sources);
        return EXIT_FAILURE;
    }
    fclose(stream);

    // load targets
    double complex * targets = NULL;
    unsigned long int t = 0;
    unsigned long int nt;
    stream = fopen(targetFileName, "r");
    nt = complexDPReadall(stream, &targets, &t);
    if (nt == 0){
        printf("Not able to load target points as double complex.\n");
        fclose(stream);
        free(sources);
        free(targets);
        return EXIT_FAILURE;
    }
    fclose(stream);

    // load scalers
    double complex * vector = NULL;
    unsigned long int v = 0;
    unsigned long int nv;
    if (vectorFileName != NULL){
        stream = fopen(vectorFileName, "r");
        nv = complexDPReadall(stream, &vector, &v);
        if (nv == 0 || nv != ns){
            printf("Unable to load scalers as double complex or \n the number of source points and scalers do not match.\n");
            fclose(stream);
            free(sources);
            free(targets);
            free(vector);
            return EXIT_FAILURE;
        }
    } else {
        // set vector to default -- all scalers are 1.0 + 0.0*I
        vector = (double complex*) malloc(ns*sizeof(double complex));
        for (unsigned long int i = 0; i < ns; i++){
            vector[i] = CMPLX(1.0, 0.0);
        }
    }

    double complex *res;
    res = (double complex *)malloc(nt*sizeof(double complex));
    for (unsigned long int i = 0; i < nt; i++){
        res[i] = CMPLX(0, 0);
    }

    fast_multipole_method_DP(res, sources, targets, vector, ns, nt, med, led, mx, mode);
    stream = fopen(resultFileName, "w");
    unsigned long int r;
    r = complexDPWriteall(stream, res, nt);
    if (r != nt){
        printf("Error in wrting to result file.\n");
        fclose(stream);
        free(sources);
        free(targets);
        free(vector);
        free(res);
        return EXIT_FAILURE;
    }
    fclose(stream);
    free(sources);
    free(targets);
    free(vector);
    free(res);
    return EXIT_SUCCESS;
}


int fmmMP(char *sourceFileName, char* targetFileName, 
char* vectorFileName, char* resultFileName,
unsigned int med, unsigned int led,
unsigned long int mx, int mode, unsigned long int prec
){
    // load sources
    mpc_t * sources = NULL;
    unsigned long int s = 0;
    unsigned long int ns = 0;
    FILE *stream;
    stream = fopen(sourceFileName, "r");
    ns = complexMPReadall(stream, &sources, &s, prec);
    if (ns == 0){
        printf("Not able to load source points as multi-precision complex.\n");
        fclose(stream);
        free(sources);
        return EXIT_FAILURE;
    }
    fclose(stream);

    // load targets
    mpc_t * targets = NULL;
    unsigned long int t = 0;
    unsigned long int nt = 0;
    stream = fopen(targetFileName, "r");
    nt = complexMPReadall(stream, &targets, &t, prec);
    if (nt == 0){
        printf("Not able to load target points as multi-precision complex.\n");
        fclose(stream);
        free(sources);
        free(targets);
        return EXIT_FAILURE;
    }
    fclose(stream);

    // load scalers
    mpc_t * vector = NULL;
    unsigned long int v = 0;
    unsigned long int nv = 0;
    if (vectorFileName != NULL){
        stream = fopen(vectorFileName, "r");
        nv = complexMPReadall(stream, &vector, &v, prec);
        if (nv == 0 || nv != ns){
            printf("Unable to load scalers as multi-precision complex or \n the number of source points and scalers do not match.\n");
            fclose(stream);
            free(sources);
            free(targets);
            free(vector);
            return EXIT_FAILURE;
        }
    } else {
        // set vector to default -- all scalers are 1.0 + 0.0*I
        vector = (mpc_t *) malloc(ns*sizeof(mpc_t));
        for (unsigned long int i = 0; i < ns; i++){
            mpc_init2(vector[i], prec);
            mpc_set_d(vector[i], 1, MPC_RNDNN);
        }
    }

    mpc_t *res;
    res = (mpc_t *)malloc(nt*sizeof(mpc_t));
    for (unsigned long int i = 0; i < nt; i++){
        mpc_init2(res[i], prec);
    }

    fast_multipole_method_MP(res, sources, targets, vector, ns, nt, med, led, mx, mode, prec);
    stream = fopen(resultFileName, "w");
    unsigned long int r;
    r = complexMPWriteall(stream, res, nt);
    if (r != nt){
        printf("Error in wrting to result file.\n");
        fclose(stream);
        for (unsigned long int i = 0; i < ns; i++){
            mpc_clear(sources[i]);
            mpc_clear(vector[i]);
        }
        for (unsigned long int i = 0; i < nt; i++){
            mpc_clear(targets[i]);
            mpc_clear(res[i]);
        }
        free(sources);
        free(targets);
        free(vector);
        free(res);
        return EXIT_FAILURE;
    }
    for (unsigned long int i = 0; i < ns; i++){
        mpc_clear(sources[i]);
        mpc_clear(vector[i]);
    }
    for (unsigned long int i = 0; i < nt; i++){
        mpc_clear(targets[i]);
        mpc_clear(res[i]);
    }
    fclose(stream);
    free(sources);
    free(targets);
    free(vector);
    free(res);
    return EXIT_SUCCESS;
}


int main(int argc, char* argv[]){

/*
    -s filename        set filename for the input of sources
    -t filename        set filename for the input of targets
    -v filename        set filename for the input of scalers. 
                       All scaler are set to (1.0 + 0.0*I) if this file is not supplied.
    -o filename        set filename to store the output. 
    -m binary          set the kernel function with 
                       0 for K_x(y) = 1/(y-x) and 1 for K_x(y) = Log (y - x).
    -d degree          set the degree (uint) of multipole expansion. 
                       The default and minimum is 15 (and maximum is 25 for now).
    -l degree          set the degree (uint) of local expansion. 
                       The default and minimum is 20 (and maximum is 30 for now).
    -x number          set the maximum source number in the finest grids, 
                       i.e. a grid will not be further sub-divided if number of sources 
                       contained is less than this number. The default and minimum is 300.
    -p digits          set the number of digits used in computation. Double precision used 53 digits.   
*/

    char* sourceFileName, *vectorFileName, *targetFileName, *resultFileName;
    sourceFileName = vectorFileName = targetFileName = resultFileName = NULL;
    int mode = 0; // default mode is 0
    unsigned int med = 15; // default multipole expansion degree is 15
    unsigned int led = 20; // default local expansion degree is 20
    unsigned int tempUint;
    unsigned long int mx = 300; // default maximum source number in the finest grids is 300 
    unsigned long int tempUlint;
    unsigned long int prec = 53; // doube precision is used for prec <= 53 

    int i = 1;
    while (i < argc){
        if (argv[i][0] == '-') {
            char* arg = argv[i];
            int sind = 0; 
            int argind = 0;
            while (arg[++sind] != '\0') {
                switch(arg[sind]){
                    case 'm':
                        if ((i + (++argind)) < argc){
                            if (sscanf(argv[i + argind], "%d", &mode)!= 1){
                                printf("Invalid argument for -m option\n");
                                return EXIT_FAILURE;
                            }
                            if (mode != 0 && mode != 1){
                                printf("Invalid argument for -m option\n"
                                "0 for K_x(y) = 1/(y-x) and \n"
                                "1 for K_x(y) = log(y-x)\n");
                                return EXIT_FAILURE;
                            }
                        } else {
                            printf("No argument for -m option\n");
                            return EXIT_FAILURE;
                        }
                        break;                        
                    case 's': 
                        if ((i + (++argind)) < argc){
                            sourceFileName = argv[i + argind];
                        } else {
                            printf("No argument for -s option\n");
                            return EXIT_FAILURE;
                        }
                        break;
                    case 't':
                        if ((i + (++argind)) < argc){
                            targetFileName = argv[i + argind];
                        } else {
                            printf("No argument for -t option\n");
                            return EXIT_FAILURE;
                        }
                        break;
                    case 'v':
                        if ((i + (++argind)) < argc){
                            vectorFileName = argv[i + argind];
                        } else {
                            printf("No argument for -v option\n");
                            return EXIT_FAILURE;
                        }
                        break;
                    case 'o':
                        if ((i + (++argind)) < argc){
                            resultFileName = argv[i + argind];
                        } else {
                            printf("No argument for -o option\n");
                            return EXIT_FAILURE;
                        }
                        break;
                    case 'd':
                        if ((i + (++argind)) < argc){
                            if (sscanf(argv[i + argind], "%u", &tempUint) != 1){
                                printf("Invalid argument of -d option\n");
                                return EXIT_FAILURE;
                            }
                            med = (tempUint >= 15)? tempUint: 15;
                            med = (med <= 25)? med: 25;
                        } else {
                            printf("No argument for -d option\n");
                            return EXIT_FAILURE;
                        }
                        break;
                    case 'l':
                        if ((i + (++argind)) < argc){
                            if (sscanf(argv[i + argind], "%u", &tempUint) != 1){
                                printf("Invalid argument of -l option\n");
                                return EXIT_FAILURE;
                            }
                            led = (tempUint >= 20)? tempUint: 20;
                            led = (led <= 30)? led: 30;
                        } else {
                            printf("No argument for -l option\n");
                            return EXIT_FAILURE;
                        }
                        break;
                    case 'x':
                        if ((i + (++argind)) < argc){
                            if (sscanf(argv[i + argind], "%lu", &tempUlint) != 1){
                                printf("Invalid argument of -x option\n");
                                return EXIT_FAILURE;
                            }
                            mx = (tempUlint >= 300)? tempUlint: 300;
                        } else {
                            printf("No argument for -x option\n");
                            return EXIT_FAILURE;
                        }
                        break;
                    case 'p':
                        if ((i + (++argind)) < argc){
                            if (sscanf(argv[i + argind], "%lu", &tempUlint) != 1){
                                printf("Invalid argument of -p option\n");
                                return EXIT_FAILURE;
                            }
                            prec = (tempUlint >= 53)? tempUlint: 53;
                        } else {
                            printf("No argument for -p option\n");
                            return EXIT_FAILURE;
                        }
                        break;
                    default:
                        printf("Invalid option. \n");
                        return EXIT_FAILURE;
                }
            }
            i += ++argind;
        } else {
            i++;
        }
    }
    if (sourceFileName == NULL) {
        printf("No input file specified for Sources.\n");
        return EXIT_FAILURE;
    }

    if (targetFileName == NULL){
        printf("No input file specified for Targets.\n");
        return EXIT_FAILURE;
    }


    if (resultFileName == NULL){
        resultFileName = "tempRes.txt";
        printf("Output filename not specified, and is set to %s\n", resultFileName);
    }
    // printf("Output filename: %s\n", resultFileName);
    // printf("Degree of multipole expansion is %u\n", med);
    // printf("Degree of local expansion is %u\n", led);
    // printf("Max number of sources in the finest grid is %lu\n", mx);


    if (prec <= 53) { // use double precision for computation
        fmmDP(sourceFileName, targetFileName, vectorFileName, resultFileName, med, led, mx, mode);
    } else {
        fmmMP(sourceFileName, targetFileName, vectorFileName, resultFileName, med, led, mx, mode, prec);
    }

    printf("Computation results are written to file %s\n", resultFileName);

    return EXIT_SUCCESS;
}

