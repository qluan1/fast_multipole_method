#ifndef mpfmm_h
#define mpfmm_h


#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#include <time.h>
#include <string.h>


/* 

from aux.c

*/

struct TreeNode {

    unsigned long int numSources;
    unsigned long int * sourcesPool;
    unsigned long int numTargets;
    unsigned long int * targetsPool;

    int isDoublePrecision;
    double radiusDouble;
    double centerRealDouble, centerImagDouble;

    unsigned long int precision;
    mpfr_t radiusMP;
    mpfr_t centerRealMP, centerImagMP;

    unsigned int multipoleExpansionDegree;
    double complex *multipoleExpansionDouble;
    mpc_t *multipoleExpansionMP;

    unsigned int localExpansionDegree;
    double complex *localExpansionDouble;
    mpc_t *localExpansionMP;

    struct TreeNode * ch1;
    struct TreeNode * ch2;
    struct TreeNode * ch3;
    struct TreeNode * ch4;
    struct TreeNode * parent;

    unsigned int numNearNeighbors;
    struct TreeNode** nearNeighbors;
    unsigned int numInteractions;
    struct TreeNode** interactions;
};

typedef struct TreeNode TreeNode_t[1];

struct TreeNode* newTreeNodeDP(struct TreeNode* par,
    double radiusDouble, double centerRealDouble, double centerImagDouble);

void freeTreeNode(struct TreeNode* node);

void printTreeNodeDP(struct TreeNode *node);

void getTreeNodeChildren(struct TreeNode *node, struct TreeNode **children);

void freeTree(struct TreeNode* root);

int isNearNeighborDP(struct TreeNode *n1, struct TreeNode *n2);

void nodeNeighborInteractionDP(struct TreeNode *node);

void generateNeighbourInteractionDP(struct TreeNode * root);

void upwardPassDP(struct TreeNode *root, 
double complex * sources, 
double complex * vector,
unsigned int multi_deg,
int maxNumSources,
int m);

void downwardPassDP(
struct TreeNode *root,
unsigned int local_deg,
int m // m = 0 -> sum of reciprocals
);

struct TreeNode * createTreeDP(
double complex *sourcesD,
double complex *vectorD,
double complex *targetsD,
unsigned long int numSrc,
unsigned long int numTgt,
unsigned long int maxs,
struct TreeNode ** nodeMap
);

void subdivideDP(struct TreeNode *node, struct TreeNode *ch1, 
struct TreeNode *ch2, struct TreeNode *ch3, struct TreeNode *ch4, 
double complex * sources, double complex *targets);

void evaluateDP(double complex *res, 
double complex *sources, unsigned long int ns,
double complex *vector,
double complex *targets, unsigned long int nt,
struct TreeNode *root, int m);

void evaluateLEDP(double complex *res, 
double complex *sources, unsigned long int ns,
double complex *vector, struct TreeNode **map,
double complex *targets, unsigned long int nt,
struct TreeNode *root, int m);

void fast_multipole_method_DP(double complex *res,
double complex *sources, double complex *targets,
double complex *vector, 
unsigned long int ns, unsigned long int nt, 
unsigned int med, unsigned int led, 
unsigned long int mx, int mode);

struct StackNode;

struct Stack;

struct StackNode* newStackNode(void * node);

struct Stack* newStack();

int stackIsEmpty(struct Stack* s);

void stackPush(struct Stack* s, void* t);

void stackPushLeft(struct Stack *s, void* t);

void* stackPop(struct Stack* s);

void* stackPopLeft(struct Stack* s);

void stackFreeStack(struct Stack* s);

unsigned long int complexDPReadall(FILE *input, 
double complex **dataptr, unsigned long int *sizeptr);

unsigned long int complexDPWriteall(FILE* output, 
double complex * data, unsigned long int n);



/*
    from treeMP.c
*/



unsigned long int complexMPReadall(FILE *input, 
mpc_t **dataptr, unsigned long int *sizeptr, unsigned long int prec);

unsigned long int complexMPWriteall(FILE* output, 
mpc_t *data, unsigned long int n);

struct TreeNode* newTreeNodeMP(struct TreeNode *par,
mpfr_t radiusMP, mpfr_t centerRealMP, mpfr_t centerImagMP,
unsigned long int prec);

struct TreeNode * createTreeMP(
mpc_t *sourcesM,
mpc_t *vectorM,
mpc_t *targetsM,
unsigned long int numSrc,
unsigned long int numTgt,
unsigned long int maxs,
struct TreeNode ** nodeMap,
unsigned long int prec
);

void subdivideMP(struct TreeNode *node, struct TreeNode *ch1, 
struct TreeNode *ch2, struct TreeNode *ch3, struct TreeNode *ch4, 
mpc_t * sources, mpc_t *targets);

int isNearNeighborMP(struct TreeNode *n1, struct TreeNode *n2);

void nodeNeighborInteractionMP(struct TreeNode *node);

void generateNeighbourInteractionMP(struct TreeNode * root);

void upwardPassMP(struct TreeNode *root,
    mpc_t *sources,
    mpc_t *vector,
    unsigned int multi_deg,
    int maxNumSources,
    unsigned long int prec,
    int mode);

void downwardPassMP(struct  TreeNode *root,
unsigned long int local_deg, unsigned long int prec, int mode);

void evaluateLEMP(mpc_t *res, 
mpc_t *sources, unsigned long int ns,
mpc_t *vector, struct TreeNode **map,
mpc_t *targets, unsigned long int nt,
struct TreeNode *root, unsigned long int prec, int mode);


void fast_multipole_method_MP(mpc_t *res,
mpc_t *sources, mpc_t *targets,
mpc_t *vector, 
unsigned long int ns, unsigned long int nt, 
unsigned int med, unsigned int led, 
unsigned long int mx, int mode, unsigned long int prec);
/* 

from expansion.c

*/



void multipoleExpansionReciprocalDP(
double complex *coeff, double complex center, unsigned int deg, 
unsigned long int numSources, unsigned long int *sourceIndices, 
double complex * sources, double complex *vec, double complex *tempS);

void multipoleExpansionLogDP(
double complex*coeff, double complex center, unsigned int deg, 
unsigned long int numSources, unsigned long int *sourcesIndices,
double complex * sources, double complex * vec);

void multipoleExpansionDP(
double complex*coeff, double complex center, unsigned int deg, 
unsigned long int numSources, unsigned long int *sourcesIndices,
double complex * sources, double complex * vec, int m, double complex *tempS);


void shiftMultipoleExpansionReciprocalDP(
double complex* coeff, double complex nCenter, double complex oCenter, 
unsigned int deg, double complex *oCoeff, double complex *temp, 
double complex *z0);

void shiftMultipoleExpansionLogDP(
double complex* coeff, double complex nCenter, double complex oCenter, 
unsigned int deg, double complex *oCoeff, double complex *temp, 
double complex *z0);

void shiftMultipoleExpansionDP(
double complex* coeff, double complex nCenter, double complex oCenter, 
unsigned int deg, double complex *oCoeff, double complex *temp, 
double complex *z0, int m);

void multipoleToLocalReciprocalDP(double complex* coeff, double complex nCenter, double complex oCenter,
unsigned int md, unsigned int ld, double complex *oCoeff, 
double complex *temp, double complex *z);

void multipoleToLocalLogDP(double complex* coeff, double complex nCenter, double complex oCenter,
unsigned int md, unsigned int ld, double complex *oCoeff, 
double complex *temp, double complex *z);

void multipoleToLocalDP(double complex* coeff, double complex nCenter, double complex oCenter,
unsigned int md, unsigned int ld, double complex *oCoeff, 
double complex *temp, double complex *z, int m);

void shiftLocalExpansionDP(double complex *coeff, double complex nCenter, double complex oCenter,
unsigned int d, double complex *oCoeff, double complex *temp);

double complex evalLocalExpansionDP(double complex z, double complex center,
unsigned int d, double complex *coeff);

double complex evalMultipoleExpansionLogDP(double complex z, double complex center,
unsigned int d, double complex *coeff);

double complex evalMultipoleExpansionReciprocalDP(double complex z, double complex center,
unsigned int d, double complex *coeff);

double complex evalMultipoleExpansionDP(double complex z, double complex center,
unsigned int d, double complex *coeff, int m);


/* 

from expansionMP.c

*/

void multipoleExpansionReciprocalMP(mpc_t *coeff, mpc_t center,
unsigned long int deg,
unsigned long int numSources, unsigned long int *sourcesIndices,
mpc_t *sources, mpc_t *vec, 
mpc_t *tempVec, unsigned long int tempLen);


void multipoleExpansionLogMP(mpc_t *coeff, mpc_t center, 
unsigned int deg,  
unsigned long int numSources, unsigned long int *sourcesIndices,
mpc_t * sources, mpc_t * vec,
mpc_t * tempVec, unsigned long int tempLen);


void multipoleExpansionMP(mpc_t *coeff, mpc_t center, 
unsigned int deg,  
unsigned long int numSources, unsigned long int *sourcesIndices,
mpc_t * sources, mpc_t * vec,
mpc_t * tempVec, unsigned long int tempLen, int m);

void shiftMultipoleExpansionReciprocalMP(mpc_t *coeff, 
mpc_t nCenter, mpc_t oCenter, 
unsigned int deg, mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose);


void shiftMultipoleExpansionLogMP(mpc_t *coeff, 
mpc_t nCenter, mpc_t oCenter, 
unsigned int deg, mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose);


void shiftMultipoleExpansionMP(mpc_t *coeff, 
mpc_t nCenter, mpc_t oCenter, 
unsigned int deg, mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose, int mode);


void multipoleToLocalReciprocalMP(mpc_t *coeff,
mpc_t nCenter, mpc_t oCenter,
unsigned long int md, unsigned long int ld,
mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose);

void multipoleToLocalLogMP(mpc_t *coeff,
mpc_t nCenter, mpc_t oCenter,
unsigned long int md, unsigned long int ld,
mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose);

void multipoleToLocalMP(mpc_t *coeff,
mpc_t nCenter, mpc_t oCenter,
unsigned long int md, unsigned long int ld,
mpc_t *oCoeff, 
mpc_t *tempVec, unsigned long int tempLen,
mpc_t *zVec, unsigned long int zLen,
mpz_t choose, int mode);

void shiftLocalExpansionMP(mpc_t *coeff,
mpc_t nCenter, mpc_t oCenter,
unsigned int d, mpc_t *oCoeff,
mpc_t *tempVec, unsigned long int tempLen,
mpc_t shift);

void evalLocalExpansionMP(mpc_t res, mpc_t z,
mpc_t center, unsigned int d, mpc_t *coeff,
mpc_t shift);

void evalMultipoleExpansionLogMP(mpc_t res, mpc_t z,
mpc_t center, unsigned int d, mpc_t *coeff);

void evalMultipoleExpansionReciprocalMP(mpc_t res, mpc_t z,
mpc_t center, unsigned int d, mpc_t *coeff);

void evalMultipoleExpansionMP(mpc_t res, mpc_t z,
mpc_t center, unsigned int d, 
mpc_t *coeff, int mode);



/* 

from const.c

*/
extern const char* CHOOSE_STRING[61][61];
#endif