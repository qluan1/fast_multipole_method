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


/* 
    Reading and writing multiple precision complex numbers from files.
*/

unsigned long int complexMPReadall(FILE *input, mpc_t **dataptr, unsigned long int *sizeptr, unsigned long int prec)
{
    mpc_t *data;
    unsigned long int size;
    unsigned long int used = 0;
    size_t read;

    if (!input || !dataptr || !sizeptr) {
        /* At least one of the parameters is NULL. */
        return 0;
    }

    if (ferror(input)) {
        /* input stream is already in error state. */
        return 0;
    }

    if (!*dataptr || !*sizeptr) {
        /* *dataptr is NULL, or *sizeptr == 0,
           so we initialize them to empty. */
        *dataptr = NULL;
        *sizeptr = 0;
    }
    data = *dataptr;
    size = *sizeptr;    

    while (1) {
        if (used >= size) {
            /* We need to grow the data array. */

            /* Simple allocation policy:
               allocate in sets of roughly 1024 vectors. */
            size = (used | 1023) + 1021;
            data = realloc(data, size * sizeof *data);
            if (!data) {
                /* Realloc failed! */
                printf("Realloc failed\n");
                return 0;
            }

            *dataptr = data;
            *sizeptr = size;

            for (unsigned long int i = used; i < size; i++){
                mpc_init2(data[i], prec);
            }
        }
        if (mpc_inp_str(data[used], input, &read, 10, MPC_RNDNN) == -1){
            break;
        }
        /* One more MP complex number read and saved */
        used++;
    }

    /* If there was an actual I/O error, or
       the file contains unread data, set used to 0. */
    if ( ferror(input) || !feof(input)){
        printf("I/O error\n");
        used = 0;
    }
    return used;    
}

unsigned long int complexMPWriteall(FILE* output, mpc_t *data, unsigned long int n){
    unsigned long int used = 0;
    for (unsigned long int i = 0; i < n; i++){
        if (mpc_out_str(output, 10, 0, data[i], MPC_RNDNN) <= 0){
            printf("error writing a complex number\n");
            return 0;
        }
        fprintf(output, "\n");
        used++;
    }
    return used;
}

/*
    Tree Structure functions for multiple precision computation.
*/

struct TreeNode* newTreeNodeMP(struct TreeNode *par,
mpfr_t radiusMP, mpfr_t centerRealMP, mpfr_t centerImagMP,
unsigned long int prec){
    struct TreeNode * temp = (struct TreeNode*)malloc(sizeof(struct TreeNode));

    temp->isDoublePrecision = 0;
    temp->precision = prec;
    mpfr_init2(temp->radiusMP, prec);
    mpfr_init2(temp->centerRealMP, prec);
    mpfr_init2(temp->centerImagMP, prec);
    mpfr_set(temp->radiusMP, radiusMP, MPFR_RNDN);
    mpfr_set(temp->centerRealMP, centerRealMP, MPFR_RNDN);
    mpfr_set(temp->centerImagMP, centerImagMP, MPFR_RNDN);

    temp->parent = par;
    temp->ch1 = NULL;
    temp->ch2 = NULL;
    temp->ch3 = NULL;
    temp->ch4 = NULL;
    temp->numSources = 0;
    temp->numTargets = 0;
    temp->numNearNeighbors = 0;
    temp->nearNeighbors = (struct TreeNode **)malloc(9*sizeof(struct TreeNode*));
    temp->numInteractions = 0;
    temp->interactions = (struct TreeNode **)malloc(27*sizeof(struct TreeNode*));
    temp->multipoleExpansionDegree = 0;
    temp->localExpansionDegree = 0;
    temp->multipoleExpansionMP = NULL;
    temp->localExpansionMP = NULL;
    return temp;
}

struct TreeNode * createTreeMP(
mpc_t *sourcesM,
mpc_t *vectorM,
mpc_t *targetsM,
unsigned long int numSrc,
unsigned long int numTgt,
unsigned long int maxs,
struct TreeNode ** nodeMap,
unsigned long int prec
){
    mpfr_t min_x, max_x, min_y, max_y, radius, x, y;
    mpfr_init2(min_x, prec);
    mpfr_init2(max_x, prec);
    mpfr_init2(min_y, prec);
    mpfr_init2(max_y, prec);
    mpfr_init2(radius, prec);
    mpfr_init2(x, prec);
    mpfr_init2(y, prec);

    mpfr_set(min_x, mpc_realref(sourcesM[0]), MPFR_RNDN);
    mpfr_set(max_x, mpc_realref(sourcesM[0]), MPFR_RNDN);
    mpfr_set(min_y, mpc_imagref(sourcesM[0]), MPFR_RNDN);
    mpfr_set(max_y, mpc_imagref(sourcesM[0]), MPFR_RNDN);


    for (unsigned long int i = 0; i < numSrc; ++i){
        if (mpfr_cmp(min_x, mpc_realref(sourcesM[i])) > 0){ // min_x > Re(sources[i])
            mpfr_set(min_x, mpc_realref(sourcesM[i]), MPFR_RNDN);
        } else if (mpfr_cmp(mpc_realref(sourcesM[i]), max_x) > 0 ){ // max_x < Re(sources[i])
            mpfr_set(max_x, mpc_realref(sourcesM[i]), MPFR_RNDN);
        }
        if (mpfr_cmp(min_y, mpc_imagref(sourcesM[i])) > 0){ // min_x > Re(sources[i])
            mpfr_set(min_y, mpc_imagref(sourcesM[i]), MPFR_RNDN);
        } else if (mpfr_cmp(mpc_imagref(sourcesM[i]), max_y) > 0 ){ // max_x < Re(sources[i])
            mpfr_set(max_y, mpc_imagref(sourcesM[i]), MPFR_RNDN);
        }        
    }

    for (unsigned long int i = 0; i < numTgt; ++i){
        if (mpfr_cmp(min_x, mpc_realref(targetsM[i])) > 0){ // min_x > Re(sources[i])
            mpfr_set(min_x, mpc_realref(targetsM[i]), MPFR_RNDN);
        } else if (mpfr_cmp(mpc_realref(targetsM[i]), max_x) > 0 ){ // max_x < Re(sources[i])
            mpfr_set(max_x, mpc_realref(targetsM[i]), MPFR_RNDN);
        }
        if (mpfr_cmp(min_y, mpc_imagref(targetsM[i])) > 0){ // min_x > Re(sources[i])
            mpfr_set(min_y, mpc_imagref(targetsM[i]), MPFR_RNDN);
        } else if (mpfr_cmp(mpc_imagref(targetsM[i]), max_y) > 0 ){ // max_x < Re(sources[i])
            mpfr_set(max_y, mpc_imagref(targetsM[i]), MPFR_RNDN);
        }        
    }

    // compute radius of the bounding box
    // and its center
    mpfr_add(x, min_x, max_x, MPFR_RNDN);
    mpfr_mul_d(x, x, 0.5, MPFR_RNDN);
    mpfr_add(y, min_y, max_y, MPFR_RNDN);
    mpfr_mul_d(y, y, 0.5, MPFR_RNDN);
    
    mpfr_sub(radius, max_x, min_y, MPFR_RNDN);
    mpfr_sub(max_y, max_y, min_y, MPFR_RNDN);
    if ( mpfr_cmp(max_y, radius) > 0){
        mpfr_set(radius, max_y, MPFR_RNDN);
    }
    mpfr_mul_d(radius, radius, 0.5, MPFR_RNDN);


    struct Stack *cur, *nxt, *tmp;
    struct TreeNode *node, *ch1, * ch2, *ch3, *ch4;


    struct TreeNode *root;
    root = newTreeNodeMP(NULL, radius, x, y, prec);
    root->numSources = numSrc;
    root->sourcesPool = (unsigned long int *)malloc(numSrc*sizeof(unsigned long int));
    for (unsigned long int i = 0; i < numSrc; i++){
        root->sourcesPool[i] = i;
    }

    root->numTargets = numTgt;
    root->targetsPool = (unsigned long int *)malloc(numTgt*sizeof(unsigned long int));
    for (unsigned long int i = 0; i < numTgt; i++){
        root->targetsPool[i] =i;
    }

    cur = newStack();
    stackPush(cur, root);
    nxt = newStack();
    while(!stackIsEmpty(cur)){
        while (!stackIsEmpty(cur)){
            node = (struct TreeNode*)stackPop(cur);
            if (node->numSources > maxs) { //continue subdividing
                mpfr_mul_d(radius, node->radiusMP, 0.5, MPFR_RNDN);
                mpfr_add(x, node->centerRealMP, radius, MPFR_RNDN);
                mpfr_add(y, node->centerImagMP, radius, MPFR_RNDN);
                ch1 = newTreeNodeMP( node, radius, x, y, prec);
                node->ch1 = ch1;
                mpfr_add(x, node->centerRealMP, radius, MPFR_RNDN);
                mpfr_sub(y, node->centerImagMP, radius, MPFR_RNDN);
                ch2 = newTreeNodeMP( node, radius, x, y, prec);
                node->ch2 = ch2;
                mpfr_sub(x, node->centerRealMP, radius, MPFR_RNDN);
                mpfr_add(y, node->centerImagMP, radius, MPFR_RNDN);
                ch3 = newTreeNodeMP( node, radius, x, y, prec);
                node->ch3 = ch3;
                mpfr_sub(x, node->centerRealMP, radius, MPFR_RNDN);
                mpfr_sub(y, node->centerImagMP, radius, MPFR_RNDN);
                ch4 = newTreeNodeMP( node, radius, x, y, prec);
                node->ch4 = ch4;
                stackPush(nxt, ch1);
                stackPush(nxt, ch2);
                stackPush(nxt, ch3);
                stackPush(nxt, ch4);
                subdivideMP(node, ch1, ch2, ch3, ch4, sourcesM, targetsM);
            } else { // considered a leaf
                // map the targets to leaf node
                for (unsigned long int i  = 0; i < node->numTargets; ++i){
                    unsigned long int idx = node->targetsPool[i];
                    nodeMap[idx] = node;
                }
                if (node->numSources == 0 && node->numTargets == 0){
                    freeTreeNode(node);
                }
            }
        }
        tmp = cur;
        cur = nxt;
        nxt = tmp;
        tmp = NULL;
    }
    stackFreeStack(cur);
    stackFreeStack(nxt);
    stackFreeStack(tmp);
    mpfr_clear(min_x);
    mpfr_clear(max_x);
    mpfr_clear(min_y);
    mpfr_clear(max_y);
    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(radius);
    return root;
}


void subdivideMP(struct TreeNode *node, struct TreeNode *ch1, 
struct TreeNode *ch2, struct TreeNode *ch3, struct TreeNode *ch4, 
mpc_t * sources, mpc_t *targets){
    // wipe the source pool of node and distribute the sources to children nodes
    ch1->sourcesPool = (unsigned long int *)malloc(node->numSources*sizeof(unsigned long int));
    ch2->sourcesPool = (unsigned long int *)malloc(node->numSources*sizeof(unsigned long int));
    ch3->sourcesPool = (unsigned long int *)malloc(node->numSources*sizeof(unsigned long int));
    ch4->sourcesPool = (unsigned long int *)malloc(node->numSources*sizeof(unsigned long int));

    for (unsigned long int i = 0; i < node->numSources; i++){
        unsigned long int idx = node->sourcesPool[i];
        //mpc_realref(sources[idx]);
        //mpc_imagref(sources[idx]);
        if( mpfr_cmp(mpc_realref(sources[idx]), node->centerRealMP) > 0) {
            if ( mpfr_cmp(mpc_imagref(sources[idx]), node->centerImagMP) > 0) {
                // put source in ch1
                ch1->sourcesPool[ch1->numSources] = idx;
                ch1->numSources += 1;
            } else {
                // put source in ch2
                ch2->sourcesPool[ch2->numSources] = idx;
                ch2->numSources += 1;
            } 
        } else {
            if (mpfr_cmp(mpc_imagref(sources[idx]), node->centerImagMP) > 0) {
                // put source in ch3
                ch3->sourcesPool[ch3->numSources] = idx;
                ch3->numSources += 1;
            } else {
                // put source in ch4
                ch4->sourcesPool[ch4->numSources] = idx;
                ch4->numSources += 1;
            }
        }
    }
    free(node->sourcesPool);
    node->sourcesPool = NULL;
    node->numSources = 0;

    ch1->sourcesPool = (unsigned long int*)realloc(ch1->sourcesPool, 
    ch1->numSources * sizeof(unsigned long int));
    ch2->sourcesPool = (unsigned long int*)realloc(ch2->sourcesPool, 
    ch2->numSources * sizeof(unsigned long int));
    ch3->sourcesPool = (unsigned long int*)realloc(ch3->sourcesPool, 
    ch3->numSources * sizeof(unsigned long int));
    ch4->sourcesPool = (unsigned long int*)realloc(ch4->sourcesPool, 
    ch4->numSources * sizeof(unsigned long int));

    ch1->targetsPool = (unsigned long int *)malloc(node->numTargets*sizeof(unsigned long int));
    ch2->targetsPool = (unsigned long int *)malloc(node->numTargets*sizeof(unsigned long int));
    ch3->targetsPool = (unsigned long int *)malloc(node->numTargets*sizeof(unsigned long int));
    ch4->targetsPool = (unsigned long int *)malloc(node->numTargets*sizeof(unsigned long int));

    for (unsigned long int i = 0; i < node->numTargets; i ++){
        unsigned long int idx = node->targetsPool[i];
        if( mpfr_cmp(mpc_realref(targets[idx]), node->centerRealMP) > 0) {
            if (mpfr_cmp(mpc_imagref(targets[idx]), node->centerImagMP) > 0) {
                // put target in ch1
                ch1->targetsPool[ch1->numTargets] = idx;
                ch1->numTargets += 1;
            } else {
                // put target in ch2
                ch2->targetsPool[ch2->numTargets] = idx;
                ch2->numTargets += 1;
            } 
        } else {
            if (mpfr_cmp(mpc_imagref(targets[idx]), node->centerImagMP) > 0) {
                // put target in ch3
                ch3->targetsPool[ch3->numTargets] = idx;
                ch3->numTargets += 1;
            } else {
                // put target in ch4
                ch4->targetsPool[ch4->numTargets] = idx;
                ch4->numTargets += 1;
            }
        }
    }
    free(node->targetsPool);
    node->numTargets = 0;
    node->targetsPool = NULL;
    
    ch1->targetsPool = (unsigned long int*)realloc(ch1->targetsPool, 
    ch1->numTargets * sizeof(unsigned long int));
    ch2->targetsPool = (unsigned long int*)realloc(ch2->targetsPool, 
    ch2->numTargets * sizeof(unsigned long int));
    ch3->targetsPool = (unsigned long int*)realloc(ch3->targetsPool, 
    ch3->numTargets * sizeof(unsigned long int));
    ch4->targetsPool = (unsigned long int*)realloc(ch4->targetsPool, 
    ch4->numTargets * sizeof(unsigned long int));    

    return;
}




int isNearNeighborMP(struct TreeNode *n1, struct TreeNode *n2){
    /* 
       check if two tree nodes are near neighbors:
       (1) same node
       (2) share an endpoint
       (3) share a side

       return 1 if they are NN,
       return 0 otherwise.
       The tree nodes must be on the same level.
       ***This function should not be used alone.***
    */
   if (n1 == n2){
       return 1;
   }

    int v, h;
    mpfr_mul_d(n1->radiusMP, n1->radiusMP, 2.5, MPFR_RNDN); // replace n1->radiusMP with threshold 
    if (mpfr_cmp(n1->centerRealMP, n2->centerRealMP) > 0 ){ // 
        mpfr_sub(n1->centerRealMP, n1->centerRealMP, n2->centerRealMP, MPFR_RNDN);
        v = (mpfr_cmp(n1->centerRealMP, n1->radiusMP) < 0)? 1 : 0;
        mpfr_add(n1->centerRealMP, n1->centerRealMP, n2->centerRealMP, MPFR_RNDN);
    } else {
        mpfr_sub(n2->centerRealMP, n2->centerRealMP, n1->centerRealMP, MPFR_RNDN);
        v = (mpfr_cmp(n2->centerRealMP, n1->radiusMP) < 0)? 1 : 0;
        mpfr_add(n2->centerRealMP, n2->centerRealMP, n1->centerRealMP, MPFR_RNDN);
    }

    if (mpfr_cmp(n1->centerImagMP, n2->centerImagMP) > 0 ){ // 
        mpfr_sub(n1->centerImagMP, n1->centerImagMP, n2->centerImagMP, MPFR_RNDN);
        h = (mpfr_cmp(n1->centerImagMP, n1->radiusMP) < 0)? 1 : 0;
        mpfr_add(n1->centerImagMP, n1->centerImagMP, n2->centerImagMP, MPFR_RNDN);
    } else {
        mpfr_sub(n2->centerImagMP, n2->centerImagMP, n1->centerImagMP, MPFR_RNDN);
        h = (mpfr_cmp(n2->centerImagMP, n1->radiusMP) < 0)? 1 : 0;
        mpfr_add(n2->centerImagMP, n2->centerImagMP, n1->centerImagMP, MPFR_RNDN);
    }
    mpfr_mul_d(n1->radiusMP, n1->radiusMP, 0.4, MPFR_RNDN);
    return (v * h);
}

void nodeNeighborInteractionMP(struct TreeNode *node){
    /*
    find nearby nodes on the same level and classify them
    into either Near Neighbor or Interaction
    */
    int numCandidates = 0;
    struct TreeNode ** candidates = (struct TreeNode**)malloc(36*sizeof(struct TreeNode *));
    struct TreeNode * par, * tmp; 
    struct TreeNode **children = (struct TreeNode **)malloc(4*sizeof(struct TreeNode *));
    par = node->parent;
    if (par != NULL){
        for (int i = 0; i < par->numNearNeighbors; ++i){
            tmp = par->nearNeighbors[i];
            getTreeNodeChildren(tmp, children);
            for (int c = 0; c < 4; c++){
                tmp = children[c];
                if (tmp != NULL){
                    candidates[numCandidates] = tmp;
                    numCandidates++;
                }             
            }     
        }

        for (int i = 0; i < numCandidates; i++){
            tmp = candidates[i];
            if (isNearNeighborMP(tmp, node) == 1){
                node->nearNeighbors[node->numNearNeighbors] = tmp;
                node->numNearNeighbors++;
            } else {
                node->interactions[node->numInteractions] = tmp;
                node->numInteractions++;
            }
        }
    } else { // node is a root
        node->nearNeighbors[node->numNearNeighbors] = node;
        node->numNearNeighbors++;
    }
    free(candidates);
    free(children);
}

void generateNeighbourInteractionMP(struct TreeNode * root){
    struct Stack * stack = newStack();
    struct TreeNode* node;
    stackPush(stack, root);
    while (!stackIsEmpty(stack)){
        node = (struct TreeNode*)stackPop(stack);
        nodeNeighborInteractionMP(node);
        if (node->ch1 != NULL) {
            stackPush(stack, node->ch1);
        }
        if (node->ch2 != NULL) {
            stackPush(stack, node->ch2);
        }
        if (node->ch3 != NULL) {
            stackPush(stack, node->ch3);
        }
        if (node->ch4 != NULL) {
            stackPush(stack, node->ch4);
        }
    }
    stackFreeStack(stack);
    return;
}

void upwardPassMP(struct TreeNode *root,
    mpc_t *sources,
    mpc_t *vector,
    unsigned int multi_deg,
    int maxNumSources,
    unsigned long int prec,
    int mode){
    /* 
        mode is used as a selector with 
        0: fmm for sums of reciprocals
        1: fmm for sums of complex logs
    */

   unsigned long int tempLen = 0;
   unsigned long int zLen = 0;
   tempLen = (tempLen < maxNumSources + 1) ? maxNumSources + 1 : tempLen;
   tempLen = (tempLen < multi_deg + 2) ? multi_deg + 2 : tempLen;
   zLen = (zLen < multi_deg + 1) ? multi_deg + 1 : zLen;

    mpc_t center, oCenter;
    mpc_t *tempVec = (mpc_t *)malloc(tempLen*sizeof(mpc_t));
    mpc_t *zVec = (mpc_t *)malloc(zLen*sizeof(mpc_t));
    mpz_t choose;

    mpc_init2(center, prec);
    mpc_init2(oCenter, prec);
    for (unsigned long int i = 0; i < tempLen; i++){
        mpc_init2(tempVec[i], prec);
    }
    for (unsigned long int i = 0; i < zLen; i++){
        mpc_init2(zVec[i], prec);
    }
    mpz_init2(choose, prec);

    struct Stack *s, *rev;
    s = newStack();
    rev = newStack();
    /* 
    create a reverse list of the tree nodes from bottom (with no leaves) to root
    compute the initial expansion at non-empty leaves
    */

    stackPush(s, root);
    struct TreeNode *node, *child;

    while (!stackIsEmpty(s)){
        node = (struct TreeNode*)stackPop(s);
        if (node->ch1 != NULL || node->ch2 != NULL || 
        node->ch3 != NULL || node->ch4 != NULL){ // not a leaf node
            stackPush(rev, node);
            if (node->ch1) stackPushLeft(s, node->ch1);
            if (node->ch2) stackPushLeft(s, node->ch2);
            if (node->ch3) stackPushLeft(s, node->ch3);
            if (node->ch4) stackPushLeft(s, node->ch4);
        } else { // is a leaf node
            mpc_set_fr_fr(center, node->centerRealMP, node->centerImagMP, MPC_RNDNN);
            node->multipoleExpansionDegree = multi_deg;
            node->multipoleExpansionMP = 
            (mpc_t *)malloc((multi_deg + 1)*sizeof(mpc_t));
            for (int i = 0; i <= multi_deg; i++){
                mpc_init2(node->multipoleExpansionMP[i], prec);
                mpc_set_d(node->multipoleExpansionMP[i], 0, MPC_RNDNN);
            }
            multipoleExpansionMP(node->multipoleExpansionMP, center, 
            multi_deg, node->numSources,
            node->sourcesPool, sources, vector, tempVec, tempLen, mode);
        }
    }

    struct TreeNode ** children = (struct TreeNode**)malloc(4*sizeof(struct TreeNode*));

    while(!stackIsEmpty(rev)){
        node = (struct TreeNode*)stackPop(rev);
        mpc_set_fr_fr(center, node->centerRealMP, node->centerImagMP, MPC_RNDNN);
        node->multipoleExpansionDegree = multi_deg;
        node->multipoleExpansionMP = 
        (mpc_t *)malloc((multi_deg + 1)*sizeof(mpc_t));
        for (int i = 0; i <= multi_deg; i++){
            mpc_init2(node->multipoleExpansionMP[i], prec);
            mpc_set_d(node->multipoleExpansionMP[i], 0, MPC_RNDNN);
        }
        getTreeNodeChildren(node, children);
        for (int i = 0; i < 4; i++){
            child = children[i];
            if (child != NULL){
                mpc_set_fr_fr(oCenter, child->centerRealMP, child->centerImagMP, MPC_RNDNN);
                shiftMultipoleExpansionMP(node->multipoleExpansionMP, 
                center, oCenter, 
                node->multipoleExpansionDegree, child->multipoleExpansionMP,
                tempVec, tempLen,
                zVec, zLen, choose, mode);              
            }
        }       
    }

    // clean up mpc variables and free allocated space
    mpc_clear(center);
    mpc_clear(oCenter);
    for (unsigned long int i = 0; i < tempLen; i++){
        mpc_clear(tempVec[i]);
    }
    for (unsigned long int i = 0; i < zLen; i++){
        mpc_clear(zVec[i]);
    }
    free(tempVec);
    free(zVec);
    free(children);
}


void downwardPassMP(struct  TreeNode *root,
unsigned long int local_deg, unsigned long int prec, int mode){

    struct Stack *s;
    struct TreeNode *node, *inter;
    unsigned int d = local_deg;
    d = (local_deg < root->multipoleExpansionDegree)? root->multipoleExpansionDegree: local_deg;

    s = newStack();
    stackPush(s, root);
    unsigned long int tempLen = d+2;
    unsigned long int zLen = d+1;

    mpc_t center, oCenter, shift;
    mpc_t *tempVec = (mpc_t *)malloc(tempLen*sizeof(mpc_t));
    mpc_t *zVec = (mpc_t *)malloc(zLen*sizeof(mpc_t));
    mpz_t choose;
    
    mpc_init2(center, prec);
    mpc_init2(oCenter, prec);
    mpc_init2(shift, prec);
    for (unsigned long int i = 0; i < tempLen; i++){
        mpc_init2(tempVec[i], prec);
    }
    for (unsigned long int i = 0; i < zLen; i++){
        mpc_init2(zVec[i], prec);
    }
    mpz_init2(choose, prec);

    while (!stackIsEmpty(s)){
        node = (struct TreeNode*)stackPop(s);
        mpc_set_fr_fr(center, node->centerRealMP, node->centerImagMP, MPC_RNDNN);

        if (node->ch1 != NULL){
            stackPushLeft(s, node->ch1);
        }
        if (node->ch2 != NULL){
            stackPushLeft(s, node->ch2);
        }
        if (node->ch3 != NULL){
            stackPushLeft(s, node->ch3);
        }
        if (node->ch4 != NULL){
            stackPushLeft(s, node->ch4);
        }        
        node->localExpansionDegree = local_deg;
        node->localExpansionMP = 
        (mpc_t *)malloc((local_deg+1)*sizeof(mpc_t));
        for (int i = 0; i <= local_deg; ++i){
            mpc_init2(node->localExpansionMP[i], prec);
            mpc_set_d(node->localExpansionMP[i], 0, MPC_RNDNN);
        }

        for (unsigned int i = 0; i < node->numInteractions; ++i){
            inter = node->interactions[i];
            mpc_set_fr_fr(oCenter, inter->centerRealMP, inter->centerImagMP, MPC_RNDNN);
            multipoleToLocalMP(node->localExpansionMP, center,
            oCenter, inter->multipoleExpansionDegree, 
            local_deg, inter->multipoleExpansionMP,
            tempVec, tempLen,
            zVec, zLen, choose, mode);
        }

        if(node->parent != NULL){
            inter = node->parent;
            mpc_set_fr_fr(oCenter, inter->centerRealMP, inter->centerImagMP, MPC_RNDNN);
            shiftLocalExpansionMP(node->localExpansionMP, center, oCenter,
            local_deg, inter->localExpansionMP,
            tempVec, tempLen, shift);
        }
    }


    mpc_clear(center);
    mpc_clear(oCenter);
    mpc_clear(shift);
    mpz_clear(choose);
    for (unsigned long int i = 0; i < tempLen; i++){
        mpc_clear(tempVec[i]);
    }
    for (unsigned long int i = 0; i < zLen; i++){
        mpc_clear(zVec[i]);
    }
    free(tempVec);
    free(zVec);
}

void evaluateLEMP(mpc_t *res, 
mpc_t *sources, unsigned long int ns,
mpc_t *vector, struct TreeNode **map,
mpc_t *targets, unsigned long int nt,
struct TreeNode *root, unsigned long int prec, int mode){
    struct TreeNode *node, *temp;
    mpc_t ts, center, shift;
    mpc_init2(ts, prec);
    mpc_init2(center, prec);
    mpc_init2(shift, prec);
    for (unsigned long int t = 0; t < nt; ++t){
        mpc_set_d(res[t], 0, MPC_RNDNN);
        node = map[t];
        mpc_set_fr_fr(center, node->centerRealMP, node->centerImagMP, MPC_RNDNN);
        evalLocalExpansionMP(ts, targets[t], center, node->localExpansionDegree, 
        node->localExpansionMP, shift);
        mpc_add(res[t], res[t], ts, MPC_RNDNN);
        for (int j = 0; j < node->numNearNeighbors; ++j){
            temp = node->nearNeighbors[j];
            if (mode == 0){
                for (unsigned long int k = 0; k < temp->numSources; ++k){
                    unsigned long int idx = temp->sourcesPool[k];
                    mpc_sub(shift, targets[t], sources[idx], MPC_RNDNN);
                    mpc_div(shift, vector[idx], shift, MPC_RNDNN);
                    mpc_add(res[t], res[t], shift, MPC_RNDNN);
                }
            } else {
                for (unsigned long int k = 0; k < temp->numSources; ++k){
                    unsigned long int idx = temp->sourcesPool[k];
                    mpc_sub(shift, targets[t], sources[idx], MPC_RNDNN);
                    mpc_log(shift, shift, MPC_RNDNN);
                    mpc_mul(shift, vector[idx], shift, MPC_RNDNN);
                    mpc_add(res[t], res[t], shift, MPC_RNDNN);
                }
            }

        }
    }
    mpc_clear(ts);
    mpc_clear(center);
    mpc_clear(shift);
}

void fast_multipole_method_MP(mpc_t *res,
mpc_t *sources, mpc_t *targets,
mpc_t *vector, 
unsigned long int ns, unsigned long int nt, 
unsigned int med, unsigned int led, 
unsigned long int mx, int mode, unsigned long int prec){
    struct TreeNode ** map = (struct TreeNode**)malloc(nt*sizeof(struct TreeNode*));
    struct TreeNode *root;
    root = createTreeMP(sources, vector, targets, ns, nt, mx, map, prec);
    upwardPassMP(root, sources, vector, med, mx, prec, mode);
    generateNeighbourInteractionMP(root);
    downwardPassMP(root, led, prec, mode);
    evaluateLEMP(res, sources, ns, vector, map, targets, nt, root, prec, mode);
    freeTree(root);
    free(map);
}