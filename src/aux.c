#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "../include/mp_fmm.h"

/* 
Tree Structures For Partition of FMM
*/



struct TreeNode* newTreeNodeDP(struct TreeNode* par,
    double radiusDouble, double centerRealDouble, double centerImagDouble){
    struct TreeNode * temp = (struct TreeNode*)malloc(sizeof(struct TreeNode));

    temp->radiusDouble = radiusDouble;
    temp->centerRealDouble = centerRealDouble;
    temp->centerImagDouble = centerImagDouble;

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
    temp->multipoleExpansionDouble = NULL;
    temp->localExpansionDouble = NULL;
    return temp;
}

void getTreeNodeChildren(struct TreeNode *node, struct TreeNode **children){
    struct TreeNode * tmp;
    tmp = node;
    for (int i = 0; i < 4; i++){
        children[i] = NULL;
    }
    children[0] = tmp->ch1;
    children[1] = tmp->ch2;
    children[2] = tmp->ch3;
    children[3] = tmp->ch4;
    return;
}


void freeTreeNode(struct TreeNode* node){
    // doest not mean to be used by itself

    if (node->parent) {
        struct TreeNode* par = node->parent;
        if (par->ch1 == node){
            par->ch1 = NULL;
        }
        else if (par->ch2 == node){
            par->ch2 = NULL;
        }
        else if (par->ch3 == node){
            par->ch3 = NULL;
        }
        else if (par->ch4 == node){
            par->ch4 = NULL;
        }
    }

    if (node->ch1) node->ch1->parent = NULL;
    if (node->ch2) node->ch2->parent = NULL;
    if (node->ch3) node->ch3->parent = NULL;
    if (node->ch4) node->ch4->parent = NULL;

    if (node->sourcesPool) free(node->sourcesPool);
    if (node->targetsPool) free(node->targetsPool);
    if (node->nearNeighbors) free(node->nearNeighbors);
    if (node->interactions) free(node->interactions);

    if (node->isDoublePrecision == 0) {
        // if (node->multipoleExpansion != NULL) {
        //     for (unsigned int i = 0; i <= node->multipoleExpansionDegree; ++i){
        //         //mpc_clear(node->multipoleExpansion[i]);
        //         mpf_clear(node->multipoleExpansion[i]);
        //     }
        //     free(node->multipoleExpansion);
        // }
        // if (node->localExpansion != NULL) {
        //     for (unsigned int i = 0; i <= node->localExpansionDegree; ++i){
        //         //mpc_clear(node->multipoleExpansion[i]);
        //         mpf_clear(node->multipoleExpansion[i]);
        //     }
        //     free(node->localExpansion);
        // }
        
    } else {
        if (node->multipoleExpansionDouble) free(node->multipoleExpansionDouble);
        if (node->localExpansionDouble) free(node->localExpansionDouble);
    }
}



void freeTree(struct TreeNode* root){
    //free the entire tree from root using level order traversal
    struct Stack* cur = newStack();
    stackPush(cur, root);
    struct Stack* nxt = newStack();
    struct Stack *tmp;
    struct TreeNode* temp;

    while(!stackIsEmpty(cur)) {
        while(!stackIsEmpty(cur)) {
            temp = (struct TreeNode*)stackPop(cur);
            if (temp->ch1 != NULL){
                stackPush(nxt, temp->ch1);
            }
            if (temp->ch2 != NULL){
                stackPush(nxt, temp->ch2);
            }
            if (temp->ch3 != NULL){
                stackPush(nxt, temp->ch3);
            }
            if (temp->ch4 != NULL){
                stackPush(nxt, temp->ch4);
            }
            freeTreeNode(temp);
            temp = NULL;
        }
        tmp = cur;
        cur = nxt;
        nxt = tmp;
    }
    stackFreeStack(cur);
    stackFreeStack(nxt);
    return;
}

int isNearNeighborDP(struct TreeNode *n1, struct TreeNode *n2){
    /* 
       check if two tree nodes are near neighbors:
       (1) same node
       (2) share an endpoint
       (3) share a side

       return 1 if they are NN,
       return 0 otherwise.
    */
    double tolerance = 2.1 * n1->radiusDouble;
    return (fabs(n1->centerRealDouble - n2->centerRealDouble) <= tolerance && 
    fabs(n1->centerImagDouble - n2->centerImagDouble) <= tolerance);
}

void nodeNeighborInteractionDP(struct TreeNode *node){
    /*
    find nearby nodes on the same level and classify them
    into either Near Neighbor or Interaction
    */
    int numCandidates = 0;
    struct TreeNode * c1, * c2, * c3, * c4;
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
            if (isNearNeighborDP(tmp, node) == 1){
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
}

void generateNeighbourInteractionDP(struct TreeNode * root){
    struct Stack * stack = newStack();
    struct TreeNode* node;
    stackPush(stack, root);
    while (!stackIsEmpty(stack)){
        node = (struct TreeNode*)stackPop(stack);
        nodeNeighborInteractionDP(node);
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
    stack = NULL;
    return;
}

void upwardPassDP(struct TreeNode *root, 
    double complex * sources, 
    double complex * vector,
    unsigned int multi_deg,
    int maxNumSources,
    int m){
    /* 
        m is used as a selector with 
        0: fmm for sums of reciprocals
        1: fmm for sums of complex logs
    */

    double complex center, oCenter;
    struct Stack *s, *rev;
    s = newStack();
    rev = newStack();
    /* 
    create a reverse list of the tree nodes from bottom (with no leaves) to root
    compute the initial expansion at non-empty leaves
    */

    stackPush(s, root);
    struct TreeNode *node, *child;


    double complex* tempS = (double complex *)malloc(maxNumSources * sizeof(double complex));
    /*
    tempS temporary space for required computation in multipole expansion
    */
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
            center = CMPLX(node->centerRealDouble, node->centerImagDouble);
            node->multipoleExpansionDegree = multi_deg;
            node->multipoleExpansionDouble = 
            (double complex *)malloc((multi_deg + 1)*sizeof(double complex));
            for (int i = 0; i <= multi_deg; i++){
                node->multipoleExpansionDouble[i] = 0.0 + 0.0*I;
            }
            multipoleExpansionDP(node->multipoleExpansionDouble, center, 
            multi_deg, node->numSources,
            node->sourcesPool, sources, vector, m, tempS);
        }
    }
    free(tempS);
    tempS = NULL;
    
    double complex *temp = (double complex*)malloc((multi_deg+1)*sizeof(double complex));
    double complex *z0   = (double complex*)malloc((multi_deg+1)*sizeof(double complex));
    struct TreeNode ** children = (struct TreeNode**)malloc(4*sizeof(struct TreeNode*));

    while(!stackIsEmpty(rev)){
        node = (struct TreeNode*)stackPop(rev);
        center = CMPLX(node->centerRealDouble, node->centerImagDouble);
        node->multipoleExpansionDegree = multi_deg;
        node->multipoleExpansionDouble = 
        (double complex *)malloc((multi_deg + 1)*sizeof(double complex));
        for (int i = 0; i <= multi_deg; i++){
                node->multipoleExpansionDouble[i] = 0.0 + 0.0*I;
        }
        getTreeNodeChildren(node, children);
        for (int i = 0; i < 4; i++){
            child = children[i];
            if (child != NULL){
                oCenter = CMPLX(child->centerRealDouble, child->centerImagDouble);
                shiftMultipoleExpansionDP(node->multipoleExpansionDouble, 
                center, oCenter, multi_deg, 
                child->multipoleExpansionDouble, temp, z0, m);              
            }
        }       
    }
    free(children);
    free(temp);
    free(z0);
    stackFreeStack(s);
    stackFreeStack(rev);
    return;
}

void downwardPassDP(
struct TreeNode *root,
unsigned int local_deg,
int m // m = 0 -> sum of reciprocals
){
    struct Stack *s;
    struct TreeNode *node, *inter;
    unsigned int d;
    if (local_deg >= root->multipoleExpansionDegree){
        d = local_deg;
    } else {
        d = root->multipoleExpansionDegree;
    }

    s = newStack();
    stackPush(s, root);

    double complex *temp = (double complex*)malloc((local_deg+1)*sizeof(double complex));
    double complex *z    = (double complex*)malloc((d+1)*sizeof(double complex));

    while (!stackIsEmpty(s)){
        node = (struct TreeNode*)stackPop(s);

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
        node->localExpansionDouble = 
        (double complex*)malloc((local_deg+1)*sizeof(double complex));
        for (int i = 0; i <= local_deg; ++i){
            node->localExpansionDouble[i] = CMPLX(0, 0);
        }

        for (unsigned int i = 0; i < node->numInteractions; ++i){
            inter = node->interactions[i];
            multipoleToLocalDP(node->localExpansionDouble,
            CMPLX(node->centerRealDouble, node->centerImagDouble),
            CMPLX(inter->centerRealDouble, inter->centerImagDouble),
            inter->multipoleExpansionDegree, 
            node->localExpansionDegree,
            inter->multipoleExpansionDouble,
            temp, z, m);
        }

        if(node->parent != NULL){
            inter = node->parent;
            shiftLocalExpansionDP(node->localExpansionDouble,
            CMPLX(node->centerRealDouble, node->centerImagDouble),
            CMPLX(inter->centerRealDouble, inter->centerImagDouble),
            local_deg,
            inter->localExpansionDouble,
            temp);
        }
    }
    free(temp);
    free(z);
    stackFreeStack(s);
    return;
}


struct TreeNode * createTreeDP(
double complex *sourcesD,
double complex *vectorD,
double complex *targetsD,
unsigned long int numSrc,
unsigned long int numTgt,
unsigned long int maxs,
struct TreeNode ** nodeMap
){
    // adaptive partitioning    
    double min_x = creal(sourcesD[0]);
    double max_x = creal(sourcesD[0]);
    double min_y = cimag(sourcesD[0]);
    double max_y = cimag(sourcesD[0]);
    double x, y, radius;

    for (unsigned long int i = 0; i < numSrc; ++i){
        x = creal(sourcesD[i]);
        y = cimag(sourcesD[i]);
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }
    for (unsigned long int i = 0; i < numTgt; ++i){
        x = creal(targetsD[i]);
        y = cimag(targetsD[i]);
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }

    // compute radius of the bounding box
    // and its center
    x = (min_x + max_x)/2.0;
    y = (min_y + max_y)/2.0;
    radius = max_x - min_x;
    if (radius < max_y - min_y) radius = max_y - min_y;
    radius *= 0.50; 
    struct Stack *cur, *nxt, *tmp;
    struct TreeNode *node, *ch1, * ch2, *ch3, *ch4;


    struct TreeNode *root;
    root = newTreeNodeDP(NULL, radius, x, y);
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
                radius = 0.5*node->radiusDouble;
                ch1 = newTreeNodeDP( node, radius, 
                node->centerRealDouble + radius, 
                node->centerImagDouble + radius);
                node->ch1 = ch1;
                ch2 = newTreeNodeDP(node, radius, 
                node->centerRealDouble + radius, 
                node->centerImagDouble - radius);
                node->ch2 = ch2;
                ch3 = newTreeNodeDP(node, radius, 
                node->centerRealDouble - radius, 
                node->centerImagDouble + radius);
                node->ch3 = ch3;
                ch4 = newTreeNodeDP(node, radius, 
                node->centerRealDouble - radius, 
                node->centerImagDouble - radius);
                node->ch4 = ch4;
                stackPush(nxt, ch1);
                stackPush(nxt, ch2);
                stackPush(nxt, ch3);
                stackPush(nxt, ch4);
                subdivideDP(node, ch1, ch2, ch3, ch4, sourcesD, targetsD);
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
    }
    return root;
}

void subdivideDP(struct TreeNode *node, struct TreeNode *ch1, 
struct TreeNode *ch2, struct TreeNode *ch3, struct TreeNode *ch4, 
double complex * sources, double complex *targets){
    // wipe the source pool of node and distribute the sources to children nodes
    ch1->sourcesPool = (unsigned long int *)malloc(node->numSources*sizeof(unsigned long int));
    ch2->sourcesPool = (unsigned long int *)malloc(node->numSources*sizeof(unsigned long int));
    ch3->sourcesPool = (unsigned long int *)malloc(node->numSources*sizeof(unsigned long int));
    ch4->sourcesPool = (unsigned long int *)malloc(node->numSources*sizeof(unsigned long int));

    for (unsigned long int i = 0; i < node->numSources; i++){
        unsigned long int idx = node->sourcesPool[i];
        double real, imag;
        real = creal(sources[idx]);
        imag = cimag(sources[idx]);
        if(real > node->centerRealDouble) {
            if (imag > node->centerImagDouble) {
                // put source in ch1
                ch1->sourcesPool[ch1->numSources] = idx;
                ch1->numSources += 1;
            } else {
                // put source in ch2
                ch2->sourcesPool[ch2->numSources] = idx;
                ch2->numSources += 1;
            } 
        } else {
            if (imag > node->centerImagDouble) {
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
        double real, imag;
        real = creal(targets[idx]);
        imag = cimag(targets[idx]);
        if(real > node->centerRealDouble) {
            if (imag > node->centerImagDouble) {
                // put target in ch1
                ch1->targetsPool[ch1->numTargets] = idx;
                ch1->numTargets += 1;
            } else {
                // put target in ch2
                ch2->targetsPool[ch2->numTargets] = idx;
                ch2->numTargets += 1;
            } 
        } else {
            if (imag > node->centerImagDouble) {
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


void evaluateDP(double complex *res, 
double complex *sources, unsigned long int ns,
double complex *vector,
double complex *targets, unsigned long int nt,
struct TreeNode *root, int m){
    double complex z;
    double real, imag;
    struct TreeNode *cur, *nxt;
    for (unsigned long int t = 0; t < nt; ++t){
        z = targets[t];
        real = creal(z);
        imag = cimag(z);
        res[t] = 0;
        nxt = root;
        while (nxt) {
            cur = nxt;
            for( int j = 0; j < cur->numInteractions; ++j){
                nxt = cur->interactions[j];
                res[t] += evalMultipoleExpansionDP(z, 
                CMPLX(nxt->centerRealDouble, nxt->centerImagDouble),
                nxt->multipoleExpansionDegree, nxt->multipoleExpansionDouble, m);
            }
            if(real > cur->centerRealDouble) {
                if (imag > cur->centerImagDouble) {
                    nxt = cur->ch1;
                } else {
                    nxt = cur->ch2;
                } 
            } else {
                if (imag > cur->centerImagDouble) {
                    nxt = cur->ch3;
                } else {
                    nxt = cur->ch4;
                }
            }
        }
        // cur is the leave node
        for (int j = 0; j < cur->numNearNeighbors; ++j){
            nxt = cur->nearNeighbors[j];
            for (int k = 0; k < nxt->numSources; ++k){
                unsigned long int idx = nxt->sourcesPool[k];
                res[t] += (m==0) ? vector[idx]/(z - sources[idx]): 
                vector[idx] * clog(z - sources[idx]);
            }
        }
    }
    return;
}

void evaluateLEDP(double complex *res, 
double complex *sources, unsigned long int ns,
double complex *vector, struct TreeNode **map,
double complex *targets, unsigned long int nt,
struct TreeNode *root, int m) {
    struct TreeNode *node, *temp;
    for (unsigned long int t = 0; t < nt; ++t){
        res[t] = 0;
        node = map[t];
        res[t] += evalLocalExpansionDP(targets[t], 
        CMPLX(node->centerRealDouble, node->centerImagDouble), 
        node->localExpansionDegree, node->localExpansionDouble);
        for (int j = 0; j < node->numNearNeighbors; ++j){
            temp = node->nearNeighbors[j];
            for (int k = 0; k < temp->numSources; ++k){
                unsigned long int idx = temp->sourcesPool[k];
                res[t] += (m==0) ? vector[idx]/(targets[t] - sources[idx]): 
                vector[idx] * clog(targets[t] - sources[idx]);
            }
        }
    }
}

void fast_multipole_method_DP(double complex *res,
double complex *sources, double complex *targets,
double complex *vector, 
unsigned long int ns, unsigned long int nt, 
unsigned int med, unsigned int led, 
unsigned long int mx, int mode){
    struct TreeNode ** map = (struct TreeNode**)malloc(nt*sizeof(struct TreeNode*));
    struct TreeNode *root;
    root = createTreeDP(sources, vector, targets, ns, nt, mx, map);
    upwardPassDP(root, sources, vector, med, mx, mode);
    generateNeighbourInteractionDP(root);
    downwardPassDP(root, led, mode);
    evaluateLEDP(res, sources, ns, vector, map, targets, nt, root, mode);
    freeTree(root);
    free(map);
}


/* --- Stack Functionality --- */

struct StackNode {
    void * key;
    struct StackNode* prev;
    struct StackNode* next;
};

struct Stack {
    struct StackNode * head, *tail;
};

struct StackNode* newStackNode(void* node) {
    struct StackNode * temp = (struct StackNode*)malloc(sizeof(struct StackNode));
    temp->key = node;
    temp->prev = NULL;
    temp->next = NULL;
    return temp;
}

struct Stack* newStack(){
    struct Stack* s = (struct Stack*)malloc(sizeof(struct Stack));
    s->head = s->tail = NULL;
    return s;
}

int stackIsEmpty(struct Stack* s){
    return s->head == NULL;
}

void stackPush(struct Stack* s, void* t){
    struct StackNode* temp = newStackNode(t);
    if (s->tail == NULL) {
        s->head = s->tail = temp;
        return;
    }

    s->tail->next = temp;
    temp->prev = s->tail;
    s->tail = temp;
    return;
}

void stackPushLeft(struct Stack *s, void* t){
    struct StackNode* temp = newStackNode(t);
    if (s->head == NULL){
        s->head = s->tail = temp;
        return;
    }

    s->head->prev = temp;
    temp->next = s->head;
    s->head = temp;
    return;
}

void* stackPop(struct Stack* s){
    if (stackIsEmpty(s)) {
        return NULL;
    } 

    void * temp = s->tail->key;
    if (s->head == s->tail) { // last node
        free(s->tail);
        s->head = s->tail = NULL;
    } else {
        s->tail = s->tail->prev;
        free(s->tail->next);
        s->tail->next = NULL;
    }
    return temp;
}

void* stackPopLeft(struct Stack* s){
    if (stackIsEmpty(s) ){
        return NULL;
    }
    void* temp = s->head->key;
    if (s->head == s->tail) {
        free(s->tail);
        s->head = s->tail = NULL;
    } else {
        s->head = s->head->next;
        free(s->head->prev);
        s->head->prev = NULL;
    }
    return temp;
}

void stackFreeStack(struct Stack* s){
    if (s == NULL){
        return;
    }
    while (stackIsEmpty(s) == 0){ // while not empty
        stackPop(s);
    }
    free(s);
    return;
}


/* reading complex numbers from file */


/* Read all double complex numbers from stream 'input',
   into a dynamically allocated array.
   (Similar to POSIX.1 getline(), but with double vectors.)
   If *dataptr is not NULL, and *sizeptr > 0,
   it will initially be used (but reallocated if needed).
   Returns the number of vectors read,
   or 0 with errno set if an error occurs.
*/
unsigned long int complexDPReadall(FILE *input, double complex **dataptr, unsigned long int *sizeptr)
{
    double complex *data;
    double tempReal, tempImag;
    unsigned long int size;
    unsigned long int used = 0;

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
        }

        if (fscanf(input, " %lf %lf", &tempReal, &tempImag) != 2)
            break;
        data[used] = CMPLX(tempReal, tempImag);
        /* One more vector read successfully. */
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

unsigned long int complexDPWriteall(FILE* output, double complex * data, unsigned long int n) {
    unsigned long int used = 0;
    for (unsigned long int i = 0; i < n; i++){
        if (fprintf(output, "%.17g %.17g\n", creal(data[i]), cimag(data[i])) <= 0){
            return 0;
        }
        used++;
    }
    return used;
}