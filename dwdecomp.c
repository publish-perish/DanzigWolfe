#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <ilcplex/cplexcheck.h>
#include <ilcplex/cplex.h>


#include "dwdecomp.h"

/*   minimize  c*x
 *   subject to  Ax = b
 *               Hx = d
 *               l <= x <= u
 *
 *  Treat the constraints with A as the complicating constraints, and
 *  the constraints with H as the subproblems.
 *
 */

#define SENSE -1

int DEBUG = 0;
int VERBOSE = 0;

CPXENVptr     env;
CPXLPptr      slp[MAXSUBS];
CPXLPptr      mlp;
CPXLPptr      rmp;

int status;
int flag;
char fname[256], sname[256];
int i, j, m, n, l, q, p, idx;
int ROWSMAST;
int ROWSSUB[MAXSUBS];

int build(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest);

int branch(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest, double *soln);

double solvebranch(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest, double *soln);

int debug(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest, char* dfname);

int main(int argc, char* argv[])
{
    int NNODES = 3, NSPECTRUMS = 1, NSESSIONS = 1;
    if(argc >= 4) {
        NNODES = atoi(argv[1]);
        NSPECTRUMS = atoi(argv[2]); 
        NSESSIONS = atoi(argv[3]);
    }
     
    int status;
    int i, j, sarray[2*NSESSIONS];
    int *links[MAXSUBS], *requests[NNODES];
    int *src, *dest;
    double soln[MAXEXTREMES];

    clock_t begin, end;
    double time_spent;

    for(i=0; i<MAXEXTREMES; i++) soln[i] = 0.0;

    // SET UP RANDOM NETWORK NODES
    for(i=0; i<MAXSUBS; i++) links[i] = (int*)calloc(MAXSUBS, sizeof(int));
    for(i=0; i<NNODES; i++) requests[i] = (int*)calloc(NNODES, sizeof(int));

    for(i=0; i<MAXSUBS; i++) {
        for(j=0; j<MAXSUBS; j++) {
            if(rand() % 100 > 60 && (i != j))
                links[i][j] = 1;
            else
                links[i][j] = 0;
        }
    }
    // SET UP RANDOM SESSION REQUESTS
    for(i=0; i<2*NSESSIONS; i++) {
        if(i<NNODES)
            sarray[i] = i;
        else
            sarray[i] = 0;
    }
	shuffle(sarray, 2*NSESSIONS, NSESSIONS);
    src = &sarray[0];
    dest = &sarray[NSESSIONS];
    printvecd(sarray, 0, 2*NSESSIONS);


    if(argc == 6) {
        printf("Setting session request...\n");
        sarray[0]  = atoi(argv[4]) - 1;
        sarray[NSESSIONS]  = atoi(argv[5]) - 1;
    }   
    // START CPLEX
    env = CPXopenCPLEX (&status);
    if ( env == NULL ) {
        char  errmsg[1024];
        fprintf (stderr, "Could not open CPLEX environment.\n");
        CPXgeterrorstring (env, status, errmsg);
        fprintf (stderr, "%s", errmsg);
        goto TERMINATE;
    }
    status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
    status = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_ON);
    status = CPXsetintparam (env, CPX_PARAM_LPMETHOD, CPX_ALG_DUAL);
    status = CPXsetintparam (env, CPXPARAM_SolutionType, CPX_BASIC_SOLN);

    printf("Problem has %d nodes, %d spectrums and %d sessions.\n\n", NNODES, NSPECTRUMS, NSESSIONS);
    printf("Transmission Ranges:\n");
    printmatd(links, 0, NNODES, 0, NNODES);
    printf("Session Requests:\n");
    for(i=0; i<NSESSIONS; i++) {
        printf("f%d : ", i+1); printf("%d ---> %d\n", src[i]+1, dest[i]+1);
    }
    printf("\n");

    printmatd(requests, 0, NNODES, 0, NNODES);
    begin = clock();
    build(VERBOSE, DEBUG, NNODES, NSPECTRUMS, NSESSIONS, links, src, dest);
    branch(VERBOSE, DEBUG, NNODES, NSPECTRUMS, NSESSIONS, links, src, dest, soln);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Run time = %f\n", time_spent);

    // Print Solution
    printf("\n\n"); printf("%*s\n", 70,"-"); printf("\n\n");
    for(i=0; i<NSESSIONS; i++)
        requests[src[i]][dest[i]] = 1;



    // PRINT SOLUTION
    printf("Solution:\n");
    checksolnmat(soln, NNODES, NSPECTRUMS, NSESSIONS);
    printf("\n\n"); printf("%*s\n", 70,"-"); printf("\n\n");

TERMINATE:
    if ( status ) fprintf (stderr,"Exit status was %d\n", status); 
    for(int n=0; n<NNODES; n++) CPXfreeprob (env, &slp[n]);
    CPXfreeprob (env, &mlp);
    CPXfreeprob (env, &rmp);
    if ( env != NULL ) {
        status = CPXcloseCPLEX (&env);
        if ( status ) {
            char  errmsg[1024];
            fprintf (stderr, "Could not close CPLEX environment.\n");
            CPXgeterrorstring (env, status, errmsg);
            fprintf (stderr, "%s", errmsg);
        }
    }
 
    for(i=0; i<MAXSUBS; i++) free(links[i]);
    for(n=0; n<NNODES; n++)  free(requests[n]); 
    
    return 0;
}

void namecols(CPXENVptr env, CPXLPptr lp, int NSESSIONS, int NSPECTRUMS, int NNODES)
{
    char varname[256];
    int status;
    for(int m=0; m<NSPECTRUMS; m++) {
        for(int i=0; i<NNODES; i++) {
            for(int j=0; j<NNODES; j++) {
                sprintf(varname, "m%ds%d%d", m+1, i+1, j+1);
                //printf("idx=%d ", i*NNODES+j);
                status = CPXchgname(env, lp, 'c', NNODES*NNODES*m+i*NNODES+j, varname);
                //printf("%s ", varname);
            }
        }
    }
    for(int l=0; l<NSESSIONS; l++){
        for(int i=0; i<NNODES; i++) {
            for(int j=0; j<NNODES; j++) {
                sprintf(varname, "l%df%d%d", l+1, i+1, j+1);
                //printf("idx=%d ",NNODES*NNODES+NNODES*NNODES*l+i*NNODES+j); 
                status = CPXchgname(env, lp, 'c', NNODES*NNODES*NSPECTRUMS+NNODES*NNODES*l+i*NNODES+j, varname);
                //printf("%s ", varname);
            }
        }
    }
    for(int l=0; l<NSESSIONS; l++){
        sprintf(varname, "rlm%d", l+1);
        //printf("idx=%d ", NNODES*NNODES+NNODES*NNODES*NSESSIONS+l);
        status = CPXchgname(env, lp, 'c', NNODES*NNODES*NSPECTRUMS+NNODES*NNODES*NSESSIONS+l, varname);
        //printf("%s ", varname);
    }
}
