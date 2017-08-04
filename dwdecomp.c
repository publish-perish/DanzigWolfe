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

double solve(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest, double *soln);

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
        //goto TERMINATE;
    }
    status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
    status = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_OFF);
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
    for(i=0; i<NNODES; i++) for(j=0; j<NNODES; j++) requests[i][j] = 0;
    for(i=0; i<NSESSIONS; i++)
        requests[src[i]][dest[i]] = 1;

    printmatd(requests, 0, NNODES, 0, NNODES);
    begin = clock();
    //build(VERBOSE, DEBUG, NNODES, NSPECTRUMS, NSESSIONS, links, src, dest);
    int NNNODES = NNODES*NNODES;
    int ROWSSUB[NNODES];
    int NCOLS = NSPECTRUMS*NNNODES+NSESSIONS*NNNODES+NSESSIONS;
    double coef = 0;
    char sense = 'E';
    ROWSMAST = 0;

    for(i=0; i<NNODES; i++) ROWSSUB[i] = 0;
    
    // CREATE LPs
    rmp = CPXcreateprob (env, &status, "rmp");
    mlp = CPXcreateprob (env, &status, "mlp");
    for(n=0; n<NNODES; n++) {
        sprintf(sname, "%s%d", "subprob", n+1);
        slp[n] = CPXcreateprob (env, &status, "sub");
    }
    // CREATE CPLEX SUBPROBLEMS
    for(n=0; n<NNODES; n++) {
        status = CPXnewcols(env, slp[n], NCOLS, NULL, NULL, NULL, NULL, NULL);
        namecols(env, slp[n], NSESSIONS, NSPECTRUMS, NNODES);
        status = CPXnewrows(env, slp[n], MAXROWS, NULL, NULL, NULL, NULL);
        status = CPXchgobjsen(env, slp[n], -1);
    }
    //-----------------------------CREATE SUBPROBLEM MATRIX---------------------------------//
    printf("Creating %dx%d subproblem matrix...\n", MAXROWS, NCOLS);
    // Routing (4)
    printf("Routing (4)\t: c%d - ", ROWSSUB[0]+1);
    for(i=0; i<NNODES; i++) {
        for(j=0; j<NNODES; j++) 
            if(links[i][j]) {
            for(l=0; l<NSESSIONS; l++) {
                    status = CPXchgcoef(env, slp[i], ROWSSUB[i], NSPECTRUMS*NNNODES+NNNODES*l+NNODES*i+j, -1);
            }
            for(m=0; m<NSPECTRUMS; m++) {
                    status = CPXchgcoef(env, slp[i], ROWSSUB[i], NNNODES*m+NNODES*i+j, 1);
            }
            status = CPXchgcoef(env, slp[i], ROWSSUB[i], -1, 0); 
            sense = 'E';
            status = CPXchgsense(env, slp[i], 1, &ROWSSUB[i], &sense);
            ROWSSUB[i] += 1;
            }
    }
    printf("c%d\n", ROWSSUB[0]);
    // Routing (5)        
    printf("Routing (5)\t: c%d - ", ROWSSUB[0]+1);
    for(i=0; i<NNODES; i++) {
        for(l=0; l<NSESSIONS; l++) 
            if(i == dest[l]) {
                for(j=0; j<NNODES; j++) 
                    if(links[i][j]) {
                        status = CPXchgcoef(env, slp[i], ROWSSUB[i], NSPECTRUMS*NNNODES+NNNODES*l+NNODES*i+j, 1);
                    }
                status = CPXchgcoef(env, slp[i], ROWSSUB[i], -1, 0); 
                sense = 'E';
                status = CPXchgsense(env, slp[i], 1, &ROWSSUB[i], &sense);
                ROWSSUB[i] += 1;
            }
    }
    double ub = 1.0;
    char lu = 'U';
    for(n=0; n<NNODES; n++) 
        for(int m=0; m<NSPECTRUMS; m++) 
            for(i=0; i<NNODES; i++) 
               for(j=0; j<NNODES; j++) {
                    idx = NNNODES*m+NNODES*i+j;
                    status = CPXchgbds (env, slp[n], 1, &idx, &lu, &ub);
                }
    if(VERBOSE) 
        for(n=0; n<NNODES; n++) 
            printf("\nSubproblem %d matrix has %d rows and %d cols.\n", n+1, ROWSSUB[n], NCOLS);
       
    // WRITE OUT SUBS
    if(DEBUG)
        for(n=0; n<NNODES; n++) {
            status = CPXdelrows (env, slp[n], ROWSSUB[n], MAXROWS-1);
            sprintf(fname, "%s%d%s", "orig_slp", n+1, ".lp");
            status = CPXwriteprob(env, slp[n], fname, NULL);
        }


    // CREATE MASTER PROBLEM
    status = CPXnewcols(env, mlp, NCOLS, NULL, NULL, NULL, NULL, NULL);
    namecols(env, mlp, NSESSIONS, NSPECTRUMS, NNODES);
    status = CPXnewrows(env, mlp, MAXROWS, NULL, NULL, NULL, NULL);
    status = CPXchgobjsen(env, mlp, -1);
    //-----------------------------CREATE CONSTRAINT MATRIX FOR RMP---------------------------------//
    printf("\nCreating master matrix...\n");
    // Routing (1)
    // TODO: FIX!
    // Self-interference (2)
    printf("Self-interference (2)   : c%d - ", ROWSMAST+1);
    for(m=0; m<NSPECTRUMS; m++) {
        for(i=0; i<NNODES; i++) {
            for(j=0; j<NNODES; j++) {
                if(links[j][i] == 1) {
                    status = CPXchgcoef(env, mlp, ROWSMAST, NNNODES*m+NNODES*j+i, 1);
                    for(q=0; q<NNODES; q++) {
                        if(links[i][q])
                            status = CPXchgcoef(env, mlp, ROWSMAST, NNNODES*m+NNODES*i+q, 1);
                    }
                }
            }
            status = CPXchgcoef(env, mlp, ROWSMAST, -1, 1); 
            sense = 'L';
            status = CPXchgsense(env, mlp, 1, &ROWSMAST, &sense);
            ROWSMAST += 1;
        }
    }
    printf("c%d\n", ROWSMAST);
    // General interference (3)
    printf("General interference (3): c%d - ", ROWSMAST+1);
    for(m=0; m<NSPECTRUMS; m++) {
        for(p=0; p<NNODES; p++) 
            for(q=0; q<NNODES; q++) 
                if(links[p][q] == 1) { // q in T_p i.e. exists link p->q
                    status = CPXchgcoef(env, mlp, ROWSMAST, NNNODES*m+NNODES*p+q, 1); // s_pq^m
                    for(i=0; i<NNODES; i++) { // q \in I_i i.e. exists i->q
                        if(links[i][q]){
                            for(j=0; j<NNODES; j++)
                                if(links[i][j]) // j in T_i i.e. exists link i->j
                                    status = CPXchgcoef(env, mlp, ROWSMAST, NNNODES*m+NNODES*i+j, 1);
                    }}
                    status = CPXchgcoef(env, mlp, ROWSMAST, -1, 1); 
                    sense = 'L';
                    status = CPXchgsense(env, mlp, 1, &ROWSMAST, &sense);
                    ROWSMAST += 1;
                }
    }
    printf("c%d\n", ROWSMAST);
    // Routing (5)        
    printf("Routing (5)             : c%d - ", ROWSMAST+1);
    // TODO: FIX!
    for(l=0; l<NSESSIONS; l++) 
        for(i=0; i<NNODES; i++)
            if(i == src[l]) {
                for(j=0; j<NNODES; j++) 
                    if(links[i][j])
                        status = CPXchgcoef(env, mlp, ROWSMAST, NSPECTRUMS*NNNODES+NNNODES*l+NNODES*j+i, 1);
                status = CPXchgcoef(env, mlp, ROWSMAST, -1, 0); 
                sense = 'E';
                status = CPXchgsense(env, mlp, 1, &ROWSMAST, &sense);
                ROWSMAST += 1;
            }
    printf("c%d\n", ROWSMAST);
    // Data rate (6)        
    printf("Data rate (6)           : c%d - ", ROWSMAST+1);
    for(l=0; l<NSESSIONS; l++) {
        if(links[src[l]][dest[l]]) {
                status = CPXchgcoef(env, mlp, ROWSMAST, NSPECTRUMS*NNNODES+NNNODES*l+NNODES*src[l]+dest[l], 1);
        }
        for(j=0; j<NNODES; j++) {
            if(links[src[l]][j]) {
                status = CPXchgcoef(env, mlp, ROWSMAST, NSPECTRUMS*NNNODES+NNNODES*l+NNODES*src[l]+j, 1);
            }
            if(links[j][dest[l]]) {
                status = CPXchgcoef(env, mlp, ROWSMAST, NSPECTRUMS*NNNODES+NNNODES*l+NNODES*j+dest[l], 1);
            }
        }
        status = CPXchgcoef(env, mlp, ROWSMAST, NSPECTRUMS*NNNODES+NSESSIONS*NNNODES+l, -1);
        status = CPXchgcoef(env, mlp, ROWSMAST, -1, 0); 
        sense = 'E';
        status = CPXchgsense(env, mlp, 1, &ROWSMAST, &sense);
        ROWSMAST += 1;
    }
    printf("c%d\n", ROWSMAST);
    // Relay nodes (7)
    printf("Relay nodes (7)         : c%d - ", ROWSMAST+1);
    for(l=0; l<NSESSIONS; l++) {
        for(i=0; i<NNODES; i++) {
            if(i != src[l] && i != dest[l]) { 
                for(j=0; j<NNODES; j++) {
                    if(links[i][j] ) { 
                        status = CPXchgcoef(env, mlp, ROWSMAST, NSPECTRUMS*NNNODES+NNNODES*l+NNODES*i+j, 1);
                        for(p=0; p<NNODES; p++)
                            if(links[p][i]) {
                                status = CPXchgcoef(env, mlp, ROWSMAST, NSPECTRUMS*NNNODES+NNNODES*l+NNODES*p+i, -1);
                            }
                    }
                }
            }
        status = CPXchgcoef(env, mlp, ROWSMAST, -1, 0); 
        sense = 'E';
        status = CPXchgsense(env, mlp, 1, &ROWSMAST, &sense);
        ROWSMAST += 1; 
        }
    }
    printf("c%d\n", ROWSMAST);
    for(int m=0; m<NSPECTRUMS; m++) 
        for(i=0; i<NNODES; i++) 
            for(j=0; j<NNODES; j++) {
                idx = NNNODES*m+NNODES*i+j;
                status = CPXchgbds (env, mlp, 1, &idx, &lu, &ub);
            }
    for(i=NCOLS-NSESSIONS; i<NCOLS; i++) {             
        status = CPXchgcoef(env, mlp, -1, i, 1.0);
    }
    status = CPXdelrows (env, mlp, ROWSMAST, MAXROWS-1);

    if(VERBOSE) 
        printf("\nConstraint matrix has %d rows and %d cols.\n", ROWSMAST, NCOLS);
 
    // WRITE OUT MASTER
    if(DEBUG)
        status = CPXwriteprob(env, mlp, "orig_mlp.lp", NULL);


    // CREATE CPLEX RMP
    status = CPXchgobjsen(env, rmp, -1);
    status = CPXnewcols(env, rmp, NSESSIONS+NNODES, NULL, NULL, NULL, NULL, NULL);
    status = CPXnewrows (env, rmp, ROWSMAST+NNODES, NULL, NULL, NULL, NULL);  
    //-----------------------------FORMULATE THE RMP---------------------------------//
    printf("\nStarting RMP formulation...\n");
    // Set the lambdas
    for(i=0; i<ROWSMAST; i++)
        for(j=0; j<NSESSIONS; j++) {
            status = CPXgetcoef (env, mlp, i, NCOLS-NSESSIONS+j, &coef);
            status = CPXchgcoef(env, rmp, i, NNODES+j, coef);
        }
    for(i=ROWSMAST; i<ROWSMAST+NNODES; i++) {
            status = CPXchgcoef(env, rmp, i, -1, 1.0); 
        for(j=0; j<NNODES; j++) {
            if(i == j+ROWSMAST) {
                status = CPXchgcoef(env, rmp, i, j, 1.0);
            }
        }
    }
    for(i=NNODES; i<NSESSIONS+NNODES; i++) { 
        status = CPXchgcoef(env, rmp, -1, i, 1.0);
    }
    for(i=0; i<ROWSMAST; i++) {
        status = CPXgetsense (env, mlp, &sense, i, i);
        status = CPXchgsense(env, rmp, 1, &i,  &sense);
    }
    for(i=0; i<ROWSMAST; i++) {
        status = CPXgetcoef(env, mlp, i, -1, &coef);
        status = CPXchgcoef(env, rmp, i, -1, coef);
    }
    // WRITE OUT RMP
    if (DEBUG) 
        status = CPXwriteprob(env, rmp, "orig_rmp.lp", NULL);

    printf("Problem has %d nodes, %d spectrums and %d sessions.\n\n", NNODES, NSPECTRUMS, NSESSIONS);


    //branch(VERBOSE, DEBUG, NNODES, NSPECTRUMS, NSESSIONS, links, src, dest, soln);
    FILE *fp;
    int NNEXTREMES = 0, maxidx, np;
    int lpstat, found, iter = 0;
    double MAXVAL = 0, MINVAL = CPX_INFBOUND;
    double newcol[ROWSMAST+NSPECTRUMS];
    double objval[NNODES+1];
    double pi[ROWSMAST+NNODES];
    double x[NNODES+1][NCOLS];
    double zeros[NCOLS];
    int     cstat[MAXEXTREMES];
printf("STARTING ITERATION\n");
    printf("\n\n"); printf("%*s\n", 20,"-");  printf("%*s\n", 20,"-"); printf("\n\n");
    
    fp = fopen("extremepts.out", "w");
    for(j=0; j<NCOLS; j++) { soln[j] = 0.0; zeros[j] = 0;}
    for(n=0; n<NNODES; n++) for(i=0; i<NCOLS; i++) x[n][i] =0.0;


while(1) {

    for(i=0; i<ROWSMAST+NNODES; i++)
        pi[i] = 0.0;
    // SOLVE SUBPROBLEMS
    if(rmp == NULL || env == NULL) {
        printf("%p %p : CPLEX POINTER DIED. Exiting.\n", env, rmp);
        return 0;
    }
    status = CPXoptimize (env, rmp);
    //status = CPXdualopt(env, rmp);
    //if ( rmpstat != CPX_STAT_OPTIMAL ) goto TERMINATE;
    for(i=0; i<NNODES; i++) 
        for(j=0; j<NCOLS; j++) x[i][j] = 0.0;
    status = CPXsolution(env, rmp, &lpstat, &objval[0], x[0], pi, NULL, NULL);
    if(VERBOSE) {
        printf("objval = %g\n", objval[0]);
        printf("x="); printvec(x[0], 0, NSESSIONS+NNEXTREMES+NNODES, NNODES);
        printf("Dual values for constraints:\n"); printvec(pi, 0, ROWSMAST, NNNODES);
        printf("Dual values for lambdas:\n"); printvec(pi, ROWSMAST, ROWSMAST+NNODES, NNNODES);
        printf("NNEXTREMES = %d\n", NNEXTREMES);
        printf("\n");
    }
    iter += 1;
    if(NNEXTREMES >= MAXEXTREMES) {
        printf("\n Out of memory.. terminating\n");
        return 0;
    }

    // UPDATE COST AND RHS OF SUBPROB Ki
    for(n=0; n<NNODES; n++) {
        for(i=0; i<NCOLS; i++) {
            x[n+1][i] = 0;
        }
    }
    if(VERBOSE) {
        printf("UPDATING COSTS OF SUBPROBLEMS...\n");
    }

    for(n=0; n<NNODES; n++) {
        for(m=0; m<NSPECTRUMS; m++) {
            for(j=NNNODES*m+NNODES*n; j<NNNODES*m+NNODES*(n+1); j++) {
                for(i=0; i<ROWSMAST; i++) {
                        status = CPXgetcoef (env, mlp, i, j, &coef);
                        x[n+1][j] -=  pi[i] * coef;
                }
            }
        }
        for(l=0; l<NSESSIONS; l++) {
            for(j=NSPECTRUMS*NNNODES+NNNODES*l+NNODES*n; j<NSPECTRUMS*NNNODES+NNNODES*l+NNODES*(n+1); j++) {
                for(i=0; i<ROWSMAST; i++) {
                        status = CPXgetcoef (env, mlp, i, j, &coef);
                        x[n+1][j] -=  pi[i] * coef;
                }
            }
        }
    }

    // Update costs
    for(n=0; n<NNODES; n++) {
        if(pi[ROWSMAST+n]>0 ) {
             for(i=0; i<NCOLS-NSESSIONS; i++)
                    if(x[n+1][i] > 0)
                    x[n+1][i] -= pi[ROWSMAST+n];
               }
    }
    for(n=0; n<NNODES; n++) {
        if(VERBOSE) {
            printf("lp%d", n+1); printvec(x[n+1], 0, NCOLS, NNNODES);
            printnonzero(x[n+1], NNODES, NSPECTRUMS, NSESSIONS);
        }
        for(j=0; j<NCOLS; j++) {
            status = CPXchgcoef(env, slp[n], -1, j, x[n+1][j]);
        }
    }
 
    if(iter == 1 && DEBUG) {
        printf("Writing out subproblems...\n");
        for(n=0; n<NNODES; n++) {
            sprintf(fname, "%s%d%s", "slp", n+1, ".lp");
            status = CPXwriteprob(env, slp[n], fname, NULL);
        }
        status = CPXwriteprob(env, rmp, "rmp.lp", NULL);
    }
    if(DEBUG) {
        sprintf(fname, "%s%d%s", "rmp", iter, ".lp");
        status = CPXwriteprob(env, rmp, fname, NULL);
    }
   
    MINVAL = CPX_INFBOUND;
    // CHECK FOR OPTIMALITY
    found = 0;

    for(n=0; n<NNODES; n++) {
        if(VERBOSE)
            printf("\nSOLVING SUBPROBLEM %d...\n", n+1);
        if(slp[n] == NULL || env == NULL) {
            printf("%p %p : CPLEX POINTER DIED. Exiting.\n", env, slp[n]);
            return 0;
        }
        status = CPXoptimize (env, slp[n]);
        status = CPXsolution (env, slp[n], &lpstat, &objval[n+1], x[n+1], NULL, NULL, NULL);
            if ( lpstat != CPX_STAT_OPTIMAL ) {
				objval[n+1] = 0;
                printf("\n\t\t\t SOLUTION IS NOT OPTIMAL!\n");
            }
        if(VERBOSE) {
            printf("objval = %g   ", objval[n+1]); //printvec(x[1], 0, NCOLS, NNNODES);
            printnonzero(x[n+1], NNODES, NSPECTRUMS, NSESSIONS);
        }
    }
    
    // Find most positive reduced Cost
    MAXVAL = 1.0e-8; np = 0;
    maxidx = 0;
    for(n=0; n<NNODES; n++) {
        if( (objval[n+1]) >= MAXVAL) { 
            if(VERBOSE)
                printf("Found extreme point on node %d\n", n+1);
            MAXVAL = objval[n+1];
            maxidx = n;
            np += 1;
        }
    }
    if(MAXVAL <= 1.0e-8) {
        printf("\nTerminating...\n");
        break;
    }

	if(VERBOSE)
    	printf("Found %d possible entering extreme points, MAXVAL = %g, choosing node %d\n", np, MAXVAL, maxidx+1);
	found = 1;
	for(i=0; i<ROWSMAST+NNODES; i++) {newcol[i] = 0;}
	// Compute the new column
	for(m=0; m<NSPECTRUMS; m++) {
		for(i=0; i<ROWSMAST; i++) {
			for(j=NNNODES*m+NNODES*maxidx; j<NNNODES*m+NNODES*(maxidx+1); j++) {
					status = CPXgetcoef (env, mlp, i, j, &coef);
					newcol[i] += coef* x[maxidx+1][j];
			}
		}
	}
	for(l=0; l<NSESSIONS; l++) {
		for(i=0; i<ROWSMAST; i++) {
			for(j=NSPECTRUMS*NNNODES+NNNODES*l+NNODES*maxidx; j<NSPECTRUMS*NNNODES+NNNODES*l+NNODES*(maxidx+1); j++) {
					status = CPXgetcoef (env, mlp, i, j, &coef);
					newcol[i] += coef * x[maxidx+1][j];
			}
		}
	}
	newcol[ROWSMAST+maxidx] = 1.0;
	if(VERBOSE) {
		printvec(newcol, 0, ROWSMAST+NNODES, NNODES);
		printf("ADDING NEW COLUMN FOR NODE %d\n\n", maxidx+1);
		printf("NNEXTREMES = %d\n", NNEXTREMES);
	}

    // Delete nonbasic cols
    status = CPXgetbase (env, rmp, cstat, NULL); 
    for(i=0; i<NNODES+NNEXTREMES; i++) {
        if(cstat[i] != 1) {
            //status = CPXdelcols(env, rmp, i, i);
            printf("Deleting col %d\n", i+1);
            //NNEXTREMES -= 1;
        }
    }
            
	status = CPXnewcols(env, rmp, 1, NULL, NULL, NULL, NULL, NULL);
	for(i=0; i<ROWSMAST+NNODES; i++) {
		status = CPXchgcoef(env, rmp, i, NNODES+NNEXTREMES, newcol[i]);
	}
	NNEXTREMES += 1;
	for(i=0; i<NNODES+NNEXTREMES; i++) {
		status = CPXchgcoef(env, rmp, -1, i, 0);
	}
	for(i=NNODES+NNEXTREMES; i<NNODES+NNEXTREMES+NSESSIONS; i++) {
		status = CPXchgcoef(env, rmp, -1, i, 1);
		
	}
	for(i=0; i<ROWSMAST; i++) {
		for(j=0; j<NSESSIONS; j++) {
			status = CPXgetcoef (env, mlp, i, NCOLS-NSESSIONS+j, &coef);
			status = CPXchgcoef(env, rmp, i, NNODES+NNEXTREMES+j, coef);
		}
	}
	status = CPXchgcoef(env, rmp, ROWSMAST+maxidx, -1, 1.0);

	if(VERBOSE) {
		printf("Storing solution for node %d\nx=", maxidx+1);
		printnonzero(x[maxidx+1], NNODES, NSPECTRUMS, NSESSIONS);
		printvec(x[maxidx+1], 0, NCOLS, NNODES);
	}
    writeextreme(fp, x[maxidx+1], NCOLS);

    if(found) {
		if(VERBOSE)
        	printf("Number of extremes = %d, iteration %d\n", NNEXTREMES, iter+1);
    }
    else {
        writeextreme(fp, zeros, NCOLS);
        printf("\nTerminating...\n");
        break;
    }
}
    fclose(fp);
    fp = fopen("extremepts.out", "r");
    // PRINT SOLUTION
	if(VERBOSE) {
    	printf("Solved with %d iterations. Number of extremes = %d..\n", iter+1, NNEXTREMES);
    	printf("x="); printvec(x[0], 0, NSESSIONS+NNEXTREMES+NNODES, NNODES);
	}
	for(i=0; i<NCOLS; i++) soln[i] = 0.0;
    for(i=NNODES; i<NNODES+NNEXTREMES; i++) {
        readextreme(fp, x[1], NCOLS);
        if(DEBUG && x[0][i] > 0) {
            printf("EXT%d = %g * ", i+1, x[0][i]);
            printvec(x[1], 0, NCOLS, NNODES);
            printnonzero(x[1], NNODES, NSPECTRUMS, NSESSIONS);
        }
        for(j=0; j<NCOLS; j++) {
            soln[j] += x[1][j]*x[0][i];
        }
    }

    fclose(fp);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Run time = %f\n", time_spent);
    // Print Solution
    printf("\n\n"); printf("%*s\n", 70,"-"); printf("\n\n");
    printf("Solved with %d iterations. Number of extremes = %d..\n", iter+1, NNEXTREMES);
    printf("Transmission Ranges:\n");
    printmatd(links, 0, NNODES, 0, NNODES);
    printf("Session Requests:\n");
    printmatd(requests, 0, NNODES, 0, NNODES);
    // PRINT SOLUTION
    printf("Solution:\n");
    checksolnmat(soln, NNODES, NSPECTRUMS, NSESSIONS);
    printf("\n\n"); printf("%*s\n", 70,"-"); printf("\n\n");

	if(DEBUG) {
        debug(1, 1, NNODES, NSPECTRUMS, NSESSIONS, links, src, dest, fname);
	}

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
