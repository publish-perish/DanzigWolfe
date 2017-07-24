#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ilcplex/cplexcheck.h>
#include <ilcplex/cplex.h>

#include "dwdecomp.h"

double solvebranch(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest, double *soln)
{

    int NNEXTREMES = 0, maxidx, np;
    int NNNODES = NNODES*NNODES;
    int NCOLS = NSPECTRUMS*NNNODES+NSESSIONS*NNNODES+NSESSIONS;
    double EXTREMECOST[NNODES];
    int lpstat, found, iter = 0;
    double MAXVAL = 0, MINVAL = CPX_INFBOUND;
    double newcol[ROWSMAST+NSPECTRUMS];// = (double*)calloc(ROWSMAST+NSPECTRUMS, sizeof(double));
    double objval[NNODES+1];// = (double*)calloc(NNODES+1, sizeof(double));
    double pi[ROWSMAST+NNODES];// = (double*)calloc(ROWSMAST+NNODES, sizeof(double));
    double subcost[NNODES][NCOLS];// = (double**)calloc(NNODES, sizeof(double*));
    double x[MAXEXTREMES][NCOLS];// = (double**)calloc(MAXEXTREMES, sizeof(double*));
    for(n=0; n<NNODES; n++) EXTREMECOST[n] = 0;
    for(n=0; n<NNODES; n++) for(i=0; i<NCOLS; i++) subcost[n][i] = 0.0;// = (double*)calloc(NCOLS, sizeof(double));
    double EXTREMES[NNODES][MAXEXTREMES][NCOLS];;
    int NEXTREMES[NNODES];
    double coef;
    for(n=0; n<NNODES; n++) {
        //EXTREMES[n] = (double**)calloc(MAXEXTREMES, sizeof(double));
        for(i=0; i<MAXEXTREMES; i++) for(j=0; j<NCOLS; j++) EXTREMES[n][i][j] = 0.0; // (double*)calloc(NCOLS, sizeof(double));
        NEXTREMES[n] = 0;
    }
    for(n=0; n<MAXEXTREMES; n++) for(i=0; i<NCOLS; i++) x[n][i] =0.0;// = (double*)calloc(NCOLS, sizeof(double));
    int cnt[NNODES]; for(i=0; i<NNODES; i++) cnt[i] = 0;
    for(j=0;j<NCOLS; j++) soln[j] = 0.0;

    for(n=0; n<NNODES; n++)
        NEXTREMES[n] += 1;

    printf("\n\n"); printf("%*s\n", 20,"-"); printf("STARTING ITERATION"); printf("%*s\n", 20,"-"); printf("\n\n");
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
    for(i=0; i<MAXEXTREMES; i++) 
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
        printf("NNEXTREMES = %d\n", NNEXTREMES);
    iter += 1;
    if(NNEXTREMES >= MAXEXTREMES) {
        printf("\n Out of memory.. terminating\n");
        return 0;
    }

    // UPDATE COST AND RHS OF SUBPROB Ki
    for(n=0; n<NNODES; n++) {
        for(i=0; i<NCOLS; i++)
            subcost[n][i] = 0;
    }
    if(VERBOSE) {
        printf("UPDATING COSTS OF SUBPROBLEMS...\n");
    }

    for(n=0; n<NNODES; n++) {
        for(m=0; m<NSPECTRUMS; m++) {
            for(j=NNNODES*m+NNODES*n; j<NNNODES*m+NNODES*(n+1); j++) {
                for(i=0; i<ROWSMAST; i++) {
                        status = CPXgetcoef (env, mlp, i, j, &coef);
                        subcost[n][j] -=  pi[i] * coef;
                }
            }
        }
        for(l=0; l<NSESSIONS; l++) {
            for(j=NSPECTRUMS*NNNODES+NNNODES*l+NNODES*n; j<NSPECTRUMS*NNNODES+NNNODES*l+NNODES*(n+1); j++) {
                for(i=0; i<ROWSMAST; i++) {
                        status = CPXgetcoef (env, mlp, i, j, &coef);
                        subcost[n][j] -=  pi[i] * coef;
                }
            }
        }
    }

    // Update costs
    for(n=0; n<NNODES; n++) {
        if(pi[ROWSMAST+n]>0 ) {
             for(i=0; i<NCOLS-NSESSIONS; i++)
                    if(subcost[n][i] > 0)
                    subcost[n][i] -= pi[ROWSMAST+n];
               }
    }
    for(n=0; n<NNODES; n++) {
        if(VERBOSE) {
            printf("lp%d", n+1); printvec(subcost[n], 0, NCOLS, NNNODES);
            printnonzero(subcost[n], NNODES, NSPECTRUMS, NSESSIONS);
        }
        for(j=0; j<NCOLS; j++) {
            status = CPXchgcoef(env, slp[n], -1, j, subcost[n][j]);
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
        if(slp[n] == NULL || env == NULL) {
            printf("%p %p : CPLEX POINTER DIED. Exiting.\n", env, rmp);
            return 0;
        }
        if(VERBOSE)
            printf("\nSOLVING SUBPROBLEM %d...\n", n+1);
        if(slp[n] == NULL || env == NULL) {
            printf("%p %p : CPLEX POINTER DIED. Exiting.\n", env, slp[n]);
        }
        status = CPXoptimize (env, slp[n]);
        status = CPXsolution (env, slp[n], &lpstat, &objval[n+1], x[n+1], NULL, NULL, NULL);
        if(VERBOSE) {
            if ( lpstat != CPX_STAT_OPTIMAL ) {
				objval[n+1] = 0;
                printf("\n\t\t\t SOLUTION IS NOT OPTIMAL!\n");
            }
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
		printf("EXTREMES[%d][%d]\n", maxidx, NEXTREMES[maxidx]);
	}
	for(j=0; j<NCOLS; j++) {
		EXTREMES[maxidx][NEXTREMES[maxidx]][j] = x[maxidx+1][j];
	}
	NEXTREMES[maxidx] += 1;

    if(found) {
		if(VERBOSE)
        	printf("Number of extremes = %d, iteration %d\n", NNEXTREMES, iter+1);
    }
    else {
        printf("\nTerminating...\n");
        break;
    }
}

    // PRINT SOLUTION
	if(VERBOSE) {
    	printf("Solved with %d iterations. Number of extremes = %d..\n", iter+1, NNEXTREMES);
    	printf("x="); printvec(x[0], 0, NSESSIONS+NNEXTREMES+NNODES, NNODES);
	}
	for(i=0; i<NCOLS; i++) soln[i] = 0.0;
    for(i=0; i<NNODES+NNEXTREMES; i++) {
        for(n=0; n<NNODES; n++) {
            status = CPXgetcoef (env, rmp, ROWSMAST+n, i, &coef);
            if(coef == 1 ) {
                if(DEBUG) {
                    printf("Soln node %d = %g * ", n+1, x[0][i]);
                    printnonzero(EXTREMES[n][cnt[n]], NNODES, NSPECTRUMS, NSESSIONS);
                    printvec(EXTREMES[n][cnt[n]], 0, NCOLS, NNODES);
                    printnonzero(EXTREMES[0][cnt[0]], NNODES, NSPECTRUMS, NSESSIONS);
                }
                for(j=0; j<NCOLS; j++) {
                    soln[j] += EXTREMES[n][cnt[n]][j]*x[0][i];
                }
            cnt[n] += 1;
            }
        }
    }

	if(DEBUG) {
    	printf("\n\n"); printf("%*s\n", 70,"-"); printf("\n\n");
		printf("Transmission Ranges:\n");
		printmatd(links, 0, NNODES, 0, NNODES);
		printf("Session Requests:\n");
		for(i=0; i<NSESSIONS; i++) {
			printf("f%d : ", i+1); printf("%d --(x%g)--> %d\n", src[i]+1, x[0][NNODES+NNEXTREMES+i], dest[i]+1);
		}
		printf("\n");
		int requests[NNODES][NNODES];
		for(i=0; i<NNODES; i++) for(j=0; j<NNODES; j++) requests[i][j] = 0;// = (int*)calloc(NNODES, sizeof(int));
		for(i=0; i<NSESSIONS; i++)
			requests[src[i]][dest[i]] = 1;
		//printmatd(requests, 0, NNODES, 0, NNODES);
		printf("SOLUTION:\n");
		checksolnmat(soln, NNODES, NSPECTRUMS, NSESSIONS);
		printnonzero(soln, NNODES, NSPECTRUMS, NSESSIONS);
    	printf("\n\n"); printf("%*s\n", 70,"-"); printf("\n\n");
	}

    
    double retval = objval[0];
    //for(i=0; i<NNODES; i++) for(n=0; n<MAXEXTREMES; n++) free(EXTREMES[i][n]); 
    //for(i=0; i<NNODES; i++) { free(subcost[i]); free(EXTREMES[i]);}
    //for(n=0; n<MAXEXTREMES; n++) {free(x[n]);}
    //TODO: WHY CAN I NOT FREE THIS MEMORY?
    /*
    free(pi); 
    free(x); 
    free(newcol); 
    free(subcost);
    free(objval); 
    */

    return retval;

}


