#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ilcplex/cplexcheck.h>
#include <ilcplex/cplex.h>

#include "dwdecomp.h"

int branch(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest, double *soln)
{

    int K, k, nsteps = 0;
    int NEXTREMES = 0, NROWS;
    int NNNODES = NNODES*NNODES;
    int NCOLS = NSPECTRUMS*NNNODES+NSESSIONS*NNNODES+NSESSIONS;
    double k1, k2, objval = 0;
    double K1[NCOLS], K2[NCOLS];
    int l = 0;
    char sense;

    for(i=0; i<NCOLS; i++) { K1[i] = 0.0; K2[i] = 0.0; }


// DEBUG
k1 = solvebranch(VERBOSE, DEBUG, NNODES, NSPECTRUMS, NSESSIONS, links, src, dest, soln);
return 0;


    // BRANCHING STEP
    for(n=0; n<NNODES; n++) {
        for(m=0; m<NSPECTRUMS; m++) {
            printf("\n\n\n"); printf("%*s\n", 30,"-"); printf("BRANCHING"); printf("%*s\n", 30,"-"); printf("\n\n\n");
            K = n;
            k = n*m;
            if(k == n) k++;
            nsteps += 1;
            NROWS = CPXgetnumrows(env, slp[K]);

            // Add row for branching coef (easiest way to do this)
            status = CPXnewrows(env, slp[K], 1, NULL, NULL, NULL, NULL);
            status = CPXchgcoef(env, slp[K], NROWS, NNNODES*K+K*NNODES+k, 1.0);

            if(VERBOSE) {
                printf("%*s\n", 70,"-"); printf("%*s\n", 20,"\t");
                printf("Branching on s%ds%d%d >= 1 idx %d\n", l+1, K+1, k+1, NNNODES*K+K*NNODES+k);
                printf("\nRESETTING RMP...\n\n");
            }
            // Reset RMP
            status = CPXreadcopyprob(env, rmp, "orig_rmp.lp", NULL);
            for(i=0; i<NCOLS; i++) K1[i] = 0.0; 
            // Set branching coef
            sense = 'G';
            status = CPXchgsense (env, slp[K], 1, &ROWSSUB[K], &sense);
            status = CPXchgcoef(env, slp[K], ROWSSUB[K], -1, 1.0);
            k1 = solvebranch(VERBOSE, DEBUG, NNODES, NSPECTRUMS, NSESSIONS, links, src, dest, K1);

            if(VERBOSE) {
                printf("K1 has objval %g\n", k1);
                printf("SOLUTION:\n"); checksolnmat(K1, NNODES, NSPECTRUMS, NSESSIONS);
                printf("%*s\n", 70,"-"); printf("%*s\n", 20,"\t");
                printf("Branching on s%ds%d%d <= 0\n", l+1, K+1, k+1);
                printf("\nRESETTING RMP...\n\n");
            }
            // Reset RMP
            status = CPXreadcopyprob(env, rmp, "orig_rmp.lp", NULL);
            for(i=0; i<NCOLS; i++) K2[i] = 0.0; 
            // Set branching coef
            sense = 'L';
            status = CPXchgsense (env, slp[K], 1, &NROWS, &sense);
            status = CPXchgcoef(env, slp[K], NROWS, -1, 0.0);
            status = CPXwriteprob(env, slp[K], "K2.lp", NULL);
            k2 = solvebranch(VERBOSE, DEBUG, NNODES, NSPECTRUMS, NSESSIONS, links, src, dest, K2);

            if(VERBOSE) {
                printf("K2 has objval %g\n", k2);
                printf("SOLUTION:\n"); checksolnmat(K2, NNODES, NSPECTRUMS, NSESSIONS);
            }

            if(k1 >= k2) {
                status = CPXchgcoef(env, slp[K], NROWS, -1, 1.0);
                if(k1 > objval) {// && (K1-(int)K1 == 0)) { \\ TODO: Test if soln is integer
                    if(VERBOSE) printf("Choosing K1\n");
                    objval = k1;
                    for(j=0; j<NCOLS; j++)
                        soln[j] = K1[j];
                }
            }
            else if(k2 > k1) {
                if(k2 > objval) { // && (K2-(int)K2 == 0)) {
                    if(VERBOSE) printf("Choosing K2\n");
                    objval = k2;
                    for(j=0; j<NCOLS; j++)
                        soln[j] = K2[j];
                }
            }
        }
    }
    printf("Solved with %d branching steps. Number of extremes = %d..\n", nsteps+1, NEXTREMES);

    return 0;
}
