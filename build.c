#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ilcplex/cplexcheck.h>
#include <ilcplex/cplex.h>

#include "dwdecomp.h"


int build(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest)
{
    
    int NNNODES = NNODES*NNODES;
    int NCOLS = NSPECTRUMS*NNNODES+NSESSIONS*NNNODES+NSESSIONS;
    int ROWSSUB[NNODES];
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
                        if(links[i][q])
                            for(j=0; j<NNODES; j++)
                                if(links[i][j]) // j in T_i i.e. exists link i->j
                                    status = CPXchgcoef(env, mlp, ROWSMAST, NNNODES*m+NNODES*i+j, 1);
                    }
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

    return 0;
}

