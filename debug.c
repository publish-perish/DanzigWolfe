#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ilcplex/cplexcheck.h>
#include <ilcplex/cplex.h>

#include "dwdecomp.h"


int debug(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest, char* dfname)
{
    FILE *fp;
    int status;
    char sname[256];
    int NNNODES = NNODES*NNODES;
    double rhs[MAXROWS];
    char sense[MAXROWS];
    int i, j, m, l, q, p;
    int NCOLS = NSPECTRUMS*NNNODES+NSESSIONS*NNNODES+NSESSIONS;
    int ROWS = 0;
    
    int *requests[NNODES];
    for(int n=0; n<NNODES; n++) requests[n] = (int*)calloc(NNODES, sizeof(int));
    for(i=0; i<NSESSIONS; i++)
        requests[src[i]][dest[i]] = 1;

    int *indices = (int*)malloc(NNNODES*NSPECTRUMS*sizeof(int)); 
    for(int i=0; i<NSPECTRUMS*NNNODES; i++) indices[i] = i;
    char *ctype = (char*)malloc(NNNODES*NSPECTRUMS*sizeof(char)); 
    for(int i=0; i<NSPECTRUMS*NNNODES; i++) ctype[i] = 'B';
    double *ub = (double*)malloc(NNNODES*NSPECTRUMS*sizeof(double)); 
    for(int i=0; i<NSPECTRUMS*NNNODES; i++) ub[i] = 1;
    char *btype = (char*)malloc(NNNODES*NSPECTRUMS*sizeof(char)); 
    for(int i=0; i<NSPECTRUMS*NNNODES; i++) btype[i] = 'U';

    // CREATE LPs
    mlp = CPXcreateprob (env, &status, "master");
    for(int n=0; n<NNODES; n++) {
        sprintf(sname, "%s%d", "subprob", n+1);
        slp[n] = CPXcreateprob (env, &status, sname);
    }

	// CREATE PROBLEM MATRIX
    double *mat[MAXROWS];
    for(i=0; i<MAXROWS; i++) mat[i] = (double*)calloc(NCOLS, sizeof(double)); 

    // Self-interference (2)
    printf("Self-interference (2)   : c%d - ", ROWS+1);
    for(m=0; m<NSPECTRUMS; m++) {
        for(i=0; i<NNODES; i++) {
            for(j=0; j<NNODES; j++) {
                if(links[j][i] == 1) {
                    mat[ROWS][NNNODES*m+NNODES*j+i] = 1;
                    for(q=0; q<NNODES; q++) {
                        if(links[i][q])
                            mat[ROWS][NNNODES*m+NNODES*i+q] = 1;
                    }
                }
            }
            rhs[ROWS] = 1;
            sense[ROWS] = 'L';
            ROWS += 1;
        }
        }
    printf("c%d\n", ROWS);
    // General interference (3)
    printf("General interference (3): c%d - ", ROWS+1);
    for(m=0; m<NSPECTRUMS; m++) {
        for(p=0; p<NNODES; p++) 
            for(q=0; q<NNODES; q++) 
                if(links[p][q] == 1) {
                    mat[ROWS][NNNODES*m+NNODES*p+q] = 1;
                    for(i=0; i<NNODES; i++) {
                        if(links[i][q]) {
                            for(j=0; j<NNODES; j++)
                                if(links[i][j]) 
                                   mat[ROWS][NNNODES*m+NNODES*i+j] = 1;
                    }}
                    rhs[ROWS] = 1;
                    sense[ROWS] = 'L';
                    ROWS += 1;
                }
    }
    printf("c%d\n", ROWSMAST);
    // Routing (4)
    printf("Routing (4)\t: c%d - ", ROWS+1);
    for(i=0; i<NNODES; i++) 
        for(j=0; j<NNODES; j++) 
            if(links[i][j]) {
                for(l=0; l<NSESSIONS; l++) {
                        //printf("setting f%d %d,%d\n",l+1, i+1, j+1);
                        mat[ROWS][NSPECTRUMS*NNNODES+NNNODES*l+NNODES*i+j] = -1;
                }
                for(m=0; m<NSPECTRUMS; m++) {
                        //printf("\t setting s%d %d,%d\n", m+1, i+1, j+1);
                        mat[ROWS][NNNODES*m+NNODES*i+j] = 1;
                }
                rhs[ROWS] = 0;
                sense[ROWS] = 'E';
                ROWS += 1;
            }
    printf("c%d\n", ROWS);
    // Routing (5)        
    printf("Routing (5)\t: c%d - ", ROWS+1);
    for(i=0; i<NNODES; i++) 
        for(l=0; l<NSESSIONS; l++) 
            if(i == dest[l]) {
                for(j=0; j<NNODES; j++) 
                    if(links[i][j])
                        mat[ROWS][NSPECTRUMS*NNNODES+NNNODES*l+NNODES*i+j] = 1;
                rhs[ROWS] = 0;
                sense[ROWS] = 'E';
                ROWS += 1;
            }
    for(l=0; l<NSESSIONS; l++) 
        for(i=0; i<NNODES; i++)
            if(i == src[l]) {
                for(j=0; j<NNODES; j++) 
                    if(links[j][i])
                        mat[ROWS][NSPECTRUMS*NNNODES+NNNODES*l+NNODES*j+i] = 1;
                rhs[ROWS] = 0;
                sense[ROWS] = 'E';
                ROWS += 1;
            }
    printf("c%d\n", ROWS);
    // Data rate (6)        
    printf("Data rate (6) \t: c%d - ", ROWS+1);
    for(l=0; l<NSESSIONS; l++) {
        //printf("\nsession %d\n", l+1);
        if(links[src[l]][dest[l]]) {
                mat[ROWS][NSPECTRUMS*NNNODES+NNNODES*l+NNODES*src[l]+dest[l]] = 1;
                //printf("\t setting src,dest = %d,%d\n", src[l]+1, dest[l]+1);
        }
        for(j=0; j<NNODES; j++) {
            if(links[src[l]][j] ){//&& (src[l] < j && j <= dest[l])) {
                mat[ROWS][NSPECTRUMS*NNNODES+NNNODES*l+NNODES*src[l]+j] = 1;
                //printf("\t setting src,j = %d,%d\n", src[l]+1, j+1);
            }
            if(links[j][dest[l]] ){//&& (src[l] >= j && j > dest[l])) {
                mat[ROWS][NSPECTRUMS*NNNODES+NNNODES*l+NNODES*j+dest[l]] = 1;
                //printf("\t setting j,dest = %d,%d\n", j+1, dest[l]+1);
            }
        }
        mat[ROWS][NSPECTRUMS*NNNODES+NSESSIONS*NNNODES+l] = -1;
        rhs[ROWS] = 0;
        sense[ROWS] = 'E';
        ROWS += 1;
        //printf("\n");
    }
    printf("c%d\n", ROWS);
    // Relay nodes (7)
    printf("Relay nodes (7)         : c%d - ", ROWS+1);
    printf("\n");
    for(l=0; l<NSESSIONS; l++) {
        //printf("\n\nsession %d\n%d --> %d\n", l+1, src[l]+1, dest[l]+1);
        for(i=0; i<NNODES; i++) {
            if(i != src[l] && i != dest[l]) { 
                for(j=0; j<NNODES; j++) {
                    if(links[i][j] ) { 
                        mat[ROWS][NSPECTRUMS*NNNODES+NNNODES*l+NNODES*i+j] = 1;
                        //printf("\t setting i,j = %d,%d\n", i+1, j+1);
                        for(p=0; p<NNODES; p++)
                            if(links[p][i]) {
                                //printf("\t setting p,i = %d,%d\n", p+1, i+1);
                                mat[ROWS][NSPECTRUMS*NNNODES+NNNODES*l+NNODES*p+i] = -1;
                                //break;
                            }
                    }
                }
            }
            rhs[ROWS] = 0;
            sense[ROWS] = 'E';
            ROWS += 1; 
        }
    }
    printf("c%d\n", ROWS);


   
    if(VERBOSE) {
        printf("\nConstraint matrix has %d rows and %d cols.\n", ROWS, NCOLS);
        printheaders(NNODES, NSPECTRUMS, NSESSIONS);
    }


    // CREATE CPLEX LP
    status = CPXchgobjsen(env, mlp, -1);
    status = CPXnewcols(env, mlp, NCOLS, NULL, NULL, NULL, NULL, NULL);
    status = CPXnewrows (env, mlp, ROWS+NNODES, NULL, NULL, NULL, NULL);                    
    

    for(i=0; i<ROWS; i++)                                                              
        for(j=0; j<NCOLS; j++) 
            status = CPXchgcoef(env, mlp, i, j, mat[i][j]);

    for(i=NCOLS-NSESSIONS; i<NCOLS; i++) {                                //!!
        status = CPXchgcoef(env, mlp, -1, i, 1.0);
    }

    for(i=0; i<ROWS; i++) {
        status = CPXchgcoef(env, mlp, i, -1, rhs[i]);
        status = CPXchgsense(env, mlp, 1, &i,  &sense[i]);
    }

    // Change problem type
    //status = CPXchgprobtype(env, mlp, 1);
    //status = CPXchgbds(env, mlp, NNNODES*NSPECTRUMS, indices, btype, ub);
    //status = CPXchgctype (env, mlp, NNNODES*NSPECTRUMS, indices, ctype);


    namecols(env, mlp, NSESSIONS, NSPECTRUMS, NNODES);
    if (DEBUG) 
        status = CPXwriteprob(env, mlp, "debug.lp", NULL);



    double *x = (double*)calloc(NCOLS, sizeof(double));
    status = CPXlpopt (env, mlp);
    //status = CPXmipopt (env, mlp);
    //status = CPXsolution(env, mlp, &mlpstat, objval, x, pi, NULL, NULL);
    status = CPXgetx(env, mlp, x, 0, NCOLS-1);
    //status = CPXsolwrite(env, mlp, "mlp.sol");

    printf("________________________________________________________________________________________________\n\n");
    printf("Problem has %d nodes, %d spectrums and %d sessions.\n\n", NNODES, NSPECTRUMS, NSESSIONS);
    printf("Transmission Ranges:\n");
    printmatd(links, 0, NNODES, 0, NNODES);
    printf("Resulting Sessions:\n");
    for(i=0; i<NSESSIONS; i++) {
        printf("f%d : ", i+1); printf("%d --(x%g)--> %d\n", src[i]+1, x[NNNODES*NSPECTRUMS+NNNODES*NSESSIONS+i], dest[i]+1);
    }
    printf("\n");
    printmatd(requests, 0, NNODES, 0, NNODES);
        
    // PRINT SOLUTION
    printf("Solution:\n");
    printvec(x, 0, NCOLS, NNODES);
    checksolnmat(x, NNODES, NSPECTRUMS, NSESSIONS);
    //printsoln(x, src, dest, NNODES, NSPECTRUMS, NSESSIONS);
    printf("________________________________________________________________________________________________\n\n\n");


    for(i=0; i<MAXROWS; i++) free(mat[i]);
    for(i=0; i<NNODES; i++) free(requests[i]);
    free(x); free(btype); free(ctype); free(indices); free(ub);

    return 0;

}



