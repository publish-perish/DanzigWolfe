#ifndef DWDECOMP_H
#define DWDECOMP_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <ilcplex/cplexcheck.h>
#include <ilcplex/cplex.h>

#include "util.h"

#define MAXROWS 10000
#define MAXSUBS 100
#define MAXSESSIONS 50
#define MAXSPECTRUMS 50
#define MAXNCOLS 10000
#define MAXEXTREMES 10000
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

extern CPXENVptr     env;
extern CPXLPptr      slp[MAXSUBS];
extern CPXLPptr      mlp;
extern CPXLPptr      rmp;

extern int status;
extern int flag;
extern char fname[256], sname[256];
extern int i, j, m, n, l, q, p, idx;
extern int ROWSMAST;
extern int ROWSSUB[MAXSUBS];
extern int LIVE[MAXSUBS];

extern int build(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest);

extern int branch(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest, double *soln);

extern double solvebranch(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest, double *soln);

extern int debug(int VERBOSE, int DEBUG, int NNODES, int NSPECTRUMS, int NSESSIONS, int **links, int *src, int *dest, char* dfname);


extern void namecols(CPXENVptr env, CPXLPptr lp, int NSESSIONS, int NSPECTRUMS, int NNODES);




#endif
