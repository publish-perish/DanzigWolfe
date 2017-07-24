#ifndef UTIL_H
#define UTIL_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define MAXPRINTED 140
extern void printsoln(double *x, int *src, int *dest, int NNODES, int NSPECTRUMS, int NSESSIONS);
extern void checksolnmat(double *x, int NNODES, int NSPECTRUMS, int NSESSIONS);
extern void print(double *vec, int N);
extern void printd(int *vec, int N);
extern void printc(char *vec, int N);
extern void printmat(double **mat, int n, int N, int m, int M);
extern void printmatd(int **mat, int n, int N, int m, int M);
extern void printvec(double *vec, int n, int N, int NNNODES);
extern void printvecd(int *vec, int n, int N);
extern void printheaders(int NNODES, int NSPECTRUMS, int NSESSIONS);
extern void shuffle(int *array, size_t array_size, size_t shuff_size);
extern void printnonzero(double *vec, int NNNODES, int NSPECTRUMS, int NSESSIONS);

#endif
