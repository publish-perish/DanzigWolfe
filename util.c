#include "util.h"

void printsoln(double *x, int *src, int *dest, int NNODES, int NSPECTRUMS, int NSESSIONS)
{
    int NNNODES = NNODES*NNODES;
    printf("\n");
    for(int l=0; l<NSESSIONS; l++) {
        printf("l%d : ", l+1); printf("%d ---> %d\t||", src[l]+1, dest[l]+1);
        if(src[l] < dest[l]) {
                for(int i=0; i<NNODES; i++)
                    for(int j=0; j<NNODES; j++) {
                        if(x[NSPECTRUMS*NNNODES+l*NNNODES+NNODES*i+j] > 0) 
                            printf(" %d ----(x%d)----> %d ", i+1, (int)x[NSPECTRUMS*NNNODES+l*NNNODES+NNODES*i+j], j+1); 
                    }
        }
        else if(src[l] > dest[l]) {
                for(int i=0; i<NNODES; i++)
                    for(int j=0; j<NNODES; j++) {
                        if(x[NSPECTRUMS*NNNODES+l*NNNODES+NNODES*i+j] > 0) 
                            printf(" %d ----(x%d)----> %d ", i+1, (int)x[NSPECTRUMS*NNNODES+l*NNNODES+NNODES*i+j], j+1); 
                    }
        }
        printf("\n\n");
    }
    printf("\n\n");
}

void print(double *vec, int N){
    for(int i=0; i<N; i++) printf("%g ", vec[i]);
    printf("\n");
}
void printd(int *vec, int N){
    for(int i=0; i<N; i++) printf("%d ", vec[i]);
    printf("\n");
}
void printc(char *vec, int N){
    for(int i=0; i<N; i++) printf("%c ", vec[i]);
    printf("\n");
}
void printmat(double **mat, int n, int N, int m, int M){
    for(int i=n; i<N; i++) {
        //if(i-n > MAXPRINTED) { printf(".\n.\n.\n"); break; }
        if(i<9) printf("c%d  [", i+1); 
        else if(i<99) printf("c%d [", i+1); 
        else printf("c%d[", i+1); 
        for(int j=m; j<M; j++) {
            if(j-m > MAXPRINTED) { printf("..."); break; }
            if(mat[i][j] >= 0) printf(" %g ", mat[i][j]);
            else printf("%g ", mat[i][j]);
        }
        printf("]\n"); 
    }
    printf("\n"); 
}
void printmatd(int **mat, int n, int N, int m, int M){
    for(int i=n; i<N; i++) {
        printf("%d [", i+1); 
        for(int j=m; j<M; j++) {
            printf("%d ", mat[i][j]);
        }
        printf("]\n"); 
    }
    printf("\n"); 
}
void printvec(double *vec, int n, int N, int NNNODES)
{
    printf("[");
    for(int i=n; i<N; i++) {
        if(i-n > MAXPRINTED) {
            printf("\n");
            return;
        }
        else if(vec[i] >= 0) printf(" %g ", fabs(vec[i]));
        else printf("%g ", vec[i]);
    }
    printf("]\n");
}
void printvecd(int *vec, int n, int N)
{
    printf("[");
    for(int i=n; i<N; i++) {
        if(i-n > MAXPRINTED) {
            printf("\n");
            return;
        }
        else if(vec[i] >= 0) printf(" %d ", abs(vec[i]));
        else printf("%d ", vec[i]);
    }
    printf("]\n");
}
void printheaders(int NNODES, int NSPECTRUMS, int NSESSIONS)
{
    int cnt = 0;
    for(int m=1; m<NSPECTRUMS+1; m++) {
        for(int i=1; i<NNODES+1; i++)
            for(int j=1; j<NNODES+1; j++){
                printf("s%d%d(%d) ", i, j, m*i*j); cnt++;
                if(cnt > MAXPRINTED) {
                    printf("\n");
                    return;
                }
            }
        printf("|");
    }
    for(int l=1; l<NSESSIONS+1; l++) {
        for(int i=1; i<NNODES+1; i++)
            for(int j=1; j<NNODES+1; j++) {
                printf("f%d%d ", i, j); cnt++;
                if(cnt > MAXPRINTED) {
                    printf("\n");
                    return;
                }
            }
        printf("|");
    }
    for(int l=1; l<NSESSIONS+1; l++)
        printf(" r%d", l);
    printf("\n");
}
void shuffle(int *array, size_t array_size, size_t shuff_size)
{   
	if (array_size > 1)  
	{   
		size_t i;
		for (i = 0; i < shuff_size - 1; i++) 
		{   
		  size_t j = i + rand() / (RAND_MAX / (array_size - i) + 1); 
		  int t = array[j];
		  array[j] = array[i];
		  array[i] = t;
		}   
	}   
}  
void printnonzero(double *vec, int NNODES, int NSPECTRUMS, int NSESSIONS) {
    printf("\n     ");
    for(int m=0; m<NSPECTRUMS; m++) {
        for(int i=0; i<NNODES; i++)
            for(int j=0; j<NNODES; j++){
                if(vec[NNODES*NNODES*m+NNODES*i+j] > 0)
                    printf("s%d_%d-%d(x%g) ", m+1, i+1, j+1, vec[NNODES*NNODES*m+NNODES*i+j]); 
            }
    }
    for(int l=0; l<NSESSIONS; l++) {
        for(int i=0; i<NNODES; i++)
            for(int j=0; j<NNODES; j++) {
                if(vec[NNODES*NNODES*NSPECTRUMS+NNODES*NNODES*l+NNODES*i+j] > 0)
                    printf(" f%d_%d-%d(x%g) ", l+1, i+1, j+1, vec[NNODES*NNODES*NSPECTRUMS+NNODES*NNODES*l+NNODES*i+j]); 
            }
    }
    printf("\n\n");
}

void checksolnmat(double *x, int NNODES, int NSPECTRUMS, int NSESSIONS)
{
    double *mat[NNODES];
    for(int n=0; n<NNODES; n++) mat[n] = (double*)calloc(NNODES, sizeof(double));
    int NNNODES = NNODES*NNODES;
    printf("\n");
    for(int l=0; l<NSESSIONS; l++) {
                for(int i=0; i<NNODES; i++)
                    for(int j=0; j<NNODES; j++) {
                        if(x[NSPECTRUMS*NNNODES+l*NNNODES+NNODES*i+j] > 0) 
                            mat[i][j] += x[NSPECTRUMS*NNNODES+l*NNNODES+NNODES*i+j];
                    }
    }
    printmat(mat, 0, NNODES, 0, NNODES);
    printf("\n");

    for(int n=0; n<NNODES; n++) free(mat[n]);
}
void writeextreme(FILE *fp, double *vec, int NCOLS)
{
    for(int i=0; i<NCOLS; i++)
        fprintf(fp, "%lf ", vec[i]);
    fprintf(fp, "\n");
}
void readextreme(FILE *fp, double *vec, int NCOLS)
{
    for(int i=0; i<NCOLS; i++)
        fscanf(fp, "%lf", &vec[i]);

}
