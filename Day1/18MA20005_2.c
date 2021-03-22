#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double *Gauss_Seidel(double **A, double *B, int n, float err)
{
    int key = 1, i, j;
    double *X, *x_answer, total;
    float allErrors;
    X = (double *)malloc(n * sizeof(double));
    x_answer = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        X[i] = 0;
    }
    key = 1;
    int counter = 0;
    while (key == 1)
    {
        counter ++;
        for (i = 0; i < n; i++)
        {
            total = B[i];
            for (j = 0; j < n; j++)
            {
                if (i != j)
                {
                    total -= A[i][j] * X[j];
                }
            }
            x_answer[i] = total / A[i][i];
            allErrors = (x_answer[i] - X[i]) / x_answer[i];
            if (allErrors < 0)
                allErrors *= (-1);
            if (allErrors > err)
            {
                key = 1;
                X[i] = x_answer[i];
            }
            else
            {
                key = 0;
            }
        }
        if(counter > 20) // failsafe if the final_answer explodes
        {
            printf("did not converge ");
            break;
        }
    }
    return x_answer;
}
void generatefinal_answer(int permutation_n_c_m[],double **A,double *B,int initial,int values,int m,int n,double err)
{
    int j=0,k=0,i;
    double *x;
    double *final_answer;
    x=(double *)malloc(m*sizeof(double));
    final_answer=(double *)malloc(n*sizeof(double));

    if (values==m)
    {
        double **A_square;
        A_square=(double ** )malloc(m*sizeof(double * ));
        for (j=0;j<m;j++)
            {
                A_square[j]=(double *)malloc(m*sizeof(double));
            }
        for (j=0;j<m;j++)
        {
            for (k=0;k<m;k++)
            {

                A_square[j][k]=A[j][permutation_n_c_m[k]];
            }
        }

        x = Gauss_Seidel(A_square,B,m,err);

        for(j=0;j<n;j++)
        {
            final_answer[j]=0;
        }
        for (k=0;k<m;k++)
        {
            final_answer[permutation_n_c_m[k]]=x[k];
        }
        for (j=0;j<n;j++)
        {
            printf("x[%d]  = %.3lf ",j,final_answer[j]);
        }
        printf("\n");
        return;
    }

    for (i=initial; i<=n-1 && n-i >= m-values; i++)
    {
        permutation_n_c_m[values] = i;
        generatefinal_answer(permutation_n_c_m,A,B,i+1,values+1,m,n,err);
    }
}
int main()
{
    int n,m, i, j, k;
    printf("Enter n : ");
    scanf("%d", &n);
    printf("Enter m : ");
    scanf("%d", &m);
    double error = 0.001;
    double **A, *B;
    printf("We have the matrix equation AX = B, where A is nxn, X is nx_answer and B is also nx_answer\n Now please enter the values of A - \n");
    A = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        A[i] = (double *)malloc(n * sizeof(double));
    }
    
    B = (double *)malloc(m * sizeof(double));
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("Enter A[%d][%d] :",i,j);
            scanf("%lf", &A[i][j]);
        }
    }
    printf("Now, enter the values of B\n");
    for (i = 0; i < m; i++)
    {
        printf("Enter B[%d] :",i);
        scanf("%lf", &B[i]);
    }
    int * permutation_n_c_m;
    permutation_n_c_m = (int *)malloc(m*sizeof(int));
    generatefinal_answer(permutation_n_c_m,A,B,0,0,m,n,error);
    
    return 0;
}