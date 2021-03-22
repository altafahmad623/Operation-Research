#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double *Gauss_Seidel(double **A, double *B, int n, float err)
{
    int key = 1, i, j;
    double *X, *X1, total;
    float allErrors;
    X = (double *)malloc(n * sizeof(double));
    X1 = (double *)malloc(n * sizeof(double));
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
            X1[i] = total / A[i][i];
            allErrors = (X1[i] - X[i]) / X1[i];
            if (allErrors < 0)
                allErrors *= (-1);
            if (allErrors > err)
            {
                key = 1;
                X[i] = X1[i];
            }
            else
            {
                key = 0;
            }
        }
        if(counter > 10) // failsafe if the solution explodes
        {
            break;
        }
    }
    return X1;
}
int main()
{
    int n, i, j, k;
    printf("Enter n :");
    scanf("%d", &n);
    double error = 0.001;
    double **A, *B;
    printf("We have the matrix equation AX = B, where A is nxn, X is nx1 and B is also nx1\n Now please enter the values of A - \n");
    A = (double **)malloc(n * sizeof(double *));
    for (i = 0; i < n; i++)
    {
        A[i] = (double *)malloc(n * sizeof(double));
    }
    
    B = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("Enter A[%d][%d] :",i,j);
            scanf("%lf", &A[i][j]);
        }
    }
    printf("Now, enter the values of B\n");
    for (i = 0; i < n; i++)
    {
        printf("Enter B[%d] :",i);
        scanf("%lf", &B[i]);
    }
    double *x = (double *)malloc(n * sizeof(double));
    x = Gauss_Seidel(A,B,n,error);
    printf("The solution is : \n");
    for ( i = 0; i < n; i++)
    {
        printf("x[%d] = %lf ",i,x[i]);
    }
    
    return 0;
}