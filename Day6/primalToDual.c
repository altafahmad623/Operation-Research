/*
Name : Altaf Ahmad
Roll no: 18MA20005
Dual Simplex method - solution
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
int main()
{
    int n, m, i, j, k;
    printf("Enter n : ");
    scanf("%d", &n);
    printf("Enter m : ");
    scanf("%d", &m);
    printf("Now, we need to insert the primal equation to dual :\n The primal is to be input as first the Z which could be maximum that will be the opposite in the dual\n");
    printf("And also, all the inequations should be in the form \n a_11 x_1 + a_12 x_2 + ... + a_1n x_n < = b_1\na_21 x_1 + a_22 x_2 + ... + a_2n x_n < = b_2 \n  ");
    printf("--------------------------------\n");
    printf(" a_m1 x_1 + a_m2 x_2 + ... + a_mn x_n < = b_m\n");
    double **A, *B;
    //printf("We have the matrix equation AX = B, where A is nxn, X is nx_answer and B is also nx_answer\n Now please enter the values of A - \n");
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
            printf("Enter A[%d][%d] :", i, j);
            scanf("%lf", &A[i][j]);
        }
    }

    printf("Now, enter the values of B\n");
    for (i = 0; i < m; i++)
    {
        printf("Enter B[%d] :", i);
        scanf("%lf", &B[i]);
    }
    printf("Now enter the coefficients of the optimality Z = c1 . x_1 + c2 . x_2 .. + cn . x_n;  so enter c1,c2...,cn\n");
    double *C = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        printf("Enter coefficient of x_%d :", i + 1);
        scanf("%lf", &C[i]);
    }
    printf("Converting it to dual we get : ");
    printf("Z = ");
    for ( i = 0; i < (m-1); i++)
    {
        printf("%0.2lf x_%d + ",B[i],i+1);
    }
    printf("%0.2lf x_%d\n",B[m-1],m );
    printf("Subject to : \n");
    for ( i = 0; i < n; i++)
    {
        for ( j = 0; j < (m-1); j++)
        {
            printf("%0.2lf x_%d + ",A[j][i],i+1);
        }
        printf("%0.2lf x_%d >= %lf \n",A[m-1][i],m ,C[i]);
        
        /* code */
    }
    
    return 0;
}