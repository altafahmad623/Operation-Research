/*
Name : Altaf Ahmad
Roll no: 18MA20005
Simplex method - solution
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
void simplex(double **A, double *B, double *C, int n, int m) // the function which solves the LPP using simpmlex method
{
    double sum;
    int i, j, k = n + m;
    // now we need to generate the simplex table first and foremost
    double **sim = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        sim[i] = (double *)malloc(k * sizeof(double));
    }
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            sim[i][j] = A[i][j];
        }
        for (j = n; j < (n + m); j++) //these are the coefficients of the slack variables
        {
            sim[i][j] = 0.0;
        }
        sim[i][n + i] = 1.0; // the slack variables are added for each equation and their coefficient is set to be 1 for each of the m equations
    }
    double *cb = (double *)malloc(((n + m) * sizeof(double))); // this stores the coefficients of the basic variables which are present in the optimality condition
    for (i = 0; i < n; i++)
    {
        cb[i] = C[i]; // all the normal variables have their same coefficients
    }
    for (j = n; j < (n + m); j++)
    {
        cb[j] = 0.0; // the slack variables have coefficients to be 0
    }
    int keycol, keyrow;
    double maxincmz = 0.0;                                            // this helps to find the minimum Z_i - C_i in each iteration
    double mininsol = 999999.0;                                       // this helps to find the minimum ratio in each iteration
    int *basv = (int *)malloc(m * sizeof(int));                       // stores the basic variables values
    double *Z = (double *)malloc(((n + m) * sizeof(double)));         //stores the values of Z = \sum (CB_i)*(a_i_j)
    double *CminusZ = (double *)malloc(((n + m) * sizeof(double)));   // calculates the Z_i - C_i
    double *sol = (double *)malloc(m * sizeof(double));               // stores the solutions of each iteration
    double *ratio = (double *)malloc(m * sizeof(double));             // calculates the ratios of each iteration
    double *keyrowval = (double *)malloc(((n + m) * sizeof(double))); // stores the key row in a separate array
    double *keycolval = (double *)malloc(m * sizeof(double));         // stores the key column in a separate array
    for (i = 0; i < m; i++)
    {
        basv[i] = i + n; // initially, the basic variables are put as the surplous variables
        sol[i] = B[i];   // the solution column is populated with the values of the B[i] in each eqation
    }
    int check = 1;
    int iter = 0;
    double solkey;
    double zsol;
    //now we are starting the iterations
    while (check)
    {
        if (iter >= 10) // checker for inifite iterations
            break;
        if (iter != 0) // apart from the first iteration, we don't have to change the columns
        {
            basv[keyrow] = keycol; //entering variable is put in the basic variable column
            for (i = 0; i < m; i++)
            {
                if (i == keyrow)
                {
                    sol[i] = sol[i] / keyrowval[keycol]; // for the pivot row
                }
                else
                {
                    sol[i] = sol[i] - (keycolval[i] * solkey) / keyrowval[keycol]; // for other rows
                }

                for (j = 0; j < (m + n); j++)
                {
                    if (i == keyrow)
                    {
                        sim[i][j] = sim[i][j] / keyrowval[keycol]; // for pivot row
                    }
                    else
                    {
                        sim[i][j] = sim[i][j] - (keycolval[i] * keyrowval[j]) / keyrowval[keycol]; // for other rows
                    }
                }
            }
        }
        check = 0;
        maxincmz = -90000.0;
        for (i = 0; i < (m + n); i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) // now we will calculate the Z values for each variable
            {
                sum += cb[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            CminusZ[i] = cb[i] - Z[i]; // Z_i - C_i values for each variable
            if (CminusZ[i] > 0)        // this checks the maximum condition, if we have to minimize, then the sign should be changed to Z_i - C_i <=0
            {
                check = 1;
            }
            if (CminusZ[i] > maxincmz) // finds the key column
            {
                keycol = i;
                maxincmz = CminusZ[i];
            }
        }
        sum = 0;
        for (i = 0; i < m; i++)
        {
            sum += cb[basv[i]] * sol[i]; // calculates the sum for the solution
        }
        zsol = sum;
        mininsol = 999999.0;
        for (i = 0; i < m; i++)
        {
            keycolval[i] = sim[i][keycol]; // now it evaluates the key colvalues
        }
        for (i = 0; i < m; i++)
        {
            ratio[i] = sol[i] / keycolval[i]; // finds the ratio in each iteration
            if(ratio[i] < 0)
            {
                continue;
            }
            if (ratio[i] < mininsol)
            {
                mininsol = ratio[i];
                keyrow = i;
            }
        }
        for (i = 0; i < (m + n); i++)
        {
            keyrowval[i] = sim[keyrow][i];
        }
        solkey = sol[keyrow];
        //now we print the simplex table for each iteration
        printf("\n\nIteration no: %d\n", iter);
        printf("CB Ci/basic_variables");
        for (i = 0; i < (m + n); i++)
        {
            printf(" %lf", cb[i]);
        }
        printf(" solution ratio\n");
        for (i = 0; i < m; i++)
        {
            printf("%lf %d ", cb[basv[i]], basv[i]+1);
            for (j = 0; j < (m + n); j++)
            {
                printf("%lf ", sim[i][j]);
            }
            printf("%lf %lf\n", sol[i], ratio[i]);
        }
        printf("    Z_i    ");
        for (i = 0; i < (m + n); i++)
        {
            printf("%lf ", Z[i]);
        }
        printf("\nZ-i - C_i  ");
        for (i = 0; i < (m + n); i++)
        {
            printf("%lf ", -CminusZ[i]);
        }
        printf("\nMinimum ratio is : %lf coming at pivot row : %d\n", mininsol, keyrow + 1);
        printf("Minimum Z-i - C_i is : %lf coming at pivot column: %d\n", -maxincmz, keycol + 1);
        printf("value of z is :%lf\n", zsol);
        iter++;
    }
    printf("\n The final optimal values are : ");
    for (i = 0; i < m; i++)
    {
        printf(" x_ %d = %lf ", basv[i] + 1, sol[i]);
    }
    printf(" And rest all are 0\n And the optimal value of Z is : %lf\n", zsol);
    //return sol;
}
int main() // menu driven program
{
    int n, m, i, j, k;
    printf("Enter n : ");
    scanf("%d", &n);
    printf("Enter m : ");
    scanf("%d", &m);
    double error = 0.001;
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

    //double *X = (double *)malloc(m * sizeof(double));
    simplex(A, B, C, n, m);

    return 0;
}