/*
Name : Altaf Ahmad
Roll no: 18MA20005
Dual Simplex method - solution
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
void dualSimplex(double **A, double *B, double *C, int n, int m)
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
        basv[i] = i + n; // initially, the basic variables are put as the slack variables
        sol[i] = B[i];   // the solution column is populated with the values of the B[i] in each eqation
    }
    int check = 1;
    int iter = 0;
    double solkey;
    double zsol;
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

        // finding and print the table
        printf("\n\t Iteration : %d\n", iter);
        printf("\n\t CB_i \t C_j ");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", cb[i]);
        }
        printf("\n \t \t BV. ");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t     x_%d", i + 1);
        }
        printf("\t Solution\n");
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            printf("-");
        }
        printf("\n");
        for (i = 0; i < m; i++)
        {
            printf("\t %0.2lf    x_%d ", cb[basv[i]], basv[i] + 1);
            for (j = 0; j < (m + n); j++)
            {
                printf("\t %lf ", sim[i][j]);
            }
            printf("\t %lf \n", sol[i]);
        }
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            printf("-");
        }
        for (i = 0; i < (m + n); i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) // now we will calculate the Z values for each variable
            {
                sum += cb[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            CminusZ[i] = cb[i] - Z[i]; // Z_i - C_i values for each variable
        }
        printf("\n\t Z_j \t ");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", Z[i]);
        }
        sum = 0;
        for (i = 0; i < m; i++)
        {
            sum += cb[basv[i]] * sol[i]; // calculates the sum for the solution
        }
        zsol = sum;
        printf("\t %lf", zsol);
        printf("\n \t C_j - Z-j ");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", CminusZ[i]);
        }
        printf("\n");
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            printf("-");
        }
        printf("\n");
        // finding the leaving variable with the minimum value of sol[i]
        mininsol = 0;
        for (i = 0; i < m; i++)
        {
            if (sol[i] < mininsol)
            {
                mininsol = sol[i];
                keyrow = i;
            }
            if (sol[i] < 0)
            {
                check = 1;
            }
        }
        if (check == 0)
        {
            break;
        }
        printf("The most negative value of Solution is coming at row %d, which is %lf. So the leaving variable is : x_%d \n", keyrow + 1, sol[keyrow], basv[keyrow] + 1);
        printf("We now have to find the entering variable, so we compute the following table :\n");
        // finding the entering variable by taking ratio of leaving variable row and C[j] - Z[j]
        for (i = 0; i < (20 * (m + n) + 5); i++)
        {
            printf("-");
        }
        printf("\n");
        printf("Variables ");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t   x_%d       ", i + 1);
        }
        printf("\n");
        for (i = 0; i < (20 * (m + n) + 5); i++)
        {
            printf("-");
        }
        printf("\n -(C_j - Z-j)");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", -1 * CminusZ[i]);
        }
        printf("\n x_%d \t", keyrow + 1);
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", sim[keyrow][i]);
        }
        double minratio = 100000;
        double ratio2;
        printf("\n Ratio \t");
        for (i = 0; i < (m + n); i++)
        {
            if (sim[keyrow][i] < 0)
            {
                ratio2 = (-1 * CminusZ[i]) / sim[keyrow][i];
                printf("\t %lf", ratio2);
                if (ratio2 < minratio)
                {
                    minratio = ratio2;
                    keycol = i;
                }
            }
            else
            {
                printf("\t  --      ");
            }
        }
        printf("\n");
        for (i = 0; i < (20 * (m + n) + 5); i++)
        {
            printf("-");
        }
        printf("\n Here, the minimum value of the Ratio is %lf. So the entering variable is x_%d \n", minratio, keycol + 1);
        for (i = 0; i < (m + n); i++)
        {
            keyrowval[i] = sim[keyrow][i];
        }
        solkey = sol[keyrow];
        for (i = 0; i < m; i++)
        {
            keycolval[i] = sim[i][keycol]; // now it evaluates the key colvalues
        }
        iter++;
    }
    if (iter < 9)
    {
        printf("\n The final optimal values are : ");
        for (i = 0; i < m; i++)
        {
            printf(" x_%d = %lf, ", basv[i] + 1, sol[i]);
        }
        printf(" And rest all are 0\n And the optimal value of Z is : %lf\n", zsol);
    }
}
int main() // menu driven program
{
    int n, m, i, j, k;
    printf("Enter n : ");
    scanf("%d", &n);
    printf("Enter m : ");
    scanf("%d", &m);
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
    dualSimplex(A, B, C, n, m);

    return 0;
}