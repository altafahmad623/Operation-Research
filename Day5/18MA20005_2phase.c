/*
Name : Altaf Ahmad
Roll no: 18MA20005
2 - Phase simplex method - solution
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int twoPhase(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin) // the function which solves the LPP using simpmlex method
{
    double sum;
    int i, j, k = n + msum;
    // now we need to generate the twoPhase table first and foremost
    double **sim = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        sim[i] = (double *)malloc(k * sizeof(double));
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            sim[i][j] = A[i][j];
        }
        for (j = n; j < (n + msum); j++) //these are the coefficients of the slack variables
        {
            sim[i][j] = 0.0;
        }
        //now here, we need to check the sign[i]
        //sim[i][n + i] = 1.0; // the slack variables are added for each equation and their coefficient is set to be 1 for each of the m equations
        if (sign[i] == 1)
        {
            sim[i][n + k] = -1.0; // these are the coefficients of the surplus variables
            k++;
            sim[i][n + k] = 1.0; // these are the coefficients of the artificial variables
            k++;
        }
        else
        {
            sim[i][n + k] = 1.0; // these are the coefficients of the artificial variables
            k++;
        }
    }
    double *cb = (double *)malloc(((n + msum) * sizeof(double))); // this stores the coefficients of the basic variables which are present in the optimality condition
    for (i = 0; i < n; i++)
    {
        cb[i] = C[i]; // all the normal variables have their same coefficients
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            cb[n + k] = 0.0; // these are the coefficients of the slack variables in the optimal function Z
            k++;
        }
        else if (sign[i] == 1)
        {
            cb[n + k] = 0.0; // these are the coefficients of the surplus variables in the optimal function Z
            k++;
            cb[n + k] = 1; // these are the coefficients of the artificial variables in the optimal function Z
            k++;
        }
        else
        {
            cb[n + k] = 1; // these are the coefficients of the artificial variables in the optimal function Z
            k++;
        }
    }
    if (maxmin == 0)
    {
        for (i = 0; i < (n + msum); i++)
        {
            cb[i] = -cb[i];
        }
    }
    // Phase 1 starts now :
    //we have to define cbdash for the new revised optimality condition Z = 0 x_1 + 0x_2 + .. + a_1 +a_2.. (artificial variables are added only)
    k = 0;
    double *cbdash = (double *)malloc((n + msum) * sizeof(double));
    for (i = 0; i < n; i++)
    {
        cbdash[i] = 0; // all the normal variables have their same coefficients
    }
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            cbdash[n + k] = 0.0; // these are the coefficients of the slack variables in the optimal function Z
            k++;
        }
        else if (sign[i] == 1)
        {
            cbdash[n + k] = 0.0; // these are the coefficients of the surplus variables in the optimal function Z
            k++;
            cbdash[n + k] = 1; // these are the coefficients of the artificial variables in the optimal function Z
            k++;
        }
        else
        {
            cbdash[n + k] = 1; // these are the coefficients of the artificial variables in the optimal function Z
            k++;
        }
    }
    if (maxmin == 0)
    {
        for (i = 0; i < (n + msum); i++)
        {
            cbdash[i] = -cbdash[i];
        }
    }
    int keycol, keyrow;
    double minimcmz = 0.0;                                               // this helps to find the minimum Z_i - C_i in each iteration
    double maxsol = 999999.0;                                            // this helps to find the minimum ratio in each iteration
    int *basv = (int *)malloc(m * sizeof(int));                          // stores the basic variables values
    double *Z = (double *)malloc(((n + msum) * sizeof(double)));         //stores the values of Z = \sum (CB_i)*(a_i_j)
    double *CminusZ = (double *)malloc(((n + msum) * sizeof(double)));   // calculates the Z_i - C_i
    double *sol = (double *)malloc(m * sizeof(double));                  // stores the solutions of each iteration
    double *ratio = (double *)malloc(m * sizeof(double));                // calculates the ratios of each iteration
    double *keyrowval = (double *)malloc(((n + msum) * sizeof(double))); // stores the key row in a separate array
    double *keycolval = (double *)malloc(m * sizeof(double));            // stores the key column in a separate array
    int *levvar = (int *)malloc((n + msum) * sizeof(int));

    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            basv[i] = k + n; // initially, the basic variables are put as the surplous variables
            k++;
        }
        else if (sign[i] == 1)
        {
            k++;
            basv[i] = k + n; // initially, the basic variables are put as the artificial variables
            levvar[k+n] = 1;
            k++;
        }
        else
        {
            basv[i] = k + n; // initially, the basic variables are put as the artificial variables
            levvar[k+n] = 1;
            k++;
        }
        sol[i] = B[i]; // the solution column is populated with the values of the B[i] in each eqation
    }
    int check = 1;
    int iter = 0;
    double solkey;
    double zsol;
    //now we are starting the iterations
    printf("Now, after modifying the equations, we get : \n");
    if (maxmin)
    {
        printf("Minimize: Z = ");
    }
    else
    {
        printf("Maximize: Z = ");
    }
    for (i = 0; i < (n + msum); i++)
    {
        printf(" %lf * x_%d +", cbdash[i], i + 1);
        levvar[i] = 0;
    }
    printf("\nsubject to\n");
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < (n + msum); j++)
        {
            printf(" %lf x_%d + ", sim[i][j], j + 1);
        }
        printf(" 0 = %lf\n", sol[i]);
    }
    printf("and \n");
    for (i = 0; i < (n + msum - 1); i++)
    {
        printf(" x_%d,", i + 1);
    }
    printf(" x_%d >= 0\n", i + 1);
    //entering loop for phase 1

    while (check)
    {
        if (iter >= 10) // checker for inifite iterations
        {
            printf("There are no feasable solutions\n");
            break;
        }
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

                for (j = 0; j < (n + msum); j++)
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
        minimcmz = 0;
        k=0;
        for (i = 0; i < (n + msum); i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) // now we will calculate the Z values for each variable
            {
                sum += cbdash[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            CminusZ[i] = cbdash[i] - Z[i]; // Z_i - C_i values for each variable
            if (CminusZ[i] < 0)            // this checks the maximum condition, if we have to minimize, then the sign should be changed to Z_i - C_i <=0
            {
                check = 1;
            }
            if (CminusZ[i] < minimcmz) // finds the key column
            {
                keycol = i;
                minimcmz = CminusZ[i];
            }
        }
        sum = 0;
        for (i = 0; i < m; i++)
        {
            sum += cbdash[basv[i]] * sol[i]; // calculates the sum for the solution
        }
        zsol = sum;
        maxsol = 0;
        for (i = 0; i < m; i++)
        {
            keycolval[i] = sim[i][keycol]; // now it evaluates the key colvalues
        }
        for (i = 0; i < m; i++)
        {
            ratio[i] = sol[i] / keycolval[i]; // finds the ratio in each iteration
            if ((sol[i] == 0.000) && (keycolval[i] < 0))
            {
                continue;
            }
            if (ratio[i] < 0.0000000)
            {
                continue;
            }
            if (sol[i] >= maxsol)
            {
                maxsol = sol[i];
                keyrow = i;
            }
        }
        for (i = 0; i < (n + msum); i++)
        {
            keyrowval[i] = sim[keyrow][i];
        }
        solkey = sol[keyrow];
        //now we print the twoPhase table for each iteration
        printf("\n\nIteration no: %d\n", iter);
        printf("CB  Ci/basic_variables  ");
        for (i = 0; i < (n + msum); i++)
        {
            printf("  %lf", cbdash[i]);
        }
        printf("  solution ratio\n");
        for (i = 0; i < m; i++)
        {
            printf("%lf  %d  ", cbdash[basv[i]], basv[i] + 1);
            for (j = 0; j < (n + msum); j++)
            {
                printf("%lf  ", sim[i][j]);
            }
            printf("%lf  %lf\n", sol[i], ratio[i]);
        }
        printf("  Z_i  ");
        for (i = 0; i < (n + msum); i++)
        {
            printf("%lf  ", Z[i]);
        }
        printf("\nZ-i - C_i  ");
        for (i = 0; i < (n + msum); i++)
        {
            if (maxmin == 0)
                printf("%lf  ", -CminusZ[i]);
            else
            {
                printf("%lf  ", CminusZ[i]);
            }
        }
        printf("\nMinimum ratio is : %lf coming at pivot row : %d\n", maxsol, keyrow + 1);
        if (maxmin == 0)
            printf("Minimum Z-i - C_i is : %lf coming at pivot column: %d\n", -minimcmz, keycol + 1);
        else
        {
            printf("Minimum C_i - Z-i is : %lf coming at pivot column: %d\n", minimcmz, keycol + 1);
        }

        printf("value of z is :%lf\n", zsol);
        iter++;
    }
    if (zsol == 0)
    {
        check = 1;
        printf("\nNow entering Phase 2\n");
        iter = 0;
        while (check)
        {
            if (iter >= 20) // checker for inifite iterations
            {
                printf("There are no feasable solutions\n");
                break;
            }
            if (iter != 0) // apart from the first iteration, we don't have to change the columns
            {
                //levvar[basv[keyrow]] = 1;
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

                    for (j = 0; j < (n + msum); j++)
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
            minimcmz = 0;
            for (i = 0; i < (n + msum); i++)
            {
                sum = 0;
                for (j = 0; j < m; j++) // now we will calculate the Z values for each variable
                {
                    sum += cb[basv[j]] * sim[j][i];
                }
                Z[i] = sum;
                CminusZ[i] = cb[i] - Z[i]; // Z_i - C_i values for each variable
                if (CminusZ[i] < 0)        // this checks the maximum condition, if we have to minimize, then the sign should be changed to Z_i - C_i <=0
                {
                    check = 1;
                }
                if (CminusZ[i] < minimcmz) // finds the key column
                {
                    keycol = i;
                    minimcmz = CminusZ[i];
                }
            }
            sum = 0;
            for (i = 0; i < m; i++)
            {
                sum += cb[basv[i]] * sol[i]; // calculates the sum for the solution
            }
            zsol = sum;
            maxsol = 0;
            for (i = 0; i < m; i++)
            {
                keycolval[i] = sim[i][keycol]; // now it evaluates the key colvalues
            }
            for (i = 0; i < m; i++)
            {
                ratio[i] = sol[i] / keycolval[i]; // finds the ratio in each iteration
                if ((sol[i] == 0.000) && (keycolval[i] < 0))
                {
                    continue;
                }
                if (ratio[i] < 0.0000000)
                {
                    continue;
                }
                if (sol[i] >= maxsol)
                {
                    maxsol = sol[i];
                    keyrow = i;
                }
            }
            for (i = 0; i < (n + msum); i++)
            {
                keyrowval[i] = sim[keyrow][i];
            }
            solkey = sol[keyrow];
            //now we print the twoPhase table for each iteration
            printf("\n\nIteration no: %d\n", iter);
            printf("CB  Ci/basic_variables  ");
            for (i = 0; i < (n + msum); i++)
            {
                printf("  %lf", cb[i]);
            }
            printf("  solution ratio\n");
            for (i = 0; i < m; i++)
            {
                printf("%lf  %d  ", cb[basv[i]], basv[i] + 1);
                for (j = 0; j < (n + msum); j++)
                {
                    printf("%lf  ", sim[i][j]);
                }
                printf("%lf  %lf\n", sol[i], ratio[i]);
            }
            printf("  Z_i  ");
            for (i = 0; i < (n + msum); i++)
            {
                printf("%lf  ", Z[i]);
            }
            printf("\nZ-i - C_i  ");
            for (i = 0; i < (n + msum); i++)
            {
                if (maxmin == 0)
                    printf("%lf  ", -CminusZ[i]);
                else
                {
                    printf("%lf  ", CminusZ[i]);
                }
            }
            printf("\nMinimum ratio is : %lf coming at pivot row : %d\n", maxsol, keyrow + 1);
            if (maxmin == 0)
                printf("Minimum Z-i - C_i is : %lf coming at pivot column: %d\n", -minimcmz, keycol + 1);
            else
            {
                printf("Maximum Z-i - C_i is : %lf coming at pivot column: %d\n", minimcmz, keycol + 1);
            }

            printf("value of z is :%lf\n", zsol);
            iter++;
        }
        printf("\n The final optimal values are : ");
        for (i = 0; i < m; i++)
        {
            printf(" x_ %d = %lf ", basv[i] + 1, sol[i]);
        }
        printf(" And rest all are 0\n And the optimal value of Z is : %lf\n", zsol);

        return 0;
    }
    else
    {
        printf("It does not have a feasable solution\n");
        return 0;
    }
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
    int *sign = (int *)malloc(m * sizeof(int));

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
    int msum = 0;
    printf("Now, please enter 0 for <= ; 1 for >= and 2 for = in each of the m= %d equations \n", m);
    for (i = 0; i < m; i++)
    {
        scanf("%d", &sign[i]);
        if (sign[i] == 1)
        {
            msum += 2;
        }
        else
        {
            msum++;
        }
    }
    printf("Now enter the coefficients of the optimality Z = c1 . x_1 + c2 . x_2 .. + cn . x_n;  so enter c1,c2...,cn\n");
    double *C = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        printf("Enter coefficient of x_%d :", i + 1);
        scanf("%lf", &C[i]);
    }
    printf("Enter 0 for Maximization and 1 for Minimization : ");
    int maxmin = 0, ck;
    scanf("%d", &maxmin);
    //double *X = (double *)malloc(m * sizeof(double));
    ck = twoPhase(A, B, C, n, m, msum, sign, maxmin);
    //printf("%d",msum);
    return 0;
}