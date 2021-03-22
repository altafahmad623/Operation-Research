/*
Name : Altaf Ahmad
Roll no: 18MA20005
bigM method - solution
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define bM 100000000 // here, the big M is defined as a sufficiently large number

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
int bigM(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin) // the function which solves the LPP using simpmlex method
{
    double sum;
    int i, j, k = n + msum;
    // now we need to generate the bigM table first and foremost
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
            sim[i][n + k] = -1.0;
            k++;
            sim[i][n + k] = 1.0;
            k++;
        }
        else
        {
            sim[i][n + k] = 1.0;
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
            cb[n + k] = 0.0;
            k++;
        }
        else if (sign[i] == 1)
        {
            cb[n + k] = 0.0;
            k++;
            if (maxmin == 0)
                cb[n + k] = -bM;
            else
            {
                cb[n + k] = bM;
            }

            k++;
        }
        else
        {
            if (maxmin == 0)
                cb[n + k] = -bM;
            else
            {
                cb[n + k] = bM;
            }
            k++;
        }
    }
    int keycol, keyrow;
    double maxincmz = 0.0;                                               // this helps to find the minimum Z_i - C_i in each iteration
    double mininsol = 999999.0;                                          // this helps to find the minimum ratio in each iteration
    int *basv = (int *)malloc(m * sizeof(int));                          // stores the basic variables values
    double *Z = (double *)malloc(((n + msum) * sizeof(double)));         //stores the values of Z = \sum (CB_i)*(a_i_j)
    double *CminusZ = (double *)malloc(((n + msum) * sizeof(double)));   // calculates the Z_i - C_i
    double *sol = (double *)malloc(m * sizeof(double));                  // stores the solutions of each iteration
    double *ratio = (double *)malloc(m * sizeof(double));                // calculates the ratios of each iteration
    double *keyrowval = (double *)malloc(((n + msum) * sizeof(double))); // stores the key row in a separate array
    double *keycolval = (double *)malloc(m * sizeof(double));            // stores the key column in a separate array
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
            basv[i] = k + n;
            k++;
        }
        else
        {
            basv[i] = k + n; // initially, the basic variables are put as the artificial variables
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
        printf(" %lf * x_%d +", cb[i], i + 1);
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
        maxincmz = -90000.0;
        for (i = 0; i < (n + msum); i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) // now we will calculate the Z values for each variable
            {
                sum += cb[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            if (maxmin == 0)
                CminusZ[i] = cb[i] - Z[i]; // Z_i - C_i values for each variable
            else
            {
                CminusZ[i] = Z[i] - cb[i];
            }

            if (CminusZ[i] > 0) // this checks the maximum condition, if we have to minimize, then the sign should be changed to Z_i - C_i <=0
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
            if ((sol[i] == 0.000) && (keycolval[i] < 0))
            {
                continue;
            }
            if (ratio[i] < 0.0000000)
            {
                continue;
            }
            if (ratio[i] <= mininsol)
            {
                mininsol = ratio[i];
                keyrow = i;
            }
        }
        for (i = 0; i < (n + msum); i++)
        {
            keyrowval[i] = sim[keyrow][i];
        }
        solkey = sol[keyrow];
        //now we print the bigM table for each iteration
        printf("\n\nIteration no: %d\n", iter);
        printf("\n\t CB_i \t C_j ");
        for (i = 0; i < (n + msum); i++)
        {
            printf("\t %lf", cb[i]);
        }
        printf("\n \t \t BV. ");
        for (i = 0; i < (msum + n); i++)
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
            for (j = 0; j < (n + msum); j++)
            {
                printf("\t %lf ", sim[i][j]);
            }
            printf("\t %lf  \n", sol[i]);
        }
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            printf("-");
        }
        printf("\n\t Z_j \t ");
        for (i = 0; i < (n + msum); i++)
        {
            printf("\t %lf", Z[i]);
        }
        printf("\t %lf", zsol);
        printf("\n \t C_j - Z-j ");
        for (i = 0; i < (n + msum); i++)
        {
            if (maxmin == 0)
                printf("\t %lf ", -CminusZ[i]);
            else
            {
                printf("\t %lf", CminusZ[i]);
            }
        }
        printf("\n");
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            printf("-");
        }
        printf("\n");
        printf("\nMinimum ratio is : %lf coming at pivot row : %d\n", mininsol, keyrow + 1);
        if (maxmin == 0)
            printf("Minimum Z-i - C_i is : %lf coming at pivot column: %d\n", -maxincmz, keycol + 1);
        else
        {
            printf("Maximum Z-i - C_i is : %lf coming at pivot column: %d\n", maxincmz, keycol + 1);
        }

        printf("value of z is :%lf\n", zsol);
        iter++;
    }
    for (i = 0; i < m; i++)
    {
        if (cb[basv[i]] < -1000000.0)
        {
            printf("The iterations have been completed and there are artificial variables in the base with values strictly greater than 0, so the problem has no solution (infeasible).\n");
            return 0;
        }
    }
    int *integer_sol = (int *)malloc(m * sizeof(int));
    double *fractional_part = (double *)malloc(m * sizeof(double));
    double maxfrac = -1;
    int maxfrac_ind;
    printf("\n The final optimal values are : ");
    for (i = 0; i < m; i++)
    {
        integer_sol[i] = (sol[i]);
        fractional_part[i] = sol[i] - integer_sol[i];
        printf(" x_ %d = %lf ", basv[i] + 1, sol[i]);
        if (maxfrac < fractional_part[i])
        {
            maxfrac = fractional_part[i];
            maxfrac_ind = i;
        }
    }
    printf(" And rest all are 0\n And the optimal value of Z is : %lf\n", zsol);
    for (i = 0; i < m; i++)
    {
        printf("Fractional part is %lf, and the integer part is %d of x_%d\n", fractional_part[i], integer_sol[i], basv[i] + 1);
    }
    printf("The maximum f_i is coming at the basic variable no:%d and it is %lf\n So, we need to add the following Gomorian Constraint :\n", maxfrac_ind, maxfrac);
    printf("%lf = ", -maxfrac);
    double *varibls = (double *)malloc((m + n) * sizeof(double));
    k = 0;
    for (i = m; i < (n + m); i++)
    {
        varibls[k] = -sim[maxfrac_ind][i];
        k++;
        printf("%lf * x_%d + ", -sim[maxfrac_ind][i], i + 1);
    }
    k = 0;
    printf(" x_%d\n", (n + msum + 1)); // now we are going to use the dual simplex method
    double **simA;
    //printf("We have the matrix equation AX = B, where A is nxn, X is nx_answer and B is also nx_answer\n Now please enter the values of A - \n");
    simA = (double **)malloc((m + 1) * sizeof(double *));
    for (i = 0; i <= m; i++)
    {
        simA[i] = (double *)malloc((m + n + 1) * sizeof(double));
    }
    k = 0;
    for (i = 0; i <= m; i++)
    {
        for (j = 0; j <= (m + n); j++)
        {
            if (i == m)
            {
                if (j < n)
                {
                    simA[i][j] = 0.0;
                }
                else if (j < (m + n))
                {
                    simA[i][j] = varibls[k];
                    k++;
                }
                else
                {
                    simA[i][j] = 1.0;
                }
            }
            else
            {
                if (j == (m + n))
                {
                    simA[i][j] = 0.0;
                }
                else
                {
                    simA[i][j] = sim[i][j];
                }
            }
        }
    }
    double *simB = (double *)malloc((m + 2) * sizeof(double));
    for (i = 0; i <= (m); i++)
    {
        if (i == (m))
        {
            simB[i] = -maxfrac;
        }
        else
        {
            simB[i] = sol[i];
        }
    }
    double *cbdash = (double *)malloc(((n + msum + 1) * sizeof(double)));
    for (i = 0; i <= (n + m); i++)
    {
        if (i == (m + n))
        {
            cbdash[i] = 0.0;
        }
        else
        {
            cbdash[i] = cb[i];
        }
    }

    printf("\n\t CB_i \t C_j ");
    for (i = 0; i <= (n + msum); i++)
    {
        printf("\t %lf", cbdash[i]);
    }
    printf("\n \t \t BV. ");
    for (i = 0; i <= (msum + n); i++)
    {
        printf("\t     x_%d", i + 1);
    }
    printf("\t Solution\n");
    for (i = 0; i < (25 * (m + n) + 20); i++)
    {
        printf("-");
    }
    printf("\n");
    for (i = 0; i <= m; i++)
    {
        if (i < m)
            printf("\t %0.2lf    x_%d ", cb[basv[i]], basv[i] + 1);
        else
            printf("\t %0.2lf    x_%d ", 0.000, (n + msum + 1));
        for (j = 0; j <= (n + msum); j++)
        {
            printf("\t %lf ", simA[i][j]);
        }
        printf("\t %lf  \n", simB[i]);
    }
    for (i = 0; i < (25 * (m + n) + 20); i++)
    {
        printf("-");
    }
    printf("\n");
    //dualSimplex(simA, simB, cbdash, (n + m + 1), (m + 1));

    return 0;
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
    ck = bigM(A, B, C, n, m, msum, sign, maxmin);
    //printf("%d",msum);
    return 0;
}