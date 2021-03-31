/*
Name: Altaf Ahmad
Roll no: 18MA20005
(a) m x n stable game solution
(b) m x n unstable game using Primal-Dual LP Method
*/
#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define bM 99999999
int maxArray(int *a, int n)
{
    int max = -bM;
    int k;
    for (int i = 0; i < n; i++)
    {
        if (a[i] > max)
        {
            max = a[i];
            k = i;
        }
    }
    return k;
}
int bigM(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin, int cnk) // the function which solves the LPP using simplex method
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
        printf(" %0.0lf * x_%d +", cb[i], i + 1);
    }
    printf("\nsubject to\n");
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < (n + msum); j++)
        {
            printf(" %0.0lf x_%d + ", sim[i][j], j + 1);
        }
        printf(" 0 = %0.0lf\n", sol[i]);
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

    printf("\n The final optimal values are : ");
    double*xx = new double[100]; 
    double* yVal = new double[100]; //  this stores the optimal strategy for player B
    double* xVal = new double[100]; //  this stores the optimal strategy for player A
    for (i = 0; i < m; i++)
    {
        printf(" x_ %d = %lf ", basv[i] + 1, sol[i]);
        xx[basv[i]] = sol[i];
        yVal[basv[i]] = (sol[i] / zsol);
        xVal[i] = (Z[i+n]/zsol);
    }
    printf(" And rest all are 0\n And the optimal value of Z is : %lf\n.Therefore, we get probabilites for B to be\n", zsol);
    for ( i = 0; i < n; i++)
    {
        cout<<"y_"<<i+1<<" = "<<yVal[i]<<"\n";
    }
    double value = (1/zsol) - cnk;
    cout<<"And the value of the game is : "<<value<<"\n\nAs for the player A, \n";
    for ( i = 0; i < m; i++)
    {
        cout<<"x_"<<i+1<<" = "<<xVal[i]<<"\n";
    }
    
    //value = (-1)* value;
    //cout<<"And the value of the game is : "<<value<<"\n";
    return 0;
}
int minArray(int *a, int n) // helper function to find the min max 
{
    int min = bM, k;
    for (int i = 0; i < n; i++)
    {
        if (a[i] < min)
        {
            min = a[i];
            k = i;
        }
    }
    return k;
}
void printTable(int **arr, int n, int m, int *rowmin, int *columnmax) // helper function to print the table 
{
    cout << "\n";
    cout << "\t\tB's Strategy\n\t";
    for (int i = 0; i < n; i++)
    {
        cout << "------------";
    }
    cout << "\n\t    |";
    for (int i = 1; i <= n; i++)
    {
        cout << "\tb_" << i;
    }
    cout << "\tRow min\n\t";
    for (int i = 0; i < n; i++)
    {
        cout << "------------";
    }
    cout << "\nA's ";
    for (int i = 0; i < m; i++)
    {
        if (i == 0)
        {
            cout << "Str";
        }
        if (i == 1)
        {
            cout << "at";
        }
        if (i == 2)
        {
            cout << "egy";
        }
        cout << "\ta_" << i + 1 << " |\t";
        for (int j = 0; j < n; j++)
        {
            cout << arr[i][j] << "\t";
        }
        cout << rowmin[i] << "\n";
    }
    cout << "\t";
    for (int i = 0; i < n; i++)
    {
        cout << "------------";
    }
    cout << "\nCol Max\t";
    for (int i = 0; i < n; i++)
    {
        cout << "\t" << columnmax[i];
    }
    cout << "\n";
}
void printTable(double **arr, int n, int m, int *rowmin, int *columnmax)
{
    cout << "\n";
    cout << "\t\tB's Strategy\n\t";
    for (int i = 0; i < n; i++)
    {
        cout << "------------";
    }
    cout << "\n\t    |";
    for (int i = 1; i <= n; i++)
    {
        cout << "\tb_" << i;
    }
    cout << "\tRow min\n\t";
    for (int i = 0; i < n; i++)
    {
        cout << "------------";
    }
    cout << "\nA's ";
    for (int i = 0; i < m; i++)
    {
        if (i == 0)
        {
            cout << "Str";
        }
        if (i == 1)
        {
            cout << "at";
        }
        if (i == 2)
        {
            cout << "egy";
        }
        cout << "\ta_" << i + 1 << " |\t";
        for (int j = 0; j < n; j++)
        {
            cout << arr[i][j] << "\t";
        }
        cout << rowmin[i] << "\n";
    }
    cout << "\t";
    for (int i = 0; i < n; i++)
    {
        cout << "------------";
    }
    cout << "\nCol Max\t";
    for (int i = 0; i < n; i++)
    {
        cout << "\t" << columnmax[i];
    }
    cout << "\n";
}
int stable_game(int **arr, int n, int m) // this function calculates the min max and max min and checks if the game is stable
{
    double **A = new double *[100];
    for (int i = 0; i < 100; i++)
    {
        A[i] = new double[100];
    }
    cout << "\nInitial Table : \n";
    int *rowmin = new int[100];
    int *columnmax = new int[100];
    int minpk = bM;
    for (int i = 0; i < m; i++)
    {
        minpk = bM;
        for (int j = 0; j < n; j++)
        {
            A[i][j] = arr[i][j];
            if (A[i][j] < minpk) // finds the minimum of all rows
            {
                minpk = A[i][j];
            }
        }
        rowmin[i] = minpk;
    }
    minpk = -bM;
    for (int j = 0; j < n; j++)
    {
        minpk = -bM;
        for (int i = 0; i < m; i++)
        {
            if (arr[i][j] > minpk) // finds the maximum of all columns
            {
                minpk = arr[i][j];
            }
        }
        columnmax[j] = minpk;
    }
    printTable(arr, n, m, rowmin, columnmax);
    int maxmin = maxArray(rowmin, m); // maxmin is the maximum element of all the minimum rows
    int minmax = minArray(columnmax, n), k1, k2; // and minmax is the minimum element of all the maximum columns
    k1 = rowmin[maxmin];
    k2 = columnmax[minmax];
    cout << "\n Max-Min = " << k1 << " , and Min-Max = " << k2 << "\n";
    if (k1 == k2) // if the game is stable
    {
        cout << "This game has a saddle point (game is stable) at (" << maxmin + 1 << " , " << minmax + 1 << ") and the value of the game = " << rowmin[maxmin] << "\n";
    }
    else // the game is not stable 
    {
        cout << "This game does not possess a saddle point (game is unstable)\nThe value of the game will lie between " << k1 << " and " << k2 << "\n";
        int c = k2 + 2;
        cout << "We are adding " << c << " to each element of the table \n"; // this is done to remove the negative elements if any
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A[i][j] += c;
            }
        }
        printTable(A, n, m, rowmin, columnmax);
        //int bigM(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin) 
        double * B = new double[100];
        for (int i = 0; i < m; i++)
        {
            B[i] = 1.0;
        }
        double * C = new double[100];
        for (int i = 0; i < n; i++)
        {
            C[i] = 1.0;
        }
        int * sign = new int[100];
        for (int i = 0; i < m; i++)
        {
            sign[i] = 0;
        }
        bigM(A,B,C,n,m,n,sign, 0, c); // now after modifying the table, we move on to the linear programming problem and solve it using simplex method
    }
}
int main()
{
    int n, m;
    cout << "Enter the number of rows(m) and columns(n), i.e, the number of strategies of A and B : \nEnter m : ";
    cin >> m;
    cout << "Enter n : ";
    cin >> n;
    int **A = new int *[100];
    for (int i = 0; i < 100; i++)
    {
        A[i] = new int[100];
    }
    cout << "\nNow, enter the values of the table : \n";
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("Enter a[%d][%d] :", i, j);
            cin >> A[i][j];
        }
    }
    int value = stable_game(A, n, m);
    return 0;
}