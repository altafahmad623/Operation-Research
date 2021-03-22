while (check)
    {
        if (iter >= 10) // checker for inifite iterations
            break;
        if (iter != 0) // apart from the first iteration, we don't have to change the columns
        {
            basvm[keyrow] = keycol; //entering variable is put in the basic variable column
            for (i = 0; i <= m; i++)
            {
                if (i == keyrow)
                {
                    solm[i] = solm[i] / keyrowval[keycol]; // for the pivot row
                }
                else
                {
                    solm[i] = solm[i] - (keycolval[i] * solkey) / keyrowval[keycol]; // for other rows
                }

                for (j = 0; j <= (m + n); j++)
                {
                    if (i == keyrow)
                    {
                        smim[i][j] = smim[i][j] / keyrowval[keycol]; // for pivot row
                    }
                    else
                    {
                        smim[i][j] = smim[i][j] - (keycolval[i] * keyrowval[j]) / keyrowval[keycol]; // for other rows
                    }
                }
            }
        }
        check = 0;

        // finding and print the table
        printf("\n\t Iteration : %d\n", iter);
        printf("\n\t CB_i \t C_j ");
        for (i = 0; i <= (m + n); i++)
        {
            printf("\t %lf", cbm[i]);
        }
        printf("\n \t \t BV. ");
        for (i = 0; i <= (m + n); i++)
        {
            printf("\t     x_%d", i + 1);
        }
        printf("\t solmution\n");
        for (i = 0; i < (25 * (m + n) + 20); i++)
        {
            printf("-");
        }
        printf("\n");
        for (i = 0; i < m; i++)
        {
            printf("\t %0.2lf    x_%d ", cbm[basvm[i]], basvm[i] + 1);
            for (j = 0; j < (m + n); j++)
            {
                printf("\t %lf ", smim[i][j]);
            }
            printf("\t %lf \n", solm[i]);
        }
        for (i = 0; i < (25 * (m + n) + 20); i++)
        {
            printf("-");
        }
        for (i = 0; i <= (m + n); i++)
        {
            sum = 0;
            for (j = 0; j <= m; j++) // now we will calculate the Z values for each variable
            {
                sum += cbm[basvm[j]] * smim[j][i];
            }
            Z[i] = sum;
            CminusZ[i] = cbm[i] - Z[i]; // Z_i - C_i values for each variable
        }
        printf("\n\t Z_j \t ");
        for (i = 0; i <= (m + n); i++)
        {
            printf("\t %lf", Z[i]);
        }
        sum = 0;
        for (i = 0; i <= m; i++)
        {
            sum += cbm[basvm[i]] * solm[i]; // calculates the sum for the solmution
        }
        zsolm = sum;
        printf("\t %lf", zsolm);
        printf("\n \t C_j - Z-j ");
        for (i = 0; i < (m + n); i++)
        {
            printf("\t %lf", CminusZ[i]);
        }
        printf("\n");
        for (i = 0; i < (25 * (m + n) + 20); i++)
        {
            printf("-");
        }
        printf("\n");
        // finding the leaving variable with the minimum value of solm[i]
        mininsolm = 0;
        for (i = 0; i < m; i++)
        {
            if (solm[i] < mininsolm)
            {
                mininsolm = solm[i];
                keyrow = i;
            }
            if (solm[i] < 0)
            {
                check = 1;
            }
        }
        if (check == 0)
        {
            break;
        }
        printf("The most negative value of solmution is coming at row %d, which is %lf. So the leaving variable is : x_%d \n", keyrow + 1, solm[keyrow], basvm[keyrow] + 1);
        printf("We now have to find the entering variable, so we compute the following table :\n");
        // finding the entering variable by taking ratio of leaving variable row and C[j] - Z[j]
        for (i = 0; i < (20 * (m + n) + 5); i++)
        {
            printf("-");
        }
        printf("\n");
        printf("Variables ");
        for (i = 0; i <= (m + n); i++)
        {
            printf("\t   x_%d       ", i + 1);
        }
        printf("\n");
        for (i = 0; i < (20 * (m + n) + 5); i++)
        {
            printf("-");
        }
        printf("\n -(C_j - Z-j)");
        for (i = 0; i <= (m + n); i++)
        {
            printf("\t %lf", -1 * CminusZ[i]);
        }
        printf("\n x_%d \t", keyrow + 1);
        for (i = 0; i <= (m + n); i++)
        {
            printf("\t %lf", smim[keyrow][i]);
        }
        double minratio = 100000;
        double ratio2;
        printf("\n Ratio \t");
        for (i = 0; i <= (m + n); i++)
        {
            if (smim[keyrow][i] < 0)
            {
                ratio2 = (-1 * CminusZ[i]) / smim[keyrow][i];
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
        for (i = 0; i <= (m + n); i++)
        {
            keyrowval[i] = smim[keyrow][i];
        }
        solkey = solm[keyrow];
        for (i = 0; i <= m; i++)
        {
            keycolval[i] = smim[i][keycol]; // now it evaluates the key colvalues
        }
        iter++;
    }
    if (iter < 9)
    {
        printf("\n The final optimal values are : ");
        for (i = 0; i < m; i++)
        {
            printf(" x_%d = %lf, ", basvm[i] + 1, solm[i]);
        }
        printf(" And rest all are 0\n And the optimal value of Z is : %lf\n", zsolm);
    }