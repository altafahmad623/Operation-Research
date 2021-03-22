#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double determinant(double[][25], double);
void cofactor(double[][25], double);
void transpose(double[][25], double[][25], double);
int main()
{
    double a[25][25], k, d;
    int i, j;
    printf("Enter the order of the Matrix : ");
    scanf("%lf", &k);
    printf("Enter the elements of %.0fX%.0f Matrix : \n", k, k);
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < k; j++)
        {
            scanf("%lf", &a[i][j]);
        }
    }
    d = determinant(a, k);
    if (d == 0)
        printf("\nInverse of Entered Matrix is not possible\n");
    else
        cofactor(a, k);
}
/*For calculating Determinant of the Matrix */
double determinant(double a[25][25], double k)
{
    double s = 1, det = 0, b[25][25];
    int i, j, m, n, c;
    if (k == 1)
    {
        return (a[0][0]);
    }
    else
    {
        det = 0;
        for (c = 0; c < k; c++)
        {
            m = 0;
            n = 0;
            for (i = 0; i < k; i++)
            {
                for (j = 0; j < k; j++)
                {
                    b[i][j] = 0;
                    if (i != 0 && j != c)
                    {
                        b[m][n] = a[i][j];
                        if (n < (k - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (a[0][c] * determinant(b, k - 1));
            s = -1 * s;
        }
    }
    return (det);
}
void cofactor(double num[25][25], double f)
{
    double b[25][25], fac[25][25];
    int p, q, m, n, i, j;
    for (q = 0; q < f; q++)
    {
        for (p = 0; p < f; p++)
        {
            m = 0;
            n = 0;
            for (i = 0; i < f; i++)
            {
                for (j = 0; j < f; j++)
                {
                    if (i != q && j != p)
                    {
                        b[m][n] = num[i][j];
                        if (n < (f - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[q][p] = pow(-1.0, (q + p)) * determinant(b, f - 1);
        }
    }
    transpose(num, fac, f);
}

/*Finding transpose of matrix*/

void transpose(double num[25][25], double fac[25][25], double r)
{
    int i, j;
    double b[25][25], inverse[25][25], d;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            b[i][j] = fac[j][i];
        }
    }
    d = determinant(num, r);
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            inverse[i][j] = b[i][j] / d;
        }
    }
    printf("\n\n\nThe inverse of matrix is : \n");
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            printf("\t%lf", inverse[i][j]);
        }

        printf("\n");
    }
}