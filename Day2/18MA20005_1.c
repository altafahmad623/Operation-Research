#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double *Gauss_Elim(double **A, int n) // row reduction or gauss elimination method - to transform the matrix into an upper triangular matrix  in row echelon form
{
    double c;
    double *X, total,sum;
    int i,j,k;
    X = (double *)malloc(n * sizeof(double));
    for(j=1; j<=n; j++) /* loop for the generation of upper triangular matrix*/
    {
        for(i=1; i<=n; i++)
        {
            if(i>j)
            {
                c=A[i][j]/A[j][j];
                for(k=1; k<=n+1; k++)
                {
                    A[i][k]=A[i][k]-c*A[j][k]; // does the required row operations to make it 0
                }
            }
        }
    }
    X[n-1]=A[n][n+1]/A[n][n]; // for the final variable, the other coefficients will be 0, so , the answer can directly be found in this way
    /* this loop is for backward substitution*/
    for(i=n-1; i>=1; i--)
    {
        sum=0;
        for(j=i+1; j<=n; j++)
        {
            sum=sum+A[i][j]*X[j-1];
        }
        X[i-1]=(A[i][n+1]-sum)/A[i][i]; // other variables
    }
    return X;
}
void generatefinal_answer(int permutation_n_c_m[],double **A,double *B,int initial,int values,int m,int n,double err,double* Z) // this finds the possible combinations nCm values by putting the rest of the variables(n-m) to be 0.
{
    int j=0,k=0,i;
    double *x;
    double *final_answer;
    x=(double *)malloc(m*sizeof(double));
    final_answer=(double *)malloc(n*sizeof(double));

    if (values==m)
    {
        double **A_square;
        A_square=(double ** )malloc((m+1)*sizeof(double * ));
        for (j=0;j<=m;j++)
            {
                A_square[j]=(double *)malloc((m+2)*sizeof(double)); // this adds a particular combination of values to the [m][m+1] matrix
            }
        for (j=0;j<m;j++)
        {
            for (k=0;k<m;k++)
            {

                A_square[j+1][k+1]=A[j][permutation_n_c_m[k]]; //adding the last column of the matrix to the [m][m+1] matrix
            }
        }
        for ( i = 1; i <= m; i++)
        {
            A_square[i][m+1]  = B[i-1];
        }
        

        //x = Gauss_Seidel(A_square,B,m,err);
        x  = Gauss_Elim(A_square,m);
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
            printf("x_%d  = %.3lf ",j+1,final_answer[j]);
        }
        int check = 1;
        for ( i = 0; i < n; i++) // checks the feasability of the solutions
        {
            if(final_answer[i] <0 )
            {
                check = 0;
                break;
            }
            else if(isnan(final_answer[i]))
            {
                check = 0;
                break;
            }
        }
        if(check)
        {
            double z_value = 0.0;
            printf(" Feasable solution exists: z = " );
            for ( i = 0; i < n; i++)
            {
                z_value+= final_answer[i]*Z[i];
            }
            printf("%lf ",z_value);
        }
        else
        {
            printf(" Not Feasable");
        }
        
        printf("\n");
        return;
    }

    for (i=initial; i<=n-1 && n-i >= m-values; i++)
    {
        permutation_n_c_m[values] = i;
        generatefinal_answer(permutation_n_c_m,A,B,i+1,values+1,m,n,err,Z);
    }
    //here we can also add whether we need to have maximum Z or minimum z and hence find the optimum value
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
    printf("Now enter the coefficients of the optimality Z = ax[1] + bx[2] .. so enter a,b,c...\n");
    double *Z;
    Z =  (double *)malloc(n * sizeof(double ));
    for ( i = 0; i < n; i++)
    {
        printf("Enter coefficient of x_%d :",i+1);
        scanf("%lf" , &Z[i]);
    }
    
    int * permutation_n_c_m;
    permutation_n_c_m = (int *)malloc(m*sizeof(int));
    generatefinal_answer(permutation_n_c_m,A,B,0,0,m,n,error,Z);
    
    return 0;
}