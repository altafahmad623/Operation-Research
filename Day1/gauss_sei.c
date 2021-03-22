#include <stdio.h>
#include <stdlib.h>


int isDiagonallyDominant(double **A,int n)
{
    int i,j;
    int flag=1;
    double sum=0,temp,temp2;
    for (i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {

            if(j!=i){

                if(A[i][j]<0) temp = (-1)*A[i][j];
                else temp = A[i][j];

                sum += temp;
            }
        }

        if(A[i][i]<0) temp2 = (-1)*A[i][i];
        else temp2 = A[i][i];

        if(temp2 < sum)
        {
            flag=0;
            return flag;
        }
    }
    return flag;
}


double * gaussSeidel(double **A,double *B,int n,float err)
{
    int key=1,i,j;
    double *X,*X1,sum;
    float calErr;

    X = (double *)malloc(n*sizeof(double));
    X1 = (double *)malloc(n*sizeof(double));

    for(i=0;i<n;i++)
    {
        X[i] = 0;
    }

    //TODO attempt to make it diagonally dominant.

    if (!isDiagonallyDominant(A,n))
    {
        key=0;
        printf("Not diagonally dominant, method fails->");
    }
    key = 1;

    while(key==1)
    {

        for(i=0;i<n;i++)
        {
            sum = B[i];
            for(j=0;j<n;j++)
            {
                if(i!=j)
                {
                  sum -= A[i][j]*X[j];
                }
            }
            X1[i] = sum/A[i][i];
            calErr = (X1[i]-X[i])/X1[i];
            if(calErr<0) calErr *= (-1);
            if(calErr>err)
            {
                key = 1;
                X[i] = X1[i];
            }
            else
            {
                key = 0;
            }
        }
    }
    return X1;
}

void generateSolution(int combination[],double **A,double *B,int start,int index,int m,int n,double err)
{
    int j=0,k=0,i;
    double *x;
    double *solution;
    x=(double *)malloc(m*sizeof(double));
    solution=(double *)malloc(n*sizeof(double));

    if (index==m)
    {
        double **A_square;
        A_square=(double ** )malloc(m*sizeof(double * ));
        for (j=0;j<m;j++)
            {
                A_square[j]=(double *)malloc(m*sizeof(double));
            }
        for (j=0;j<m;j++)
        {
            for (k=0;k<m;k++)
            {

                A_square[j][k]=A[j][combination[k]];
            }
        }

        x = gaussSeidel(A_square,B,m,err);

        for(j=0;j<n;j++)
        {
            solution[j]=0;
        }
        for (k=0;k<m;k++)
        {
            solution[combination[k]]=x[k];
        }
        printf("[");
        for (j=0;j<n;j++)
        {
            printf("%.3lf ",solution[j]);
        }
        printf("]\n");
        return;
    }

    for (i=start; i<=n-1 && n-i >= m-index; i++)
    {
        combination[index] = i;
        generateSolution(combination,A,B,i+1,index+1,m,n,err);
    }
}

int main()
{

	int n,m,i,j,*combination;
    float err=0.001;
	double **A,*B;

	printf("Enter the value of N: ");
	scanf("%d",&n);
    printf("Enter the value of M: ");
    scanf("%d",&m);
    printf("Enter the value of error permissible: ");
    scanf("%f",&err);

   	A = (double **)malloc(m *sizeof(double *));
    for (i=0; i<m; i++)
    {
    	A[i] = (double *)malloc(n * sizeof(double));
    }

    B = (double *)malloc(m*sizeof(double));

    combination = (int *)malloc(m*sizeof(int));

    printf("Enter MATRIX A\n");
    for(i=0;i<m;i++)
    {
    	for(j=0;j<n;j++)
    	{
    		printf("A[%d][%d]=",i,j);
    		scanf("%lf",&A[i][j]);
    	}
    }

    printf("Enter MATRIX B\n");
    for(i=0;i<m;i++)
    {
    	printf("B[%d]=",i);
		scanf("%lf",&B[i]);
    }
    printf("Printing Basic solutions \n");
    generateSolution(combination,A,B,0,0,m,n,err);

	return 0;
}