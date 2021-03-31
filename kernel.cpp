#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define bM 1000000.0
int main()
{
    int n,sum;
    cout << "Enter the number of rows and columns, i.e, the number of workers and tasks : ";
    cin >> n;
    int **A = new int *[100];
    for (int i = 0; i < 100; i++)
    {
        A[i] = new int[100];
    }
    int **out = new int *[100];
    for (int i = 0; i < 100; i++)
    {
        out[i] = new int[100];
    }
    int **ker = new int *[100];
    for (int i = 0; i < 100; i++)
    {
        ker[i] = new int[100];
    }
    cout << "\nNow, enter the values of the table : \n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("Enter A[%d][%d] :", i, j);
            cin >> A[i][j];
        }
    }
    cout << "\nNow, enter the kernel table : \n";
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            printf("Enter ker[%d][%d] :", i, j);
            cin >> ker[i][j];
        }
    }
    cout<<"\n\n";
    for (int i = 1; i < (n-1); i++)
    {
        for (int j = 1; j < (n-1); j++)
        {
            sum = (ker[0][0]*A[i-1][j-1])+(ker[0][1]*A[i-1][j])+(ker[0][2]*A[i-1][j+1]);
            sum +=(ker[1][0]*A[i][j-1])+(ker[1][1]*A[i][j])+(ker[1][2]*A[i][j+1]);
            sum +=(ker[2][0]*A[i+1][j-1])+(ker[2][1]*A[i+1][j])+(ker[2][2]*A[i+1][j+1]);
            out[i][j] = sum;
            cout<<" "<<out[i][j];
        }
        cout<<"\n";
    }
    
    return 0;
}