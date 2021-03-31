#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define bM 1000000.0
int main()
{
    int n;
    cout << "Enter the number of rows and columns, i.e, the number of workers and tasks : ";
    cin >> n;
    int **a1 = new int *[10];
    for (int i = 0; i < 10; i++)
    {
        a1[i] = new int[10];
    }
    for (int i = 0; i < n; i++)
    {
        for(int j =0 ; j < n ; j++)
        {
            cin>>a1[i][j];
        }
    }
     int **a2 = new int *[10];
    for (int i = 0; i < 10; i++)
    {
        a2[i] = new int[10];
    }
    for (int i = 0; i < n; i++)
    {
        for(int j =0 ; j < n ; j++)
        {
            cin>>a2[i][j];
        }
    }
     int **a3 = new int *[10];
    for (int i = 0; i < 10; i++)
    {
        a3[i] = new int[10];
    }
    for (int i = 0; i < n; i++)
    {
        for(int j =0 ; j < n ; j++)
        {
            cin>>a3[i][j];
        }
    }
    int sum = 0;
    for (int i = 0; i < n; i++)
    {
        cout<<"\n";
        for(int j =0 ; j < n ; j++)
        {
            sum = a1[i][j] + a2[i][j] + a3[i][j];
            cout<<" "<<sum;
        }
    }
    return 0;
}
