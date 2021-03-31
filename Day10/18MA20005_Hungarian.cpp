/*
Name: Altaf Ahmad
Roll no: 18MA20005
Assingment problem using the  The Hungarian algorithm
*/
#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define bM 1000000.0
void printTable(int **arr, int n) // helper function to print the table 
{
    cout << "\n";
    for (int i = 0; i < n; i++)
    {
        cout << "---------";
    }
    cout << "\n  |";
    for (int i = 1; i <= n; i++)
    {
        cout << "\t" << i;
    }
    cout << "\n";
    for (int i = 0; i < n; i++)
    {
        cout << "---------";
    }
    cout << "\n";
    for (int i = 0; i < n; i++)
    {
        cout << i + 1 << " |\t";
        for (int j = 0; j < n; j++)
        {
            cout << arr[i][j] << "\t";
        }
        cout << "\n";
    }
    cout << "\n";
    for (int i = 0; i < n; i++)
    {
        cout << "---------";
    }
    cout << "\n";
}
int hungarian(int **arr, int n) // method to implement the  hungarian method
{
    int **A = new int *[100];
    for (int i = 0; i < 100; i++)
    {
        A[i] = new int[100];
    }
    cout << "\nInitial Table : \n";
    printTable(arr, n);
    int *rowmin = new int[100];
    int *columnmin = new int[100];
    int minpk = bM;
    for (int i = 0; i < n; i++)
    {
        minpk = bM;
        for (int j = 0; j < n; j++)
        {
            A[i][j] = arr[i][j];
            if (A[i][j] < minpk)
            {
                minpk = A[i][j];
            }
        }
        rowmin[i] = minpk;
    }
    //subtracting row minimum
    cout<<"The row minima is : ";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] -= rowmin[i];
            
        }
        cout<<"\t" << rowmin[i];
    }
    cout << "\nAfter subtracting row minima we get\n";
    printTable(A, n);
    //subtracting column minima
    cout<<"The column minima is : ";
    for (int j = 0; j < n; j++)
    {
        minpk = bM;
        for (int i = 0; i < n; i++)
        {
            if (A[i][j] < minpk)
            {
                minpk = A[i][j];
            }
        }
        cout << "\t" << minpk;
        rowmin[j] = minpk;
    }
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            A[i][j] -= rowmin[j];
        }
    }
    cout << "\n\nAfter subtracting column minima we get\n";
    printTable(A, n);
    int iter = 0;
    int **operator1 = new int *[100];
    for (int i = 0; i < 100; i++)
    {
        operator1[i] = new int[100];
    }

    while (1)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                operator1[i][j] = 0;
            }
        }
        //covering all the zeroes with minimum number of lines
        //first row scanning is done

        if (iter++ > 3)
            return 0;
        int numzeroes, numsquares = 0;
        vector<int> horline;
        vector<int> verline;
        int **checkmatrix = new int *[100];
        for (int i = 0; i < 100; i++)
        {
            checkmatrix[i] = new int[100];
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                checkmatrix[i][j] = 0;
            }
        }
        for (int i = 0; i < n; i++)
        {
            numzeroes = 0;
            for (int j = 0; j < n; j++)
            {
                if (checkmatrix[i][j])
                {
                    continue;
                }
                if (A[i][j] == 0)
                {
                    numzeroes++;
                }
            }
            if (numzeroes == 1)
            {

                for (int j = 0; j < n; j++)
                {
                    if (checkmatrix[i][j])
                    {
                        continue;
                    }
                    if (A[i][j] == 0)
                    {
                        numsquares++;
                        operator1[i][j] = 1;
                        verline.push_back(j);
                        for (int lm = 0; lm < n; lm++)
                        {
                            checkmatrix[lm][j] = 1;
                        }
                        break;
                    }
                }
            }
        }
        //now do the column scanning
        for (int j = 0; j < n; j++)
        {
            numzeroes = 0;
            for (int i = 0; i < n; i++)
            {
                if (checkmatrix[i][j])
                    continue;
                if (A[i][j] == 0)
                {
                    numzeroes++;
                }
            }
            if (numzeroes == 1)
            {
                for (int i = 0; i < n; i++)
                {
                    if (checkmatrix[i][j])
                    {
                        continue;
                    }
                    if (A[i][j] == 0)
                    {
                        numsquares++;
                        operator1[i][j] = 1;
                        horline.push_back(i);
                        for (int lm = 0; lm < n; lm++)
                        {
                            checkmatrix[i][lm] = 1;
                        }
                        break;
                    }
                }
            }
        }
        numzeroes = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (checkmatrix[i][j])
                {
                    continue;
                }
                if (A[i][j] == 0)
                    numzeroes++;
            }
        }
        if (numzeroes)
        {
            
            for (int i = 0; i < n; i++)
            {
                //if(i == 2)
                 //   continue;
                numzeroes = 0;
                for (int j = 0; j < n; j++)
                {
                    if (checkmatrix[i][j])
                    {
                        continue;
                    }
                    if (A[i][j] == 0)
                    {
                        numzeroes++;
                        numsquares++;
                        operator1[i][j] = 1;
                        verline.push_back(j);
                        for (int lm = 0; lm < n; lm++)
                        {
                            checkmatrix[lm][j] = 1;
                        }
                        break;
                    }
                }
            }
            //now do the column scanning
            for (int j = 0; j < n; j++)
            {
                numzeroes = 0;
                for (int i = 0; i < n; i++)
                {
                    if (checkmatrix[i][j])
                        continue;
                    if (A[i][j] == 0)
                    {
                        numzeroes++;
                        numsquares++;
                        operator1[i][j] = 1;
                        horline.push_back(i);
                        for (int lm = 0; lm < n; lm++)
                        {
                            checkmatrix[i][lm] = 1;
                        }
                        break;
                    }
                }
            }
            //numsquares = n;
        }
        cout << "The check matrix after covering all the zeroes is : \n";
        for (int i = 0; i < n; i++)
        {

            for (int j = 0; j < n; j++)
            {
                cout << checkmatrix[i][j] << "\t";
            }
            cout << "\n";
        }

        cout << "\n";
        if (numsquares == n)
        {
            cout << "Optimality has reached\nThe operator is placed at the places where there is 1\n";
            printTable(operator1, n);
            int total = 0;
            cout << "That is \n";
            for (int i = 0; i < n; i++)
            {
                cout << "For operator " << i + 1 << ", job = ";
                for (int j = 0; j < n; j++)
                {
                    if (operator1[i][j])
                    {
                        total += arr[i][j];
                        cout << j + 1 << "\n";
                    }
                }
            }
            cout << "And the total cost is : " << total << "\n";
            return total;
        }
        else // otherwise, the minimum non covered value is chosen 
        {
            int cx, cy;
            minpk = bM;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (checkmatrix[i][j])
                        continue;
                    if (A[i][j] < minpk)
                    {
                        minpk = A[i][j]; // this finds the minimum non covered value
                    }
                }
            }
            cout<<"The minimum uncovered element is : "<<minpk<<"\nThis will be subtracted and added in the required cells now\n";
            for (int i = 0; i < horline.size(); i++)
            {
                for (int j = 0; j < verline.size(); j++)
                {
                    cout << "verline[j]][horline[i] = " << verline[j] << " ," << horline[i] << " \n";
                    A[horline[i]][verline[j]] += minpk; // the places where hor and vertical lines are drawn, the minimum value is added
                }
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (checkmatrix[i][j]) // while it is reduced in the non covered cases
                        continue;
                    A[i][j] -= minpk;
                }
            }
            printTable(A, n);
        }
    }
}
int main()
{
    int n;
    cout << "Enter the number of rows and columns, i.e, the number of workers and tasks : ";
    cin >> n;
    int **A = new int *[100];
    for (int i = 0; i < 100; i++)
    {
        A[i] = new int[100];
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
    int cost = hungarian(A, n);
}