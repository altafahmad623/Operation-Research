#include <bits/stdc++.h>
using namespace std;

#define trace(x) cerr << #x << ": " << x << " " << endl;

int pivotRow, pivotCol;
double pivotElem;
bool typeOfLPP;
int p,q;
bool maxi = 0;

void calculateC(vector<double> &basisCoef, vector<vector<double> > & Ar, vector<double> &C, int n, int m)
{
    for(int i=0;i<n+p+q;i++)
    {
        C[i] = 0;
        for(int j=0;j<m;j++)
        {
            C[i] += Ar[j][i]*basisCoef[j];
        }
    }
    cout<<"\tC[i]: ";
    for(int i=0;i<n+p+q;i++)
    {
        cout<<C[i]<<" ";
    }
    cout<<endl;
}

vector<double> calculateDiffZC(vector<double> Z, vector<double> C, int n, int m)
{
    vector<double> diff(n+p+q);
    for(int i=0;i<n+p+q;i++)
    {
        diff[i] = Z[i]-C[i];
    }
    cout<<"\tZ-C[i]: ";
    for(int i=0;i<n+p+q;i++)
    {
        cout<<diff[i]<<" ";
    }
    cout<<endl;
    return diff;
}

int findPivotCol(vector<double> diffZC, int n, int m)
{
    if(typeOfLPP)
    {
        double maxDif = diffZC[0];
        int ind = 0;
        for(int i=1;i<n+p+q;i++)
        {
            if(diffZC[i] > maxDif)
            {
                maxDif = diffZC[i];
                ind = i;
            }
        }
        return ind;
    }
    else
    {
        double minDif = diffZC[0];
        int ind = 0;
        for(int i=1;i<n+p+q;i++)
        {
            if(diffZC[i] < minDif)
            {
                minDif = diffZC[i];
                ind = i;
            }
        }
        return ind;
    }
}

int findPivotRow(vector<vector<double> > & Ar, vector<double> B, int n, int m)
{
    int ind = -1;
    double minRatio = 1e9;
    if(minRatio)
    for(int i=0;i<m;i++)
    {
        if(B[i]/Ar[i][pivotCol] >= 0 && B[i]/Ar[i][pivotCol] < minRatio)
        {
            minRatio = B[i]/Ar[i][pivotCol];
            ind = i;
        }
    }
    return ind;
}

void updateSimplexTable(vector<vector<double> > &Ar, vector<double> &B, double pivotElem, int n, int m)
{
    vector<vector<double> > Ar2(m);
    for(int i=0;i<m;i++)
    {
        Ar2[i] = vector<double>(n+p+q);
    }
    vector<double> B2(m);
    for(int i=0;i<m;i++)
    {
        if(i == pivotRow)
        {
            for(int j=0;j<n+p+q;j++)
            {
                Ar2[i][j] = Ar[i][j]/pivotElem;
            }
            B2[i] = B[i]/pivotElem;
        }
        else
        {
            for(int j=0;j<n+p+q;j++)
            {
                Ar2[i][j] = Ar[i][j] - Ar[i][pivotCol]*Ar[pivotRow][j]/pivotElem;
            }
            B2[i] = B[i]-B[pivotRow]*Ar[i][pivotCol]/pivotElem;
        }
    }
    B = B2;
    Ar = Ar2;
}

double findObjectiveValue(vector<double> basisCoef, vector<double> B, int m)
{
    double val=0;
    for(int i=0;i<m;i++)
    {
        val += basisCoef[i]*B[i];
    }
    if(maxi)    val *= -1;
    return val;
}

bool performIteration(vector<vector<double> > &Ar, vector<double> &B, vector<int> &basis, vector<double> &basisCoef, vector<double> &Z, vector<double> &C, int n, int m)
{
    calculateC(basisCoef, Ar, C, n, m);

    vector<double> diffZC = calculateDiffZC(Z, C, n, m);

    pivotCol = findPivotCol(diffZC, n, m);

    if(typeOfLPP && diffZC[pivotCol]<1e-6)
    {
        return 1;
    }
    else if(!typeOfLPP && diffZC[pivotCol]>-1e-6)
    {
        return 1;
    }

    pivotRow = findPivotRow(Ar, B, n, m);
    if(pivotRow == -1)
    {
        return 1;
    }

    double minRatio = B[pivotRow]/Ar[pivotRow][pivotCol];
    basis[pivotRow] = pivotCol+1;
    basisCoef[pivotRow] = Z[pivotCol];

    pivotElem = Ar[pivotRow][pivotCol];

    updateSimplexTable(Ar, B, pivotElem, n, m);

    return 0;
}

void iterationPrint(vector<vector<double> > &Ar, vector<double> B, vector<int> basis, vector<double> basisCoef, int m, int n, vector<double> Z)
{
    cout<<"\n\tPrinting Basis Variables: ";
    for(int i=0;i<m;i++)
    {
        if(basis[i]<=n)
        {
            cout<<"X"<<basis[i]<<" ";
        }
        else if(basis[i] <= n+p)
        {
            cout<<"S"<<basis[i]-n<<" ";
        }
        else
        {
            cout<<"A"<<basis[i]-n-p<<" ";
        }
    }
    cout<<endl;

    cout<<"\n\tPrinting Basis Coefficients: ";
    for(int i=0;i<m;i++)
    {
        cout<<basisCoef[i]<<" ";
    }
    cout<<endl;

    cout<<"\n\tPrinting Basis Matrix: \n";
    for(int i=0;i<m;i++)
    {
        cout<<"\t";
        if(basis[i]<=n)
        {
            cout<<"X"<<basis[i];
        }
        else if(basis[i] <= n+p)
        {
            cout<<"S"<<basis[i]-n;
        }
        else
        {
            cout<<"A"<<basis[i]-n-p;
        }
        cout<<"\t";
        for(int j=0;j<n+p+q;j++)
        {
            cout<< Ar[i][j]<<" ";
        }
        cout<<B[i]<<" ";
        cout<<B[i]/Ar[i][pivotCol]<<endl;
    }

    cout<<"\n\tPrinting objective value: ";
    cout<<findObjectiveValue(basisCoef, B, m);
    cout<<"\n\n";
}

void printLPP(vector<vector<double> > &Ar, vector<double> B, vector<int> slack, vector<int> artificial, vector<double> Z, bool typeofLPP, vector<string> sign, int n, int m)
{
    int p = slack.size();
    cout<<"LPP:- \n";
    cout<<(typeOfLPP?"Max: ":"Min: ")<<"Z = ";
    cout<<Z[0]<<"X1";
    for(int i=1;i<n+p+q-1;i++)
    {
        cout<<(Z[i]>=0?"+":"")<<Z[i];
        if(i<n)
        {
            cout<<"X"<<i+1;
        }
        else if(i<n+p)
        {
            cout<<"S"<<i-n+1;
        }
        else
        {
            cout<<"A"<<i-n-p+1;
        }
    }
    cout<<endl;
    for(int i=0;i<m;i++)
    {
        cout<<"\t"<<Ar[i][0]<<"X1";
        for(int j=1;j<n;j++)
        {
            cout<<(Ar[i][j]>=0?"+":"")<<Ar[i][j]<<"X"<<j+1;
        }
        for(int j=0;j<p;j++)
        {
            if(slack[j]==i)
            {
                cout<<(Ar[i][n+j]>=0?"+":"")<<Ar[i][n+j]<<"S"<<j+1;
            }
            if(slack[j]==-1-i)
            {
                cout<<(Ar[i][n+j]>=0?"+":"")<<Ar[i][n+j]<<"S"<<j+1;
            }
        }
        for(int j=0;j<q;j++)
        {
            if(artificial[j] == i)
            {
                cout<<(Ar[i][n+p+j]>=0?"+":"")<<Ar[i][n+p+j]<<"A"<<j+1;
            }
        }
        cout<<" "<<sign[i]<<" "<<B[i]<<endl;
    }
    cout<<endl;
}

int main()
{
    int n,m;

    cout<<"Enter the number of equations m\n";
    cin>>m;

    cout<<"Enter the number of variables n\n";
    cin>>n;
    vector<vector<double> > Ar(m);
    for(int i=0;i<m;i++)
    {
        Ar[i] = vector<double>(n);
    }

    vector<double> B(m);
    vector<double> Z(n);
    vector<double> C(n+2*m);
    vector<int> slack;
    vector<int> artificial;
    vector<string> sign(m);
    maxi = 0;

    cout<<"Enter the type of LPP. 0 for minimization type, 1 for maximization type\n";
    cin>>typeOfLPP; // 0 - min type; 1 - max type

    cout<<"Enter the coefficients of the LPP\n";
    for(int i=0;i<n;i++)
    {
        cin>>Z[i];
        if(typeOfLPP)
        {
            maxi = 1;
            Z[i] *= -1;
        }
    }
    typeOfLPP = 0;

    for(int i=0;i<m;i++)
    {
        cout<<"Enter the coefficients of the "<<i+1<<"th equation\n";
        for(int j=0;j<n;j++)
        {
            cin>>Ar[i][j];
        }

        cout<<"Enter the sign of the inequality (>= or <= or =)\n";
        cin>>sign[i];
        if(sign[i] == "<=")
        {
            slack.push_back(i);
        }
        else if (sign[i] == ">=")
        {
            slack.push_back(-i-1);
            artificial.push_back(i);
        }
        else
        {
            artificial.push_back(i);
        }

        cout<<"Enter the RHS value\n";
        cin>>B[i];
    }
    p = slack.size();
    q = artificial.size();

    for(auto it:slack)
    {
        Z.push_back(0);
        if(it>=0)
        {
            for(int i=0;i<m;i++)
            {
                if(i == it)
                {
                    Ar[i].push_back(1);
                }
                else
                {
                    Ar[i].push_back(0);
                }
            }
        }
        else
        {
            for(int i=0;i<m;i++)
            {
                if(-1-i == it)
                {
                    Ar[i].push_back(-1);
                }
                else
                {
                    Ar[i].push_back(0);
                }
            }
        }
    }

    for(auto it:artificial)
    {
        if(typeOfLPP)
        {
            Z.push_back(-1);
        }
        else
        {
            Z.push_back(1);
        }
        for(int i=0;i<m;i++)
        {
            if(i == it)
                Ar[i].push_back(1);
            else
                Ar[i].push_back(0);
        }
    }

    printLPP(Ar, B, slack, artificial, Z, typeOfLPP, sign, n, m);

    vector<int> basis(m);
    vector<double> basisCoef(m);
    for(int i=0;i<m;i++)
    {
        basis[i]= -1;
        for(int j=0;j<q;j++)
        {
            if(artificial[j] == i)
            {
                basis[i] = n+p+j+1;
                basisCoef[i] = Z[n+p+j];
                break;
            }
        }
        if(basis[i]!= -1)   continue;

        for(int j=0;j<p;j++)
        {
            if(slack[j] == i || slack[j]== -1-i)
            {
                basis[i] = n+j+1;
                basisCoef[i] = Z[n+j];
                break;
            }
        }
    }
    cout<<"PHASE 1\n\n";

    vector<double> Z_0(Z.size());
    for(int i=0;i<n+p;i++)
    {
        Z_0[i] = 0.0;
    }
    for(int i=n+p;i<n+p+q;i++)
    {
        Z_0[i] = 1.0;
    }
    iterationPrint(Ar, B, basis, basisCoef, m, n, Z_0);
    bool optimal=0;
    int iter=1;
    while(!optimal)
    {
        optimal = performIteration(Ar, B, basis, basisCoef, Z_0, C, n, m);
        if(optimal)
        {
            // cout<<"Optimal Value achieved!!\n";
            break;
        }
        cout<<"\nPerforming iteration "<<iter++<<endl;
        iterationPrint(Ar, B, basis, basisCoef, m, n, Z_0);
    }
    bool arti = 0;
    for(int i=0;i<m;i++)
    {
        if(basis[i] > n+p)
        {
            arti=1;
        }
    }
    if(arti)
    {
        cout<<"Artificial variable found in the final table. The solution does not converge.\n";
        // cout<<"Value of Z = "<<findObjectiveValue(basisCoef, B, m)<<endl;
        return 0;
    }

    q=0;    // removing artificial variables
    cout<<"\nPHASE 2\n";
    for(int i=0;i<m;i++)
    {
        basisCoef[i] = Z[basis[i]-1];
    }
    iterationPrint(Ar, B, basis, basisCoef, m, n, Z);

    optimal=0;
    iter=1;
    while(!optimal)
    {
        optimal = performIteration(Ar, B, basis, basisCoef, Z, C, n, m);
        cout<<"\nPerforming iteration "<<iter++<<endl;
        iterationPrint(Ar, B, basis, basisCoef, m, n, Z);
        if(optimal)
        {
            cout<<"Optimal Value achieved!!\n";
            break;
        }
    }

    cout<<"\nFinal Solution: \n";
    for(int i=0;i<m;i++)
    {
        if(basis[i]<=n)
        {
            cout<<"X"<<basis[i];
        }
        else if(basis[i]<=n+p)
        {
            cout<<"S"<<basis[i]-n;
        }

        cout<<" = "<<B[i]<<" ";
    }
    cout<<" rest all 0.0\n";
    cout<<"Optimum value of Z = "<<findObjectiveValue(basisCoef, B, m)<<endl;
}
