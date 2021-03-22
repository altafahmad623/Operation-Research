#include<bits/stdc++.h>
using namespace std;

void basicFeasSol(int x, int y, int n, int m, vector<vector<int>> &val, vector<int> sup, vector<int> dem)
{
    if(x==n || y==m)
    {
        return;
    }
    if(x==n-1)
    {
        val[x][y] = dem[y];
        sup[x] -= dem[y];
        dem[y] = 0;
        basicFeasSol(x,y+1,n,m,val,sup,dem);
        return;
    }
    if(y==m-1)
    {
        val[x][y] = sup[x];
        dem[y] -= sup[x];
        sup[x] = 0;
        basicFeasSol(x+1,y,n,m,val,sup,dem);
        return;
    }

    if(sup[x]>dem[y])
    {
        val[x][y] = dem[y];
        sup[x] -= dem[y];
        dem[y] = 0;
        basicFeasSol(x,y+1,n,m,val,sup,dem);
    }
    else
    {
        val[x][y] = sup[x];
        dem[y] -= sup[x];
        sup[x] = 0;
        basicFeasSol(x+1,y,n,m,val,sup,dem);
    }
}

void printVal(vector<vector<int>> val)
{
    cout<<"Printing the ?????? matrix\n";
    for(auto row:val)
    {
        for(auto it:row)
        {
            cout<<it<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;
}

int optimalCost(vector<vector<int>> val, vector<vector<int>> cost, int n, int m)
{
    int sum = 0;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            sum += val[i][j]*cost[i][j];
        }
    }
    return sum;
}

void printPoint(pair<int,int> pt)
{
    cout<<"("<<pt.first+1<<", "<<pt.second+1<<")";
}

// vector<pair<int,int>> findCycle(vector<vector<int>> val, pair<int,int> pvt, int n, int m)
// {
//     vector<pair<int,int>> cycle;
//     set<pair<int,int>> vis;
//     cycle.push_back(pvt);

//     while(!vis.count(pvt))
//     {
//         pair<int,int> cur = cycle[cycle.size()-1];
//         if(cur.first > 0 && val[cur.first-1][cur.second])
//         {
//             if(cycle.size() != 2 && cur.first-1 == pvt.first && cur.second == pvt.second)
//             {
                
//             }
//         }
//     }


// }

vector<pair<int,int>> findCycle(vector<vector<int>> val, pair<int,int> pvt, int n, int m)
{
    vector<pair<int,int>> cycle;
    cycle.push_back(pvt);

    if(pvt.first>0 && pvt.second>0)
    {
        if(val[pvt.first][pvt.second-1] && val[pvt.first-1][pvt.second] && val[pvt.first-1][pvt.second-1])
        {
            cycle.push_back({pvt.first,pvt.second-1});
            cycle.push_back({pvt.first-1,pvt.second-1});
            cycle.push_back({pvt.first-1,pvt.second});
            return cycle;
        }
    }
    if(pvt.first>0 && pvt.second<m-1)
    {
        if(val[pvt.first][pvt.second+1] && val[pvt.first-1][pvt.second] && val[pvt.first-1][pvt.second+1])
        {
            cycle.push_back({pvt.first,pvt.second+1});
            cycle.push_back({pvt.first-1,pvt.second+1});
            cycle.push_back({pvt.first-1,pvt.second});
            return cycle;
        }
    }
    if(pvt.first<n-1 && pvt.second<m-1)
    {
        if(val[pvt.first][pvt.second+1] && val[pvt.first+1][pvt.second] && val[pvt.first+1][pvt.second+1])
        {
            cycle.push_back({pvt.first,pvt.second+1});
            cycle.push_back({pvt.first+1,pvt.second+1});
            cycle.push_back({pvt.first+1,pvt.second});
            return cycle;
        }
    }
    if(pvt.first<n-1 && pvt.second>0)
    {
        if(val[pvt.first][pvt.second-1] && val[pvt.first+1][pvt.second] && val[pvt.first+1][pvt.second-1])
        {
            cycle.push_back({pvt.first,pvt.second-1});
            cycle.push_back({pvt.first+1,pvt.second-1});
            cycle.push_back({pvt.first+1,pvt.second});
            return cycle;
        }
    }

    if(cycle.size() == 1)
    {
        cout<<"Cycle not found!!";
    }
    return cycle;
}

void updateTable(vector<pair<int,int>> cycle, vector<vector<int>> &val)
{
    int vt = val[cycle[1].first][cycle[1].second];
    for(int i=3;i<cycle.size();i+=2)
    {
        vt = min(vt,val[cycle[i].first][cycle[i].second]);
    }
    int sign = -1;
    for(int i=0;i<cycle.size();i++)
    {
        val[cycle[i].first][cycle[i].second] -= sign*vt;
        sign *= -1;
    }
}

int main()
{
    int n,m;

    cout<<"Enter the number of sources: ";
    cin>>n;
    cout<<"Enter the number of destinations: ";
    cin>>m;

    vector<int> sup(n);
    vector<int> dem(m);

    cout<<"Enter the supply from each source: ";
    for(int i=0;i<n;i++)
    {
        cin>>sup[i];
    }
    cout<<"Enter the demand of each destination: ";
    for(int i=0;i<m;i++)
    {
        cin>>dem[i];
    }

    int sum=0;
    for(auto it:sup)    sum+=it;
    for(auto it:dem)    sum-=it;

    if(sum)
    {
        cout<<"The given problem is infeasible.\n";
        return 0;
    }

    vector<vector<int>> cost(n,vector<int>(m));;
    vector<vector<int>> val(n,vector<int>(m,0));;
    vector<vector<int>> C(n,vector<int>(m,0));;

    for(int i=0;i<n;i++)
    {
        cout<<"Enter the cost of transportation from source "<<i+1<<" to all destinations\n\n";
        for(int j=0;j<m;j++)
        {
            cin>>cost[i][j];
        }
    }

    cout<<"Performing North-west corner method to find basic feasible solution\n";

    basicFeasSol(0,0,n,m,val,sup,dem);

    printVal(val);

    bool optimal=0;
    int iter=1;

    while(!optimal)
    {
        cout<<"Iteration "<<iter++<<endl<<endl;

        vector<int> u(n,0), v(m,0);
        vector<bool> us(n,0), vs(m,0);

        u[0] = 0;
        us[0] = 1;
        int cnt=1;
        while(cnt < n+m)
        {
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<m;j++)
                {
                    if(val[i][j])
                    {
                        if(us[i])
                        {
                            v[j] = cost[i][j]-u[i];
                            vs[j]=1;
                            cnt++;
                        }
                        else if(vs[j])
                        {
                            u[i] = cost[i][j]-v[j];
                            us[i]=1;
                            cnt++;
                        }
                    }
                }
            }
        }

        cout<<"\tPrinting values of u[i]: ";
        for(auto it:u)  cout<<it<<" ";
        cout<<"\n\n";

        cout<<"\tPrinting values of v[j]: ";
        for(auto it:v)  cout<<it<<" ";
        cout<<"\n\n";

        optimal = 1;
        int curMax = 0;
        pair<int,int> ind={-1,-1};

        for(int i=0;i<n;i++)
        {
            for(int j=0;j<m;j++)
            {
                if(!val[i][j])
                    C[i][j] = u[i]+v[j]-cost[i][j];
                else
                    C[i][j] = 0;
                if(C[i][j]>0)
                {
                    optimal=0;
                    if(C[i][j] > curMax)
                    {
                        curMax = C[i][j];
                        ind = {i,j};
                    }
                }
            }
        }

        cout<<"\tPrinting penalty matrix:\n";
        for(int i=0;i<n;i++)
        {
            cout<<"\t";
            for(int j=0;j<m;j++)
            {
                cout<<C[i][j]<<"\t";
            }
            cout<<endl;
        }
        cout<<endl;

        if(optimal) break;

        cout<<"\tNew basis point is "<<ind.first+1<<", "<<ind.second+1<<endl;

        vector<pair<int,int>> cycle = findCycle(val, ind, n, m);
        cout<<"\tCycle found: ";
        for(auto it:cycle)
        {
            printPoint(it);
            cout<<"->";
        }
        printPoint(ind);
        cout<<"\n\n";

        updateTable(cycle, val);
        printVal(val);
    }

    cout<<"Optimality reached! \n";

    printVal(val);
    cout<<"Minimum transportation cost is "<<optimalCost(val,cost,n,m)<<endl;

}