#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define bM 100000000
struct node
{
    double val;
    bool isInList;
    int index;
};
void findVeqn(vector<pair<node, node>> totEquations, int n, double *b, node *x)
{

    int iter = 0, indx;
    //we know that x[0] = 0, after that, we have to find out x[i], for i = 1,2,..., n
    while (iter++ < 10)
    {

        for (int i = 0; i < n; i++)
        {
            if (totEquations[i].first.isInList)
            {
                totEquations[i].second.val = b[i] - totEquations[i].first.val;
                totEquations[i].second.isInList = true;
                indx = totEquations[i].second.index;
                x[indx - 1].val = totEquations[i].second.val;
                for (int j = 0; j < n; j++)
                {
                    if (indx == totEquations[j].first.index)
                    {
                        totEquations[j].first.val = totEquations[i].second.val;
                        totEquations[j].first.isInList = true;
                        x[indx - 1].val = totEquations[i].second.val;
                    }
                    else if (indx == totEquations[j].second.index)
                    {
                        totEquations[j].second.val = totEquations[i].second.val;
                        totEquations[j].second.isInList = true;
                        x[indx - 1].val = totEquations[i].second.val;
                    }
                }
            }
            else if (totEquations[i].second.isInList)
            {
                totEquations[i].first.val = b[i] - totEquations[i].second.val;
                totEquations[i].first.isInList = true;
                indx = totEquations[i].first.index;
                x[indx -1].val = totEquations[i].first.val;
                for (int j = 0; j < n; j++)
                {
                    if (indx == totEquations[j].first.index)
                    {
                        totEquations[j].first.val = totEquations[i].first.val;
                        totEquations[j].first.isInList = true;
                        x[indx -1].val = totEquations[i].first.val;
                    }
                    else if (indx == totEquations[j].second.index)
                    {
                        totEquations[j].second.val = totEquations[i].first.val;
                        totEquations[j].second.isInList = true;
                        x[indx -1].val = totEquations[i].first.val;
                    }
                }
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        cout<<"x_"<<i<<" = "<<x[i].val<<"\n";
        if (totEquations[i].first.isInList && totEquations[i].second.isInList)
            cout << "x_" << totEquations[i].first.index << " + x_" << totEquations[i].second.index << " = " << totEquations[i].first.val << " + " << totEquations[i].second.val << "\n";
    }
}
int main()
{

    node *x = new node[6];
    for (int i = 0; i < 6; i++)
    {
        x[i].isInList = false;
        x[i].val = 0.0;
    }
    pair<node, node> eqn;
    vector<pair<node, node>> totEquations;
    eqn.first.index = 0;
    eqn.first.isInList = true;
    eqn.first.val = 0.0;
    eqn.second.index = 3;
    eqn.second.isInList = false;
    totEquations.push_back(eqn);

    eqn.first.index = 0;
    eqn.first.isInList = true;
    eqn.first.val = 0.0;
    eqn.second.index = 4;
    eqn.second.isInList = false;
    totEquations.push_back(eqn);

    eqn.first.index = 1;
    eqn.first.isInList = false;
    eqn.second.index = 4;
    eqn.second.isInList = false;
    totEquations.push_back(eqn);

    eqn.first.index = 1;
    eqn.first.isInList = false;
    eqn.second.index = 6;
    eqn.second.isInList = false;
    totEquations.push_back(eqn);

    eqn.first.index = 2;
    eqn.first.isInList = false;
    eqn.second.index = 5;
    eqn.second.isInList = false;
    totEquations.push_back(eqn);

    eqn.first.index = 2;
    eqn.first.isInList = false;
    eqn.second.index = 6;
    eqn.second.isInList = false;
    totEquations.push_back(eqn);

    double b[] = {10, 5, 2, 6, 4, 8};
    findVeqn(totEquations, 6, b, x);
    return 0;
}