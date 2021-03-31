#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;
#pragma region DEBUG
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
typedef
tree<
  pair<int,int>,
  null_type,
  less<pair<int,int>>,
  rb_tree_tag,
  tree_order_statistics_node_update>
ordered_set;
void __print(int x) {cerr << x;}
void __print(long x) {cerr << x;}
void __print(long long x) {cerr << x;}
void __print(unsigned x) {cerr << x;}
void __print(unsigned long x) {cerr << x;}
void __print(unsigned long long x) {cerr << x;}
void __print(float x) {cerr << x;}
void __print(double x) {cerr << x;}
void __print(long double x) {cerr << x;}
void __print(char x) {cerr << '\'' << x << '\'';}
void __print(const char *x) {cerr << '\"' << x << '\"';}
void __print(const string &x) {cerr << '\"' << x << '\"';}
void __print(bool x) {cerr << (x ? "true" : "false");}

template<typename T, typename V>
void __print(const pair<T, V> &x) {cerr << '{'; __print(x.first); cerr << ','; __print(x.second); cerr << '}';}
template<typename T>
void __print(const T &x) {int f = 0; cerr << '{'; for (auto &i: x) cerr << (f++ ? "," : ""), __print(i); cerr << "}";}
void _print() {cerr << "]\n";}
template <typename T, typename... V>
void _print(T t, V... v) {__print(t); if (sizeof...(v)) cerr << ", "; _print(v...);}
#ifndef ONLINE_JUDGE
#define debug(x...) cerr << "[" << #x << "] = ["; _print(x)
#else
#define debug(x...)
#endif
#pragma endregion DEBUG
#define MOD 1000000007
#define MAX 500005
#define endl "\n"
vector <lli> fact(MAX,1);
vector <lli> invfact(MAX,1);
void roamin(vector<int> &rowmin, vector<vector<int> > &A){
    int n = rowmin.size();
    for(int i = 0; i<rowmin.size(); i++){
        rowmin[i] = *min_element(A[i].begin(),A[i].end());
    }
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            A[i][j] -= rowmin[i];
        }   
    }
}
void colon(vector<int> &colmin, vector<vector<int> > &A){
    int n = colmin.size();
    for(int i = 0; i<n; i++){
        int minm = 1e9;
        for(int j = 0; j<n; j++){
            minm = min(A[j][i],minm);
        }
        colmin[i] = minm;
        for(int j = 0; j<n; j++){
            A[j][i] -= minm;
        }
    }
}
bool check(vector <int> col, vector<int> row, vector<vector<int>> arr){
    int n = arr.size();
    vector <vector<int> > ticked(n,vector<int>(n,0));
    int noZeroes = 0;
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            if(arr[i][j] == 0)
                noZeroes++;
        }
    }
    for(int i = 0; i<row.size(); i++){
        for(int j = 0; j<n; j++){
            ticked[row[i]][j] = true;
        }
    }
    for(int i = 0; i<col.size(); i++){
        for(int j = 0; j<n; j++){
            ticked[j][col[i]] = true;
        }
    }
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            if(ticked[i][j] && arr[i][j] == 0){
                noZeroes--;
            }
        }
    }
    return noZeroes == 0;
    
}
void decompose(vector <int> &a, int beats){
    int idx = 0;
    a.clear();
    while(beats > 0){
        if(beats%2)
            a.push_back(idx);
        idx++;
        beats >>= 1;
    }

}
void prettyPrint(vector<vector<int> > arr){
    int n = arr.size();
    int m = arr[0].size();
    cout<<"____________________________"<<endl;
    for(int i = 0; i<n; i++){
        printf("Worker %d: ", i+1);
        for(int j = 0; j<n; j++){
            cout<<arr[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"____________________________"<<endl;
}
bool check(vector <pair<int,int> > arr){
    for(int i = 0; i<arr.size(); i++){
        for(int j = i+1; j<arr.size(); j++){
            if(arr[i].first == arr[j].first || arr[j].second == arr[i].second)
                return false;
        }
    }
    return true;
}
void hungarian(vector<vector<int> > arr, int n){
    bool flag = 69^420;
    vector <vector<int> > A = arr;
    cout<<"Initial Table: "<<endl;
    prettyPrint(A);
    int minpk = INT_MAX;
    vector <int> rowmin(n);
    vector <int> columnmin(n);
    roamin(rowmin,A);
    colon(columnmin,A);
    cout<<"After subtracting of both row and column minimum:"<<endl;
    prettyPrint(A);
    int iter = 0;
    vector <vector<int> > operator1(69,vector<int>(69,0));
    while(flag){
        if(iter++ > 3)
            return;
        int numzeroes = 0, numsquares = 0;
        vector <int> horline, verline;
        int minm = 500;
        vector <vector<int> > checkmatrix(n,vector<int>(n,0));
        int upper = (1ll<<n);
        for(int i = 0; i<upper; i++){
            vector <int> row;
            decompose(row,i);
            for(int j = 0; j<upper; j++){
                vector <int> col;
                decompose(col,j);
                if(check(col,row,A)){
                    if(row.size() + col.size() < minm){
                        minm = row.size() + col.size();
                        horline = row;
                        verline = col;
                    }
                }
            }
        }
        if(horline.size() + verline.size() == n){
            vector <pair<int,int> > zeroes;
            for(int i = 0; i<n; i++){
                for(int j = 0; j<n; j++){
                    if(A[i][j] == 0)
                        zeroes.push_back({i,j});
                }
            }
            int maxm = (1<<zeroes.size());
            vector <pair<int,int> > finx;
            for(int i = 0; i<maxm; i++){
                vector <int> idx;
                decompose(idx,i);
                if(idx.size() != n)
                    continue;
                finx.clear();
                for(auto i : idx){
                    finx.push_back(zeroes[i]);
                }
                if(check(finx))
                    break;
            }
            int optimalvalue = 0;
            for(auto i : finx){
                optimalvalue += arr[i.first][i.second];
            }
            cout<<"The optimal value is "<<optimalvalue<<endl;
            return;
        }
        for(int i = 0; i<horline.size(); i++){
            for(int j = 0; j<n; j++){
                checkmatrix[horline[i]][j] = true;
            }
        }
        for(int i = 0; i<verline.size(); i++){
            for(int j = 0; j<n; j++){
                checkmatrix[j][verline[i]] = true;
            }
        }
        int minmx = INT_MAX;
        for(int i = 0; i<n; i++){
            for(int j = 0; j<n; j++){
                if(checkmatrix[i][j])
                    continue;
                minmx = min(minmx, A[i][j]);
            }
        }
        for(int i = 0; i<n; i++){
            for(int j = 0; j<n; j++){
                if(checkmatrix[i][j])
                    continue;
                A[i][j] -= minmx;
            }
        }
        for(int i = 0; i<horline.size(); i++){
            for(int j = 0; j<verline.size(); j++){
                A[horline[i]][verline[j]] += minmx;
            }
        }
        cout<<"The checkmatrix after covering all the zeroes is"<<endl;
        for(int i = 0; i<n; i++){
            for(int j = 0; j<n; j++)
                cout<<checkmatrix[i][j]<<" ";
            cout<<endl;
        }
        cout<<endl;
        cout<<"The matrix after uncovered subt covered add"<<endl;
        prettyPrint(A);
    }
}
int main(){
    int n;
    cout<<"Enter the number of workers and tasks"<<endl;
    cin>>n;
    vector <vector<int> > A(n,vector<int>(n));
    cout<<"Enter the values of the table row wise"<<endl;
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            cin>>A[i][j];
        }
    }
    hungarian(A,n);
}