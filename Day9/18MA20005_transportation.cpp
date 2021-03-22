/*
Name: Altaf Ahmad
Roll no: 18MA20005
Transportation Problem - Initial solution using North West Corner method and Least Cost Cell method and optimisation using MOdified DIstribution method
*/
#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define bM 100000000
struct coordinates // defining a structure to get the coordinates of the variables in the table
{
    int x, y;
};
struct threePoint // this structure helps in finding the cycle in MO-DI method
{
    coordinates pos1;
    coordinates neg1, neg2;
};
struct node // this structure helps to solve the equations 
{
    double val;
    bool isInList;
    int index;
};
struct penVal // this stores the penalty values along with the coordinates to change the cycle
{
    double penal;
    coordinates p;
};
bool compareElements(penVal a, penVal b) // this helps to sort the penalty values so that we may be able to process the cycles properly
{

    return a.penal < b.penal;
}
void findVeqn(vector<pair<node, node>> totEquations, int n, double *b, node *x) // this function finds the solution of the (m+n-1) equations which come from the table
{ // what we do here is maintain a list of known variables and unknown variables. We then iterate over the equations to get the value of the unknown variables one by one
    int iter = 0, indx;
    //we know that x[0] = 0, after that, we have to find out x[i], for i = 1,2,..., n
    while (iter++ < n)
    {
        for (int i = 0; i < n; i++)
        {
            x[i].val = 0;
            //if (totEquations[i].first.isInList && totEquations[i].second.isInList)
            //   cout << "x_" << totEquations[i].first.index << " + x_" << totEquations[i].second.index << " = " << totEquations[i].first.val << " + " << totEquations[i].second.val << "\n";
        }
        for (int i = 0; i < n; i++)
        {
            if (totEquations[i].first.isInList) // if one of the values is in the list then the other is got from the equatino
            {
                totEquations[i].second.val = b[i] - totEquations[i].first.val;
                totEquations[i].second.isInList = true; // and then it is put inside the known variable list
                indx = totEquations[i].second.index;
                //x[indx - 1].val = totEquations[i].second.val; 
                for (int j = 0; j < n; j++)
                {
                    if (indx == totEquations[j].first.index)
                    {
                        totEquations[j].first.val = totEquations[i].second.val;
                        totEquations[j].first.isInList = true;
                        //x[indx - 1].val = totEquations[i].second.val;
                    }
                    else if (indx == totEquations[j].second.index)
                    {
                        totEquations[j].second.val = totEquations[i].second.val;
                        totEquations[j].second.isInList = true;
                        //x[indx - 1].val = totEquations[i].second.val;
                    }
                }
            }
            else if (totEquations[i].second.isInList) //similarly for the second variable of the equation
            {
                totEquations[i].first.val = b[i] - totEquations[i].second.val;
                totEquations[i].first.isInList = true;
                indx = totEquations[i].first.index;
                //x[indx - 1].val = totEquations[i].first.val;
                for (int j = 0; j < n; j++)
                {
                    if (indx == totEquations[j].first.index)
                    {
                        totEquations[j].first.val = totEquations[i].first.val;
                        totEquations[j].first.isInList = true;
                        //x[indx - 1].val = totEquations[i].first.val;
                    }
                    else if (indx == totEquations[j].second.index)
                    {
                        totEquations[j].second.val = totEquations[i].first.val;
                        totEquations[j].second.isInList = true;
                        //x[indx - 1].val = totEquations[i].first.val;
                    }
                }
            }
        }
    }
    for (int i = 0; i < n; i++) // after this, we just have to put the values of the variables in the list to be output
    {
        //cout << "x_" << i << " = " << x[i].val << "\n";
        if (totEquations[i].first.isInList && totEquations[i].second.isInList)
        {
            //cout << "x_" << totEquations[i].first.index << " + x_" << totEquations[i].second.index << " = " << totEquations[i].first.val << " + " << totEquations[i].second.val << "\n";
            x[totEquations[i].first.index].val = totEquations[i].first.val;
            x[totEquations[i].second.index].val = totEquations[i].second.val;
        }
    }
}

void printTable(double *a, double *b, double **c, int m, int n, double **x) // helper function to print thte tables one by one
{
    cout << "\n";
    for (int i = 0; i < n; i++)
    {
        cout << "--------------";
    }
    cout << "\n\t \t Destination\n";
    for (int i = 0; i < n; i++)
    {
        cout << "--------------";
    }
    cout << "\n\t";
    for (int i = 0; i < n; i++)
    {
        cout << "\t" << i + 1;
    }
    cout << "\t Supply\n";
    for (int i = 0; i < n; i++)
    {
        cout << "--------------";
    }
    cout << "\n";
    for (int i = 0; i < m; i++)
    {
        cout << "x\t|";
        for (int j = 0; j < n; j++)
        {
            cout << "\t" << x[i][j];
            //printf("\t%0.1lf", x[i][j]);
        }
        cout << "\n";
        cout << i + 1 << "\t|";
        for (int j = 0; j < n; j++)
        {
            cout << "\t" << c[i][j];
        }
        cout << "\t" << a[i] << "\n\n";
    }
    for (int i = 0; i < n; i++)
    {
        cout << "--------------";
    }
    cout << "\n";
    cout << "Demands :";
    for (int i = 0; i < n; i++)
    {
        cout << "\t" << b[i];
    }
    cout << "\n";
}
threePoint circle(double **c, int m, int n, double **x, coordinates theta) // this function helps to find the loop inside the table to update the table in MODi method
{
    // look to the right
    threePoint point; // this will be the returning variable which points to the three points of the edges of the loop
    for (int j = theta.y + 1; j < n; j++) // it looks to the right of the theta cell.
    {
        if (x[theta.x][j] > 0.0) //If there is a cell with value on this side then it looks for cells above nad below it which should have value, that is it is allocated or not
        {
            for (int i = 1; i < (m - theta.x); i++) // this looks for cell above the theta cell
            {
                if ((x[theta.x + i][j] > 0) && (x[theta.x + i][theta.y] > 0)) // and if found it returns it 
                {
                    point.neg1.x = theta.x;
                    point.neg1.y = j;
                    point.neg2.x = theta.x + i;
                    point.neg2.y = theta.y;
                    point.pos1.x = theta.x + i;
                    point.pos1.y = j;
                    return point;
                }
            }
            for (int i = 1; i <= theta.x; i++) // this looks for cells below the theta cell. And if found it returns it's coordinates
            {
                if ((x[theta.x - i][j] > 0) && (x[theta.x - i][theta.y] > 0))
                {
                    point.neg1.x = theta.x;
                    point.neg1.y = j;
                    point.neg2.x = theta.x - i;
                    point.neg2.y = theta.y;
                    point.pos1.x = theta.x - i;
                    point.pos1.y = j;
                    //cout<<" mm";
                    return point;
                }
            }
        }
    }

    // look to the left
    // similar case for left also 
    for (int j = theta.y - 1; j >= 0; j--)
    {
        if (x[theta.x][j] > 0.0)
        {
            for (int i = 1; i < (m - theta.x); i++)
            {
                if (x[theta.x + i][j] > 0 && x[theta.x + i][theta.y] > 0)
                {
                    point.neg1.x = theta.x;
                    point.neg1.y = j;
                    point.neg2.x = theta.x + i;
                    point.neg2.y = theta.y;
                    point.pos1.x = theta.x + i;
                    point.pos1.y = j;
                    return point;
                }
            }
            for (int i = 1; i <= theta.x; i++)
            {
                if (x[theta.x - i][j] > 0 && x[theta.x - i][theta.y] > 0)
                {
                    point.neg1.x = theta.x;
                    point.neg1.y = j;
                    point.neg2.x = theta.x - i;
                    point.neg2.y = theta.y;
                    point.pos1.x = theta.x - i;
                    point.pos1.y = j;
                    return point;
                }
            }
        }
    }
    // look to the top
    // here it first looks on the top and then searches on the right and left
    for (int i = theta.x + 1; i < m; i++)
    {
        if (x[i][theta.y] > 0.0)
        {
            for (int j = 1; j < (n - theta.y); j++)
            {
                if (x[i][theta.y + j] > 0 && x[theta.x][theta.y + j] > 0)
                {
                    point.neg1.x = theta.x;
                    point.neg1.y = theta.y + j;
                    point.neg2.x = i;
                    point.neg2.y = theta.y;
                    point.pos1.x = i;
                    point.pos1.y = theta.y + j;
                    return point;
                }
            }
            for (int j = 1; j <= theta.y; j++)
            {

                if (x[i][theta.y - j] > 0 && x[theta.x][theta.y - j] > 0)
                {
                    point.neg1.x = theta.x;
                    point.neg1.y = theta.y - j;
                    point.neg2.x = i;
                    point.neg2.y = theta.y;
                    point.pos1.x = i;
                    point.pos1.y = theta.y - j;
                    return point;
                }
            }
        }
    }

    // look to the bottom
    //and finally it looks at the bottom first and then looks on the right and left
    for (int i = theta.x - 1; i >= 0; i--)
    {
        if (x[i][theta.y] > 0.0)
        {
            for (int j = 1; j < (n - theta.y); j++)
            {
                if (x[i][theta.y + j] > 0 && x[theta.x][theta.y + j] > 0)
                {
                    point.neg1.x = theta.x;
                    point.neg1.y = theta.y + j;
                    point.neg2.x = i;
                    point.neg2.y = theta.y;
                    point.pos1.x = i;
                    point.pos1.y = theta.y + j;
                    return point;
                }
            }
            for (int j = 1; j <= theta.y; j++)
            {
                if (x[i][theta.y - j] > 0 && x[theta.x][theta.y - j] > 0)
                {
                    point.neg1.x = theta.x;
                    point.neg1.y = theta.y - j;
                    point.neg2.x = i;
                    point.neg2.y = theta.y;
                    point.pos1.x = i;
                    point.pos1.y = theta.y - j;
                    return point;
                }
            }
        }
    }
    point.neg1.x = -1; // finally if it doesn't find a loop it returns this
    return point;
}
int MODIoptim(double *a, double *b, double **c, int m, int n, double **x) // Finds the optimum solution after the basic solution is found out by the North west corner or least cost cell method depending on the problem
{
    double *u = new double[m];
    double *v = new double[n];
    int kkk = 0;
    int count = 0;
    int pq = 0;
    double extravar;
    int **xx = new int *[100];
    for (int i = 0; i < 100; i++)
    {
        xx[i] = new int[100];
    }
    bool isSelect[100][100];

    double **x3 = new double *[100];
    for (int i = 0; i < 100; i++)
    {
        x3[i] = new double[100];
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            x3[i][j] = x[i][j];
            isSelect[i][j] = false;
            xx[i][j] = (x[i][j] + 0.5);
            if (x[i][j] > 0.0)
            {
                isSelect[i][j] = true;
                count++;
            }
        }
    }
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << "\t" << isSelect[i][j];
        }
        cout << "\n";
    }

    double pen;
    double maxpen;
    coordinates p;
    coordinates q;
    cout << "Total number of select variables are " << count << "\n";
    q.x = -1;
    q.y = -1;
    if (count != (m + n - 1))
    {
        cout << "MO DI method cannot be done on this \n";
        return 0;
    }
    double *b2 = new double[100];
    node *xstar = new node[100];
    threePoint pts;
    int count1;
    int t = 1;
    penVal *penalty = new penVal[200];
    //________________________________________________
    int iter = 0;
    while (1)
    {

        pair<node, node> eqn;
        vector<pair<node, node>> totEquations; // this vector stores the equations that need to be solved to find out u_0, u_1, ..., u_m , v_0, v_1, ...,v_n
        count1 = 0;
        cout << "Select variables are \n"; // this table shows the allocated variables
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cout << "\t" << isSelect[i][j];
            }
            cout << "\n";
        }
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (isSelect[i][j] == 1)
                {
                    cout << "u_" << i << " + v_" << j << " = " << c[i][j] << "\n";
                    eqn.first.index = i;
                    if (i == 0)
                    {
                        eqn.first.isInList = true; // we are taking the value of u_0 = 0. and thus solving using this.
                        eqn.first.val = 0.0;
                    }
                    else
                    {
                        eqn.first.isInList = false;
                    }
                    eqn.second.index = m + j;

                    eqn.second.isInList = false; // rest all are unknown, so it is false here

                    totEquations.push_back(eqn);
                    b2[count1] = c[i][j];
                    count1++;
                }
            }
        }
        //printTable(a,b,c,m,n,x);
        if (count1 != (m + n - 1)) // degeneracy case. Can be resolved if found later
        {
            cout << "MO DI method cannot be done on this \n";
            return 0;
        }

        findVeqn(totEquations, count, b2, xstar); // finds the solution of the equations
        u[0] = 0;
        cout << "m = " << m << " , n = " << n << " m + n - 1 = " << count << "\n";
        cout << "Taking u_0 as 0, we get\n";
        for (int i = 1; i <= count; i++)
        {
            if (i < (m))
            {
                u[i] = xstar[i].val;
                cout << "u_" << i << " = " << u[i] << "\n";
            }
            else
            {
                v[i - m] = xstar[i].val;
                cout << "v_" << i - m << " = " << v[i - m] << "\n";
            }
        }

        t = 1;
        maxpen = -3000;
        t = 1;
        if (iter++ > 10)
            return 0;
        kkk = 0;

        for (int i = 0; i < m; i++) // now we check for the penalty values and find the maximum penalty
        {
            for (int j = 0; j < n; j++)
            {
                pen = u[i] + v[j] - c[i][j];
                if (isSelect[i][j] == 0)
                {
                    penalty[kkk].penal = pen;
                    penalty[kkk].p.x = i;
                    penalty[kkk].p.y = j;
                    //cout<<i <<" = i, j = "<<j <<" . k = "<<kkk<<"\n";
                    kkk++;
                }
                cout << "penalty[i][j] = u[" << i << "] + v[" << j << "] - c [" << i << "][" << j << "] = " << pen << "\n";
                if (pen > 0.0)
                {
                    t = 0;
                    if (pen > maxpen)
                    {
                        maxpen = pen;
                        p.x = i;
                        p.y = j;
                    }
                }
            }
        }

        if (t || (iter > 9)) // if the penalty values are less then equal to 0 then the optimum solution is found 
        {
        label:
            cout << "Optimality has reached\n";
            printTable(a, b, c, m, n, x3);
            double sum = 0.0;
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    sum += c[i][j] * x3[i][j];
                }
            }
            cout << "The cost from the optimal solution of MODI method is : ";
            printf("%0.1lf\n", sum);
            return 0;
        }
        // other wise we find the maximum penalty element and find the cycles
        sort(penalty, penalty + ((m * n) - count), compareElements);
        for (int i = 0; i < ((m * n) - count); i++)
        {
            cout << "Penalty at (" << penalty[i].p.x << " , " << penalty[i].p.y << ") = " << penalty[i].penal << " ";
        }
        cout << "\n";
        pq = 0;

        while (1)
        {
            pq++;
            p = penalty[(m * n) - count - pq].p; // this is the penalty cell

            cout << "The point p lies on the point " << p.x << " , " << p.y << " with values = " << c[p.x][p.y] << " and penalty =  " << penalty[(m * n) - count - pq].penal << "\n";
            pts = circle(c, m, n, x3, p); // finds the three points in the loop
            if (penalty[(m * n) - count - pq].penal < 0) 
            {
                cout << " pen = " << penalty[(m * n) - count - pq].penal << "\n This may not be the optimal solution but it's a noice solution\n";
                goto label;
            }
            if (pts.neg1.x == -1) // if no loop is found it looks for the next positive penalty to minimize
            {
                continue;
            }
            else // otherwise the table is updated
            {
                if (penalty[(m * n) - count - pq].penal < 0)
                {
                    cout << " pen = " << penalty[(m * n) - count - pq].penal << "\n";
                    goto label;
                }
                extravar = min(x3[pts.neg2.x][pts.neg2.y], x3[pts.neg1.x][pts.neg1.y]);
                cout << "The removing variable's value is " << extravar << "\n";
                if (abs(extravar) < 0.01)
                    continue;
                break;
            }
        }
        //cout<<"The negative valued sides are "<<x[pts.neg1.x][pts.neg1.y] <<" and "<< x[pts.neg2.x][pts.neg2.y] <<"\n";
        double k = min(x3[pts.neg2.x][pts.neg2.y], x3[pts.neg1.x][pts.neg1.y]); // here the table values are updated
        x3[pts.neg2.x][pts.neg2.y] -= k;
        x3[pts.neg1.x][pts.neg1.y] -= k;
        x3[p.x][p.y] += k;
        isSelect[p.x][p.y] = true;
        x3[pts.pos1.x][pts.pos1.y] += k;
        isSelect[pts.pos1.x][pts.pos1.y] = true;
        if (x3[pts.neg2.x][pts.neg2.y] > x3[pts.neg1.x][pts.neg1.y])
        {
            isSelect[pts.neg1.x][pts.neg1.y] = false;
        }
        else
        {
            isSelect[pts.neg2.x][pts.neg2.y] = false;
        }

        xx[pts.neg2.x][pts.neg2.y] = (x3[pts.neg2.x][pts.neg2.y] + 0.5);
        xx[pts.neg1.x][pts.neg1.y] = (x3[pts.neg1.x][pts.neg1.y] + 0.5);
        xx[p.x][p.y] = (x3[p.x][p.y] + 0.5);
        xx[pts.pos1.x][pts.pos1.y] = (x3[p.x][p.y] + 0.5);

        double sum = 0.0;
        for (int i = 0; i < m; i++) // finds the cost of the final solution
        {
            for (int j = 0; j < n; j++)
            {
                sum += c[i][j] * x3[i][j];
            }
        }
        cout << "The cost from the solution now is: ";// 
        printf("%0.1lf\n", sum);
        cout << "Updating table we get";
        printTable(a, b, c, m, n, x3);
    }
    //_____________________________________
    return 0;
}
void northWestCorner(double *a, double *b, double **c, int m, int n, double **x) // finds the bfs using north west corner method
{
    double *aa = new double[m];
    double *bb = new double[n];
    for (int i = 0; i < m; i++)
    {
        aa[i] = a[i];
    }
    for (int i = 0; i < n; i++)
    {
        bb[i] = b[i];
    }

    int col = 0, row = 0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            x[i][j] = 0; // initialises decision variables as 0 
        }
    }
    int iter = 0;
    while ((col != (n)) && (row != (m))) // starts from the north western cells and gradually moves down cells one by one
    {
        if (iter++ > (m + n))
        {
            break;
        }
        if (aa[row] > bb[col]) // if the supply value is larger, then it does this
        {
            x[row][col] = bb[col];
            aa[row] -= bb[col];
            bb[col] = 0.0;
            col++;
            printTable(aa, bb, c, m, n, x);
        }
        else // other wise this is done
        {
            x[row][col] = aa[row];
            bb[col] -= aa[row];
            aa[row] = 0.0;
            row++;
            printTable(aa, bb, c, m, n, x);
        }
    }
    double sum = 0.0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            sum += c[i][j] * x[i][j];
        }
    }
    cout << "The cost from the basic feasible solution of North West Corner method is : " << sum << "\n";
}
coordinates minInC(double **c, int m, int n) // helper function to find the minimum c in the lCM method
{
    double absolute = bM;
    coordinates min;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (c[i][j] < absolute)
            {
                min.x = i;
                min.y = j;
                absolute = c[i][j];
            }
        }
    }
    return min;
}
void LeastCostCellMethod(double *a, double *b, double **c, int m, int n, double **x) // finds the BFS using LCM method
{
    double *aa = new double[m];
    double *bb = new double[n];
    double **cc = new double *[m];
    int row, col;
    for (int i = 0; i < m; i++)
    {
        cc[i] = new double[n];
    }
    double delivery = 0;
    for (int i = 0; i < m; i++)
    {
        aa[i] = a[i];
    }
    for (int i = 0; i < n; i++)
    {
        bb[i] = b[i];
        delivery += b[i]; // here the delivery stores the vales of the total shipment left
    }
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            x[i][j] = 0; // initialises the decision variables as 0
            cc[i][j] = c[i][j]; // cost is kept in another array so that the final value is retained
        }
    }20
    coordinates min;
    int iter = 0;
    while (delivery > 0.0001) // loops untill all the deliveries are done
    {
        if (iter++ > 10)
            break;
        min = minInC(cc, m, n); // finds the minimum cell
        row = min.x;
        col = min.y;
        if (aa[min.x] > bb[min.y]) // and then the cell calues are changed accordingly
        {
            for (int i = 0; i < m; i++)
            {
                cc[i][min.y] = bM;
            }
            x[row][col] = bb[col];
            aa[row] -= bb[col];
            bb[col] = 0.0;
            printTable(aa, bb, c, m, n, x);
        }
        else
        {
            for (int j = 0; j < n; j++)
            {
                cc[min.x][j] = bM;
            }
            x[row][col] = aa[row];
            bb[col] -= aa[row];
            aa[row] = 0.0;
            row++;
            printTable(aa, bb, c, m, n, x);
        }

        delivery = 0.0;
        for (int i = 0; i < n; i++)
        {
            delivery += bb[i];
        }
        cout << "Delivery = " << delivery << "\n"; 
    }

    double sum = 0.0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            sum += c[i][j] * x[i][j];
        }
    }
    cout << "The cost from the basic feasible solution of Least Cost Cell method is : " << sum << "\n";
}
int main()
{
    int n, m;
    cout << "Enter m (The number of sources): ";
    cin >> m;
    cout << "Enter n (The number of destinations): ";
    cin >> n;
    double *a = new double[m]; // stores the supply of the sources
    double *b = new double[n]; // stores the demands of the destinations
    double **x = new double *[100]; // stores the decision variables
    for (int i = 0; i < 100; i++) 
    { 
        x[i] = new double[100];
    }
    double asum = 0.0, bsum = 0.0;
    cout << "Now, enter the a_i's. That is the supply at the ith source : \n";
    for (int i = 0; i < m; i++)
    {
        cout << "Supply at source " << i + 1 << " :";
        cin >> a[i];
        asum += a[i];
    }
    cout << "Now, enter the b_i's. That is the demand at the ith destination : \n";
    for (int i = 0; i < n; i++)
    {
        cout << "Demand at destination " << i + 1 << " :";
        cin >> b[i];
        bsum += b[i];
    }
    if (abs(bsum - asum) < 0.001)
    {
        cout << "Transportation problem is balanced\n";
    }
    else
    {
        cout << "The problem is not balanced\n";
        return 0;
    }
    double **c = new double *[m]; // stores the cost of transportation from source i to destination j as c[i][j]
    for (int i = 0; i < m; i++)
    {
        c[i] = new double[n];
    }
    cout << "Now, we have to enter cost of commodity i to be shipped to the destination j from the source i as c[i][j]\n";
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("Enter c[%d][%d] :", i, j);
            cin >> c[i][j];
            x[i][j] = 0;
        }
    }

    printTable(a, b, c, m, n, x);
    northWestCorner(a, b, c, m, n, x);
    printTable(a, b, c, m, n, x);
    LeastCostCellMethod(a, b, c, m, n, x);
    cout << "Phase 1 is complete and now we have to move to phase II. This is the table we got after phase I \n";
    printTable(a, b, c, m, n, x);
    MODIoptim(a, b, c, m, n, x);
}