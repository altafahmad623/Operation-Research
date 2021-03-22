#include <bits/stdc++.h>
#include <string>
using namespace std;
int delta(char c, int init_state)
{
    if(init_state == 0)
    {
        if(c == '0')
        {
            return init_state;
        }
        else
        {
            return (init_state+1);
        }
    }
    else if (init_state == 1 )
    {
        if(c == '0')
        {
            return (1+init_state);
        }
        else
        {
            return (init_state-1);
        }
    }
    else
    {
        if(c == '0')
        {
            return (init_state-1);
        }
        else
        {
            return (init_state);
        }
    }
    
}
int main()
{
    int state = 0;
    string input;
    getline(cin,input);
    int n = input.length();
    for (int i = 0; i < n; i++)
    {
        state = delta(input[i],state);
    }
    if(state == 0)
    {
        cout<<"Accepted"<<endl;
    }
    else
    {
        cout<<"Not Accepted"<<endl;
    }
    
}