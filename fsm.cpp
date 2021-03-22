#include <bits/stdc++.h>
#include <string>
using namespace std;
int delta(char c, int init_state)
{
    if(c == 'b')
    {
        return init_state;
    }
    else if(c == 'a')
    {
        if(init_state == 3)
        {
            return init_state;
        }
        else
        {
            return (1 + init_state);
        }
        
    }
    return 0;
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
    if(state == 3)
    {
        cout<<"Accepted"<<endl;
    }
    else
    {
        cout<<"Not Accepted"<<endl;
    }
    
}