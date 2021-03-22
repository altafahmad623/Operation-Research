"""
Implementation of Revised Simplex Method

Name : Altaf Ahmad
Roll no: 18MA20005

"""

import numpy as np
from tabulate import tabulate # this is a library to print the simplex tables. If not installed, it can be easily installed using " pip install tabulate"

def inverse(A): # here we define the function to calculate the inverse which will be used to find B^{-1}
    row = len(A) # order of the square matrix,. i.e. =n if it is a nxn matrix
    
    M = [] # this will find out the adjacency matrix and then we will find the inverse using this
    for i in range(row):
        M.append([]) # adding one row in each line
        for j in range(row):
            M[i].append(A[i][j]) # first we are adding each element as itself
        for j in range(row):
            if i==j:
                M[i].append(1) # after that we are adding the Identity matrix
            else:
                M[i].append(0)
                                
    for i in range(row-1,0,-1): # here we traverse it down one by one and exchange the cells, thta is creating the transpose
        if M[i-1][0] < M[i][0]:
            for j in range(2*row):
                temp = M[i][j]
                M[i][j] = M[i-1][j]
                M[i-1][j] = temp
                
    for i in range(row):
        for j in range(row):
            if j != i:
                if(M[i][i] == 0):
                    for k in range(i+1,row):
                        if(M[k][i]!=0):
                            for l in range(2*row):
                                M[i][l]+=M[k][l]
                            break
                temp = M[j][i] / M[i][i]
                for k in range(2*row):
                    M[j][k] -= M[i][k] * temp
           
    for i in range(row): # now we will divide the elements by the determinant to get to the inverse
        temp = M[i][i]
        for j in range(2*row):
            M[i][j] = M[i][j] / temp
                        
    I = [] # finally, this will store the inverse
    for i in range(row):
        I.append([])
        for j in range(row):
            I[i].append(M[i][j+row])
            
    return I
def multiply(A,B): # function to multiply two matrices
    
    C = []

    for i in range(len(A)):
        C.append([])
        for j in range(len(B[0])):
            C[i].append(0)
            for k in range(len(B)):
                C[i][j] += A[i][k] * B[k][j]

    return C
n = int(input("Enter n (no. of variables) : "))
m = int(input("Enter m (no. of equations) : "))
C_j = np.zeros(n+m)
#a_ij = np.empty(n+m,m)
header = ['CB', 'BV'] # this is a helper list to print the table 
for i in range(n+m):
    stri = 'x_'+str(i+1)
    header.append(stri)
header.append('b')
for i in range(n):
    print('Enter the coefficient of x'+str(i+1)+ ' in Z :')
    x = float(input())
    C_j[i] = x
a_ij = np.zeros((m,m+n))
for i in range(m):
    for j in range(n):
        print('Enter the coefficient of x'+str(j+1)+' in equation '+str(i+1))
        x = float(input())
        a_ij[i,j] = x
    a_ij[i,(n+i)] = 1
X_B = np.zeros(m)
for i in range(m):
    print('Enter the value of b' +str(i+1) + ' in AX= b')
    x = float(input())
    X_B[i] = x
B = np.eye(m) # basis
C_b = np.zeros(m) # coefficients of the basic variables
B_inv = inverse(B) # calculates the inverse of the Basis matrix using our function
y = np.dot(C_b,B_inv) # dual of the equation
bv = np.arange(m) # adds and stores the basic variables
for i in range(m):
    bv[i] = i+n
nbv = np.arange(n) # stores the non basic variables
print('The initial table is :')
rows, cols = (m, m+n + 3) # prints the initial table as asked in part (a)
arr = [[0 for i in range(cols)] for j in range(rows)] 
for i in range(rows):
    arr[i][1] = bv[i] + 1
    arr[i][0] = C_j[bv[i]]
for i in range(m):
    for j in range(m+n):
        arr[i][j+2] = a_ij[i][j]
for i in range(m):
    arr[i][cols-1] = X_B[i]
table = tabulate(arr, headers =header, tablefmt ='orgtbl')
print(table)
print('\n')
print('The basic variables are : ' + str(bv+1) ) # prints the basic and non basic variables as asked in (b) for the ith iteration
print('With the coefficients C_b : ' + str(C_b))
print('And X_B = ' + str(X_B) )
print('The non-basic variables are : '+ str(nbv + 1))
print('B : \n'+str(B))
print('B^{-1} : \n' + str(B_inv))
print('y = C_b.B^{-1} = ' + str(y))
id1 =0
nbcid = 0
max = 0
for i in range(n): # loop to calculate the values of C_j - Z_j and also maximum C_j - Z_j
    print('P' + str(nbv[i] + 1) + ' = '+ str(a_ij[:,nbv[i]]))
    print('Now, C' + str(nbv[i] + 1) + ' - Z' + str(nbv[i] + 1) + ' = C' + str(nbv[i] + 1) + ' - yP' + str(nbv[i] + 1) )
    print('=' + str(  C_j[nbv[i]] - np.dot(y, np.transpose(a_ij[:,nbv[i]] ) ) ))
    if (C_j[nbv[i]] - np.dot(y, np.transpose(a_ij[:,nbv[i]]) )) > max :
        max = C_j[nbv[i]] - np.dot(y, np.transpose(a_ij[:,nbv[i]]) )
        id1 = nbv[i]
        nbcid = i
print('\nThe maximum value of Cj - Zj is :' + str(max) )
print('\nSo the entering variable is ' + str(id1+1) )
check = 1
iterations = 0
while (check == 1): # loop for the iterations of the revised simplex technique
    check = 0
    iterations +=1
    print('\n Iteration : '+str(iterations))
    rows, cols = (m, m+n + 3) 
    arr = [[0 for i in range(cols)] for j in range(rows)]  # prints the tables here
    for i in range(rows):
        arr[i][1] = bv[i] + 1
        arr[i][0] = C_j[bv[i]]
    for i in range(m):
        for j in range(m+n):
            arr[i][j+2] = a_ij[i][j]
    for i in range(m):
        arr[i][cols-1] = X_B[i]
    table = tabulate(arr, headers =header, tablefmt ='orgtbl')
    print(table)
    print('\n')
    print('Now, Pbar_j = B^{-1} Pj = ' + str( np.dot(B_inv,  a_ij[:,id1])) )
    print('So, we divide X_B by Pbar_j to get')
    arr = np.divide(X_B,np.dot(B_inv,  a_ij[:,id1])) # now we find Pbar here to find the leaving variable
    print(np.divide(X_B,np.dot(B_inv,  a_ij[:,id1])))
    id2 = np.where(arr > 0, arr, np.inf).argmin()
    minim = np.where(arr > 0, arr, np.inf).min()
    print('Here the minimum value of the matrix is ' + str(minim) )
    print('So, the leaving variable is ' + str(bv[id2]+1))
    nbv[nbcid] = bv[id2] # enters the entering variable into the
    bv[id2] = id1 # 					list of basic variables			
    X_B = (X_B) - minim*np.dot(B_inv,  a_ij[:,id1])
    X_B[id2] = minim
    C_b[id2] = C_j[bv[id2]]
    B[:,id2] = a_ij[:,id1] # we now follow the same steps and then
    B_inv = inverse(B)     # go on wiht the iterations
    y = np.dot(C_b,B_inv)  	
    print('The basic variables are : ' + str(bv+1) )
    print('With the coefficients C_b : ' + str(C_b))
    print('And X_B = ' + str(X_B) )
    print('The non-basic variables are : '+ str(nbv + 1))
    print('B : \n'+str(B))
    print('B^{-1} : \n' + str(B_inv))
    print('y = C_b.B^{-1} = ' + str(y))
    id1 =0
    nbcid = 0
    max = 0
    for i in range(n):
        print('P' + str(nbv[i] + 1) + ' = '+ str(a_ij[:,nbv[i]]))
        print('Now, C' + str(nbv[i] + 1) + ' - Z' + str(nbv[i] + 1) + ' = C' + str(nbv[i] + 1) + ' - yP' + str(nbv[i] + 1) )
        if (C_j[nbv[i]] - np.dot(y, np.transpose(a_ij[:,nbv[i]] ))) > 0:
            check = 1
        print('=' + str(  C_j[nbv[i]] - np.dot(y, np.transpose(a_ij[:,nbv[i]] ) ) ))
        if (C_j[nbv[i]] - np.dot(y, np.transpose(a_ij[:,nbv[i]]) )) > max :
            max = C_j[nbv[i]] - np.dot(y, np.transpose(a_ij[:,nbv[i]]) )
            id1 = nbv[i]
            nbcid = i
    if check == 0:
        break
    print('\nThe maximum value of Cj - Zj is :' + str(max) )
    print('\nSo the entering variable is ' + str(id1+1))
    if iterations > 10:
        break
if iterations > 10:
    print('It is unbounded')
else: 
    print('The solution is :' )
    sumZ = 0
    for i in range(m):
        print('x'+str(bv[i]+1) + ' = ' + str(X_B[i]))
        sumZ += C_j[bv[i]] * X_B[i]
    print('And the rest of the values are 0.')
    print('And the value of Z is : ' + str(sumZ))
