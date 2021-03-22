
import numpy as np

n = int(input("Enter n (no. of variables) : "))
m = int(input("Enter m (no. of equations) : "))
C_j = np.zeros(n+m)
#a_ij = np.empty(n+m,m)
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
B_inv = np.linalg.inv(B)
y = np.dot(C_b,B_inv) # dual of the equation
bv = np.arange(m)
for i in range(m):
    bv[i] = i+n
nbv = np.arange(n)
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
    print('=' + str(  C_j[nbv[i]] - np.dot(y, np.transpose(a_ij[:,nbv[i]] ) ) ))
    if (C_j[nbv[i]] - np.dot(y, np.transpose(a_ij[:,nbv[i]]) )) > max :
        max = C_j[nbv[i]] - np.dot(y, np.transpose(a_ij[:,nbv[i]]) )
        id1 = nbv[i]
        nbcid = i
print('\nThe maximum value of Cj - Zj is :' + str(max) )
print('\nSo the entering variable is ' + str(id1+1) )
check = 1
iterations = 0
while (check == 1):
    check = 0
    print('Now, Pbar_j = B^{-1} Pj = ' + str( np.dot(B_inv,  a_ij[:,id1])) )
    print('So, we divide X_B by Pbar_j to get')
    arr = np.divide(X_B,np.dot(B_inv,  a_ij[:,id1]))
    print(np.divide(X_B,np.dot(B_inv,  a_ij[:,id1])))
    id2 = np.where(arr > 0, arr, np.inf).argmin()
    minim = np.where(arr > 0, arr, np.inf).min()
    print('Here the minimum value of the matrix is ' + str(minim) )
    print('So, the leaving variable is ' + str(bv[id2]+1))
    nbv[nbcid] = bv[id2]
    bv[id2] = id1
    X_B = (X_B) - minim*np.dot(B_inv,  a_ij[:,id1])
    X_B[id2] = minim
    C_b[id2] = C_j[bv[id2]]
    B[:,id2] = a_ij[:,id1]
    B_inv = np.linalg.inv(B)
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
print('The solution is :' )
sumZ = 0
for i in range(m):
    print('x'+str(bv[i]+1) + ' = ' + str(X_B[i]))
    sumZ += C_j[bv[i]] * X_B[i]
print('And the value of Z is : ' + str(sumZ))




