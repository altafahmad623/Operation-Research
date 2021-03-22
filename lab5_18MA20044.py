import numpy as np
MAX_ITERS = 1000
BIG_M = 1000000000

def two_phase_solve(m, n, y, c_B, coeffs_obj, row_vars, lpp_type, extra):
    print ("STARTING FIRST PHASE")
    k = 0
    new_coeffs_obj = [0.0]*m
    for elem in extra:
        if elem == 3:
            new_coeffs_obj.append(-1.0)
        else:
            new_coeffs_obj.append(0.0)
    
    for i in range(n):
        c_B[i] = new_coeffs_obj[row_vars[i]]
    print(c_B)

    while (k<1000):
        k += 1
        print("################################################")
        print("Iteration:", k)
        print ("y is: ")
        print(y)
        print("z_j for all j is:")
        z = np.matmul(c_B, y)
        print(z)

        print("z_j - c_j for all j is")
        diff = (z[:z.shape[0]-1] - new_coeffs_obj)
        print (diff)

        if (np.amin(diff)>=0):
            break
        colselect = (np.where(diff == np.amin(diff)))[0][0]
        # print([int(position[0][0])])
        ratio = y[:, y.shape[1]-1] / y[:, colselect]
        print("Ratios are:")
        print(ratio)
        for i,elem in enumerate(ratio):
            if elem < 0:
                ratio[i] = np.Infinity
        print("min +ve ratio is:")
        print(np.amin(ratio))

        rowselect = (np.where(ratio == np.amin(ratio)))[0][0]

        print ("pivot element position is:")
        print ((rowselect+1, colselect+1))
        print ("pivot element is:")
        pivot = y[rowselect, colselect]
        print(pivot)

        for i in range(y.shape[0]):
            if i != rowselect:
                y[i,:] = y[i,:] - y[rowselect,:] * y[i,colselect] / pivot

        y[rowselect,:] = y[rowselect,:] / pivot
        c_B[rowselect] = new_coeffs_obj[colselect]
        row_vars[rowselect] = colselect
        print ("variables used are: ")
        print(row_vars)
    
    if k == MAX_ITERS:  
        print("unbounded solution, terminating")

    else:
        print("optimal solution is:")
        solution = []
        for i in range(y.shape[1]-1):
            if i not in row_vars:
                solution.append(0)
            else:
                solution.append(y[row_vars.index(i), y.shape[1]-1])

        feasible = True
        for i in range(m, len(solution)):
            if extra[i-m] == 3 and solution[i] != 0:
                feasible = False
        
        if feasible == False:
            print("SOLUTION IS INFEASIBLE")
        else:
            optimum = 0.0
            for i in range(m):
                print ("value of   x", i+1, "  = ", solution[i])
                optimum += solution[i]*new_coeffs_obj[i]

            print("optimal value of objective function is")
            print(optimum * (0.5 - lpp_type) * 2)
    
    print ("STARTING SECOND PHASE")
    new_coeffs_obj = coeffs_obj
    for i in range(n):
        c_B[i] = new_coeffs_obj[row_vars[i]]
    print(c_B)
    k=0
    while (k<1000):
        k += 1
        print("################################################")
        print("Iteration:", k)
        print ("y is: ")
        print(y)
        print("z_j for all j is:")
        z = np.matmul(c_B, y)
        print(z)

        print("z_j - c_j for all j is")
        diff = (z[:z.shape[0]-1] - new_coeffs_obj)
        print (diff)

        if (np.amin(diff)>=0):
            break
        colselect = (np.where(diff == np.amin(diff)))[0][0]
        # print([int(position[0][0])])
        ratio = y[:, y.shape[1]-1] / y[:, colselect]
        print("Ratios are:")
        print(ratio)
        for i,elem in enumerate(ratio):
            if elem < 0:
                ratio[i] = np.Infinity
        print("min +ve ratio is:")
        print(np.amin(ratio))

        rowselect = (np.where(ratio == np.amin(ratio)))[0][0]

        print ("pivot element position is:")
        print ((rowselect+1, colselect+1))
        print ("pivot element is:")
        pivot = y[rowselect, colselect]
        print(pivot)

        for i in range(y.shape[0]):
            if i != rowselect:
                y[i,:] = y[i,:] - y[rowselect,:] * y[i,colselect] / pivot

        y[rowselect,:] = y[rowselect,:] / pivot
        c_B[rowselect] = new_coeffs_obj[colselect]
        row_vars[rowselect] = colselect
        print ("variables used are: ")
        print(row_vars)
    
    if k == MAX_ITERS:  
        print("unbounded solution, terminating")

    else:
        print("optimal solution is:")
        solution = []
        for i in range(y.shape[1]-1):
            if i not in row_vars:
                solution.append(0)
            else:
                solution.append(y[row_vars.index(i), y.shape[1]-1])

        feasible = True
        for i in range(m, len(solution)):
            if extra[i-m] == 3 and solution[i] != 0:
                feasible = False
        
        if feasible == False:
            print("SOLUTION IS INFEASIBLE")
        else:
            optimum = 0.0
            for i in range(m):
                print ("value of   x", i+1, "  = ", solution[i])
                optimum += solution[i]*new_coeffs_obj[i]

            print("optimal value of objective function is")
            print(optimum * (0.5 - lpp_type) * 2)

#######################################################
# Q1
print("______________________________________________________________________________")
print("Q1")
y = np.array(   [[ 1,  4,  2, -1,  1,  0,  0,  5],
                [ 3,  1,  2,  0,  0, -1,  1,  4]], dtype=float
            )
c_B = np.array([-1, -1])
new_coeffs_obj = [-2.0, -9.0, -1.0,  0.0, -1,  0.0, -1]
row_vars = [4, 6]
extra = [1, 3, 1, 3]
two_phase_solve(3, 2, y, c_B, new_coeffs_obj, row_vars, 1, extra)
#######################################################
#######################################################
# Q2
print("______________________________________________________________________________")
print("Q2")
y = np.array(   [[ 3,  1, -1,  1,  0,  0,  0,  0, 27],
                [ 1,  1,  0,  0, -1,  1,  0,  0, 21],
                [ 1,  2,  0,  0,  0,  0, -1,  1, 30]], dtype=float
            )
c_B = np.array([-BIG_M, -BIG_M, -BIG_M])
new_coeffs_obj = [-4.0, -2.0, 0.0, -BIG_M,  0.0, -BIG_M, 0.0, -BIG_M]
row_vars = [3, 5, 7]
extra = [1, 3, 1, 3, 1, 3]
two_phase_solve(2, 3, y, c_B, new_coeffs_obj, row_vars, 1, extra)

# #######################################################
# #######################################################
# # Q3
print("______________________________________________________________________________")
print("Q3")
y = np.array(   [[ 3,  1,  1,  0,  0,  0,  3],
                [ 4,  3,  0, -1,  1,  0,  6],
                [ 1,  2,  0,  0,  0,  1,  4]], dtype=float
            )
c_B = [-BIG_M, -BIG_M, 0]
new_coeffs_obj = [-2.0, -1.0, -1.0, 0.0, -BIG_M, 0.0]
row_vars = [2, 4, 5]
extra = [3, 1, 3, 2]
two_phase_solve(2, 3, y, c_B, new_coeffs_obj, row_vars, 0, extra)

#######################################################
#######################################################
# Q4
print("______________________________________________________________________________")
print("Q4")
y = np.array(   [[ 2,  1, -1,  1,  0,  0,  2],
                [ 1,  3,  0,  0,  1,  0,  3],
                [ 0,  1,  0,  0,  0,  1,  4]], dtype=float
            )
c_B = [-BIG_M, 0, 0]
new_coeffs_obj = [3, -1, 0, -BIG_M, 0, 0]
row_vars = [3, 4, 5]
extra = [1, 3, 2, 2]
two_phase_solve(2, 3, y, c_B, new_coeffs_obj, row_vars, 0, extra)

#######################################################
#######################################################
# Q5
print("______________________________________________________________________________")
print("Q5")
y = np.array(   [[ 1,  2,  3,  0,  1,  0,  0, 15],
                [ 2,  1,  5,  0,  0,  1,  0, 20],
                [ 1,  2,  1,  1,  0,  0,  1, 10]], dtype=float
            )
c_B = [-BIG_M, -BIG_M, -BIG_M]
new_coeffs_obj = [1, 2, 3, -1, -BIG_M, -BIG_M, -BIG_M]
row_vars = [4, 5, 6]
extra = [3, 3, 3]
two_phase_solve(4, 3, y, c_B, new_coeffs_obj, row_vars, 0, extra)

#######################################################
#######################################################
# Q6
print("______________________________________________________________________________")
print("Q6")
y = np.array(   [[ 1, -2,  3,  1,  0,  2],
                [ 3,  2,  4,  0,  1,  1]], dtype=float
            )
c_B = [-BIG_M, -BIG_M]
new_coeffs_obj = [2, 1, 3, -BIG_M, -BIG_M]
row_vars = [3, 4]
extra = [3, 3]
two_phase_solve(3, 2, y, c_B, new_coeffs_obj, row_vars, 0, extra)

#######################################################

print("______________________________________________________________________________")
print ("Enter number of equations: ")
n = int(input())

coeff_matrix = []
type_constraint = []
type_variable = []
b_array = []

for i in range(n):
    print("Enter coefficients of inequality/equality a_i0, a_i1, a_i2, ... a_im  where  a_i0 x_0 + a_i1 x_1 + ... a_im x_m = b_i : ")
    coeff_array = [ float(x) for x in input().split()]
    coeff_matrix.append(coeff_array)
    print("Enter b_i  where  a_i0 x_0 + a_i1 x_1 + ... a_im x_m = b_i : ")
    b_elem = float(input())
    b_array.append(b_elem)
    print("Enter the type of constraint:")
    print("\n 1 = greater than or equal to \n 0 = equality \n 2 = less than or equal to")
    constraint = int(input())
    type_constraint.append(constraint)


coeffs_obj = []
print("Enter the coeffs for objective function:")
coeffs_obj = [ float(x) for x in input().split()]

print("Enter maximisation (0) or minimisation (1):")
lpp_type = int(input())

new_coeffs_obj = []
for coeff in coeffs_obj:
    if lpp_type == 1:
        new_coeffs_obj.append(-coeff)
    else:
        new_coeffs_obj.append(coeff)

m = len(coeff_array)
for i in range(m):
    type_variable.append(-1)

extra = []
belongsto = []

for i, coeff in enumerate(coeff_matrix):
    if (type_constraint[i] == 0):
        extra.append(3)
        belongsto.append(i)
    if (type_constraint[i] == 1):
        extra.append(1)
        extra.append(3)
        belongsto.append(i)
        belongsto.append(i)
    if (type_constraint[i] == 2):
        extra.append(2)
        belongsto.append(i)

type_constraint = type_constraint + extra 

new_coeff_matrix = []

for j, coeffs in enumerate(coeff_matrix):
    for i, elem in enumerate(extra):
        if j == belongsto[i]:
            if elem == 3:
                coeffs.append(1)
            elif elem == 2:
                coeffs.append(1.0)
            elif elem == 1:
                coeffs.append(-1.0)
        else:
            coeffs.append(0)
    new_coeff_matrix.append(coeffs)

for i, elem in enumerate(extra):
    if elem == 3:
        new_coeffs_obj.append(-BIG_M)
    else:
        new_coeffs_obj.append(0)

A_matrix = np.array(new_coeff_matrix)

print ("A matrix is:")
print (A_matrix)
print ("c vector is:")
print (np.vstack(np.array(new_coeffs_obj)))
print ("b vector is:")
b_vector = np.vstack(np.array(b_array))
print (b_vector)

print("The A matrix is:")
print (A_matrix)

y = np.zeros((A_matrix.shape[0], A_matrix.shape[1] + 1))
y[:,:A_matrix.shape[1]] = A_matrix
y[:, A_matrix.shape[1]:A_matrix.shape[1]+1] = b_vector

row_vars = [0]*n
for i in range(n):
    artificial = -1
    slacksurplus = -1
    for j, elem in enumerate(extra):
        if i == belongsto[j]:
            if elem == 3:
                artificial = j
            else:
                slacksurplus = j
    if artificial > -1:
        row_vars[i] = m + artificial
    else:
        row_vars[i] = m + slacksurplus

print("coefficients of Basic variables is: ")
c_B = np.array([0]*n)
for i in range(n):
    c_B[i] = new_coeffs_obj[row_vars[i]]
print(c_B)

two_phase_solve(m, n, y, c_B, new_coeffs_obj, row_vars, lpp_type, extra)

