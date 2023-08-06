import numpy as np

u = np.array([1, 2, 3])
v = np.array([4, 5, 6])
w = np.array([u[0]+v[0], u[1]+v[1], u[2]+v[2]]) # i.e. w = [5,7,9]

print( u+v == w )

u = np.array([1, 2, 3])
v = np.array([2*u[0], 2*u[1], 2*u[2]]) # i.e. v = [2,4,6]

print( v == 2*u )

A = np.array( [ [1,2,3], [4,5,6] ] )

print(A[1][1])
print(A[1,1])
print(A[0][:])
print(A[:][0])
print(A[0,:])
print(A[:,0])