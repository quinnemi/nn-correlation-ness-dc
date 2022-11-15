import numpy as np

M = np.array([
    [1,2,3],
    [4,5,6],
    [7,8,9]
])

b = np.array([10,20,30])

#x,y,z = np.linalg.solve(M,b)
#print((x,y,z))
#print(1*x + 2*y + 3*z)
#print(4*x + 5*y + 6*z)
#print(7*x + 8*y + 9*z)

a = np.longdouble(1e-309)
b = 1e-309
print(str(a))
print(b)
print(np.finfo(np.longdouble))