import numpy as np
import matplotlib.pyplot as plt

L = 4
d1 = np.arange(0, L**3, 1)
d3 = np.reshape(d1, (L,L,L))

#print(d1)
#print(d3)

d1 = np.reshape(d3, L**3)
#print(d1)

# nearest neighbours:
# x: [z][y][x-1], [z][y][x+1]
# y: [z][y-1][x], [z][y+1][x]
# z: [z-1][y][x], [z+1][y][x]

i = (0,0,0)

def idx(i):
    global L 
    return i % L 

def getNN(a, i):
    x,y,z = i
    return np.array([
        a[z][y][x-1], a[z][y][idx(x+1)],
        a[z][y-1][x], a[z][idx(y+1)][x],
        a[z-1][y][x], a[idx(z+1)][y][x]
    ])
allNN = []
for x in range(L):
    for y in range(L):
        for z in range(L):
            allNN.append(getNN(d3,(x,y,z)))
allNN = np.array(allNN)
M = np.zeros((L**3,L**3))
i = 0
for nnArray in allNN:
    for index in nnArray:
        M[i][index] = 1
    i+=1
#print(M)
i = (0,0,0)
#plt.matshow(M)
#plt.show()
print(d3)
d3 = 1/(d3+1)
print(d3)
print(sorted(np.reshape(d3, L**3)))
print(d3[i[0]][i[1]][i[2]])