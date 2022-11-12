"""
Calculates the matrix K (K * kappa = h), c.f. notes on correlation function, eq. 25a)
Used to calculate kappa_i which determine the conductivity

Elements of K can be grouped into three cases:
- K1: diagonal matrix elements (indices i,j refer to equal sites)
- K2: matrix element denoting nearest neighbor j of site i
- K3: matrix element denoting next nearest neighbor j of site i
"""

import numpy as np
from EnergyLandscape import EnergyLandscape
import matplotlib.pyplot as plt
import timeit

L = 3
c = 0.1
beta = 1
minJumpBarrier = 1
attemptFreq = 1

energyLandscape = EnergyLandscape(L, pattern="checker") # 3D matrix
jumpRates = energyLandscape.getJumpRates("CSP", attemptFreq, beta, minJumpBarrier=minJumpBarrier) # 3D matrix
chemPot = energyLandscape.getChemPot(c, beta) # scalar value
eqOccNum = energyLandscape.getEqOccupationNumbers(c, beta) # 3D matrix


def K() -> np.array:
    idxList = []
    for z in range(L):
        for y in range(L):
            for x in range(L):
                idxList.append((x,y,z))
    idxList = np.array(idxList)

    K = np.zeros((L**3, L**3))
    for i in idxList:
        nn = getNNidx(i)
        nnn = getNNNidx(i)

        K[tupleToIndex(i)][tupleToIndex(i)] = logK1 = K1(i)
        for j in nn:
            K[tupleToIndex(i)][tupleToIndex(j)] = logK2 = K2(i, j)
            # print(logK2)
        for j in nnn:
            K[tupleToIndex(i)][tupleToIndex(j)] = logK3 = K3(i, j)
            # print(logK3)

        # print(logK1)
    
    return K
    
def K1(i: tuple) -> float:
    """Calculation of a single entry of K1 (i,j are equal sites)
        CSP model for jump rates is implied"""
    global jumpRates

    term1 = -6*acc(jumpRates, i)

    term2 = 0
    for k in getNNidx(i):
        term2_1 = 0
        for l in getNNidx(k):
            if not eqIdx(l, i):
                term2_1 += acc(jumpRates, l) * acc(eqOccNum, l)

        term2 += ((acc(jumpRates, i) - acc(jumpRates, k)) / D(i,k)) * term2_1

    return term1 + term2

def K2(i: tuple, j: tuple) -> float:
    """Calculation of a single entry of K2 (i,j are nearest neighbors)
        CSP model for jump rates is implied"""

    term1 = acc(jumpRates, j)

    term2_1 = sum([acc(jumpRates, k) * acc(eqOccNum, k)
                    for k in getNNidx(i) if not eqIdx(k, j)])
    term2 = (acc(jumpRates, i) - acc(jumpRates, j)) / D(i,j) * term2_1

    term3_1 = sum([(acc(jumpRates, i) - acc(jumpRates, k)) / D(i,k) * acc(eqOccNum, k)
                    for k in getNNidx(i) if not eqIdx(k, j)])
    term3 = acc(jumpRates, j) * (1-acc(eqOccNum, i)) / (1-acc(eqOccNum, j)) * term3_1

    return term1 + term2 + term3

def K3(i: tuple, j: tuple) -> float:
    """Calculation of a single entry of K2 (i,j are next nearest neighbors)
        CSP model for jump rates is implied"""
    
    term1_1 = sum([(acc(jumpRates, i) - acc(jumpRates, k)) / D(i,j) * acc(jumpRates, j) * (1-acc(eqOccNum, k))
                    for k in getNNidx(i)])
    term1 = acc(eqOccNum, i) / (1-acc(eqOccNum, j)) * term1_1

    return term1

def D(pos1: tuple, pos2: tuple) -> float:
    """Denominator term in eq. 17s
        pos1: i, pos2: j"""
    
    return 1/acc(eqOccNum, pos1) * (
        sum([acc(jumpRates, p) for p in getNNidx(pos1) if not eqIdx(p, pos2)])
    )   +  1/acc(eqOccNum, pos2) * (
        sum([acc(jumpRates, p) for p in getNNidx(pos2) if eqIdx(p, pos1)])
    )


############################## Supporting functions ##############################


def acc(array: np.array, tuple: tuple) -> float:
    """Access 3D array via 3-tuple"""
    return array[tuple[0]][tuple[1]][tuple[2]]

def eqIdx(i: tuple, j: tuple) -> bool:
    """Check indexes for equality"""
    return np.all([a == b for a,b in list(zip(i,j))])

def pIdx(j: int) -> int:
        """Returns periodic index"""
        global L 
        return j % L 

def getNN(a: np.array, i: tuple) -> np.array:
    """Returns entries of nearest neighbors in a periodic 3D array 'a' at position i"""
    
    x,y,z = i
    return np.array([
        a[z][y][x-1], a[z][y][pIdx(x+1)],
        a[z][y-1][x], a[z][pIdx(y+1)][x],
        a[z-1][y][x], a[pIdx(z+1)][y][x]
    ])

def getNNidx(i: tuple) -> np.array:
    """Returns indexes of nearest neighbors in a periodic 3D array of size L at position i"""
    
    x,y,z = i
    return np.array([
        (x-1,y,z), (pIdx(x+1),y,z),
        (x,y-1,z), (x,pIdx(y+1),z),
        (x,y,z-1), (x,y,pIdx(z+1))
    ])

def getNNNidx(i: tuple) -> np.array:
    """Returns indexes of next nearest neighbors in a periodic 3D array of size L at position i"""

    allNNN = np.array([getNNidx(j) for j in getNNidx(i)])
    allNNN = np.reshape(allNNN, (36, 3))
    filteredNNN = []
    for j in allNNN:
        if (not eqIdx(i, j)): # exclude central site
            # only allow unique indexes
            jInFilteredNNN = False
            for k in filteredNNN:
                if eqIdx(j, k):
                    jInFilteredNNN = True
                    break
            if not jInFilteredNNN:
                filteredNNN.append(j)
    
    return np.array(filteredNNN)

def isInList(i: tuple, list: np.array):
    return np.any([eqIdx(i, list[j]) for j in range(len(list))])

def tupleToIndex(t: tuple) -> int:
    return t[0] + t[1]*L + t[2]*L**2

#print(timeit.timeit(K, number=10))
mat = K()
print(np.linalg.det(mat))
plt.matshow(mat)
plt.show()