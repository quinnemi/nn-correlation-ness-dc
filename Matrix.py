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
from numba import njit
import configparser
import timeit
from pathlib import Path
#
config = configparser.ConfigParser()
config.read('config.ini')
L = int(config['SYSTEM PARAMETERS']['L'])
C = float(config['SYSTEM PARAMETERS']['c'])
BETA = float(config['SYSTEM PARAMETERS']['beta'])
MINJUMPBARRIER = float(config['SYSTEM PARAMETERS']['minJumpBarrier'])
ATTEMPTFREQ = float(config['SYSTEM PARAMETERS']['attemptFreq'])

#energyLandscape = EnergyLandscape(L, pattern="checker") # 3D matrix
#jumpRates = energyLandscape.getJumpRates("CSP", ATTEMPTFREQ, BETA, minJumpBarrier=MINJUMPBARRIER) # 3D matrix
#eqOccNum = energyLandscape.getEqOccupationNumbers(C, BETA) # 3D matrix

@njit(cache=True)
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

        #print(logK1)
    
    return K

@njit(cache=True)
def acc(array: np.array, tuple: tuple) -> float:
    """Access 3D array via 3-tuple"""
    i = tuple[0]
    j = tuple[1]
    k = tuple[2]
    a = array[int(i)][int(j)][int(k)]
    
    return a

@njit(cache=True)
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

@njit(cache=True)
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

@njit(cache=True)
def K3(i: tuple, j: tuple) -> float:
    """Calculation of a single entry of K2 (i,j are next nearest neighbors)
        CSP model for jump rates is implied"""
    
    term1_1 = sum([(acc(jumpRates, i) - acc(jumpRates, k)) / D(i,j) * acc(jumpRates, j) * (1-acc(eqOccNum, k))
                    for k in getNNidx(i)])
    term1 = acc(eqOccNum, i) / (1-acc(eqOccNum, j)) * term1_1

    return term1

@njit(cache=True)
def D(pos1: tuple, pos2: tuple) -> float:
    """Denominator term in eq. 17s
        pos1: i, pos2: j"""
    
    return 1/acc(eqOccNum, pos1) * (
        sum([acc(jumpRates, p) for p in getNNidx(pos1) if not eqIdx(p, pos2)])
    )   +  1/acc(eqOccNum, pos2) * (
        sum([acc(jumpRates, p) for p in getNNidx(pos2) if eqIdx(p, pos1)])
    )


############################## Supporting functions ##############################

@njit(cache=True)
def eqIdx(i: tuple, j: tuple) -> bool:
    """Check indexes for equality"""
    res = [0, 0, 0]
    for n in range(3):
        res[n] = (i[n] == j[n])
    
    return res[0] and res[1] and res[2]

@njit(cache=True)
def pIdx(j: int) -> int:
        """Returns periodic index"""
        global L 
        return int(j % L)

@njit(cache=True)
def getNN(a: np.array, i: tuple) -> np.array:
    """Returns entries of nearest neighbors in a periodic 3D array 'a' at position i"""
    
    x,y,z = i
    return np.array([
        a[z][y][pIdx(x-1)], a[z][y][pIdx(x+1)],
        a[z][pIdx(y-1)][x], a[z][pIdx(y+1)][x],
        a[pIdx(z-1)][y][x], a[pIdx(z+1)][y][x]
    ])

@njit(cache=True)
def getNNidx(i: tuple) -> np.array:
    """Returns indexes of nearest neighbors in a periodic 3D array of size L at position i"""
    
    x,y,z = i
    x = int(x)
    y = int(y)
    z = int(z)
    return np.array([
        (pIdx(x-1),y,z), (pIdx(x+1),y,z),
        (x,pIdx(y-1),z), (x,pIdx(y+1),z),
        (x,y,pIdx(z-1)), (x,y,pIdx(z+1))
    ], dtype=np.int32)

@njit(cache=True)
def getNNNidx(i: tuple) -> np.array:
    """Returns indexes of next nearest neighbors in a periodic 3D array of size L at position i"""

    # gather all candidates for nnn by collecting nn of nn of i
    allNNN = np.zeros((36, 3), dtype=np.int32)
    idx = 0
    for nnOfI in getNNidx(i):
        for nnOfnn in getNNidx(nnOfI):
            allNNN[idx] = nnOfnn
            idx += 1

    #print(f'allNNN: {allNNN}')

    #allNNN = np.array([getNNidx(j) for j in getNNidx(i)])
    #allNNN = np.reshape(allNNN, (36, 3))
    filteredNNN = np.ones((18, 3), dtype=np.int32)*9999
    idxJ = 0
    for j in allNNN:
        if (not eqIdx(i, j)): # exclude central site

            # only allow unique indexes
            jInFilteredNNN = False
            for k in filteredNNN:
                if eqIdx(j, k):
                    jInFilteredNNN = True
                    break

            if not jInFilteredNNN:
                filteredNNN[idxJ] = j
                idxJ += 1 
    return filteredNNN

@njit(cache=True)
def tupleToIndex(t: tuple) -> int:
    return t[0] + t[1]*L + t[2]*L**2

#print(timeit.timeit(K, number=100))
#mat = K()
#print(np.linalg.det(mat))
#plt.matshow(mat)
#plt.show()

#with open(Path(__file__).parent.parent.parent / "data" / "K_numba", 'w+') as file:
#    file.write(str(mat))

a = np.array([[[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]],
              [[10.,11.,12.],[13.,14.,15.],[16.,17.,18.]],
              [[19.,20.,21.],[22.,23.,24.],[25.,26.,27.]]])
i = (1,0,2)
j = (1,0,0)
k = (-3,4,5)
print(f'acc: {a, i}')
print(f'eqIdx: {eqIdx(i,i)}, {eqIdx(i,j)}')
print(f'pIdx: {pIdx(k[0])}, {pIdx(k[1])}, {pIdx(k[2])}')
print(f'getNNidx: {getNNidx(i)}')
print(f'getNNNidx: {getNNNidx(i)}')
print(f'tupleToIndex: {tupleToIndex(i)}')