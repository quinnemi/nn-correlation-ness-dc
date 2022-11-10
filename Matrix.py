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

L = 3
c = 0.1
beta = 1
minJumpBarrier = 1
attemptFreq = 1

energyLandscape = EnergyLandscape(L, pattern="checker") # 3D matrix
jumpRates = energyLandscape.getJumpRates("CSP", attemptFreq, beta, minJumpBarrier=minJumpBarrier) # 3D matrix
chemPot = energyLandscape.getChemPot(c, beta) # scalar value
eqOccNum = energyLandscape.getEqOccupationNumbers(c, beta) # 3D matrix

def K1Single(i: tuple):
    """Calculation of a single entry of K1, CSP model for jump rates is implied"""
    global jumpRates

    term1 = -6*acc(jumpRates, i)
    term2 = 0
    for k in getNN(jumpRates, i, getIdx=True):
        term2_1 = 0
        for l in getNN(jumpRates, k, getIdx=True):
            if not np.all([l, i]):
                term2_1 += acc(jumpRates, l) * acc(eqOccNum, l)

        term2 += ((acc(jumpRates, i) - acc(jumpRates, k)) / D(i,k)) * term2_1

    return term1 + term2

def D(pos1: tuple, pos2: tuple):
    """Denominator term in eq. 17
        pos1: i, pos2: j"""
    
    return 1/acc(eqOccNum, pos1) * (
        sum([acc(jumpRates, p) for p in getNN(jumpRates, pos1, getIdx=True) if not np.all([p, pos2])])
    )   +  1/acc(eqOccNum, pos2) * (
        sum([acc(jumpRates, p) for p in getNN(jumpRates, pos2, getIdx=True) if np.all([p, pos1])])
    )
    
def acc(array, tuple):
    """Access 3D array via 3-tuple"""
    return array[tuple[0]][tuple[1]][tuple[2]]

def getNN(a, i, getIdx=False):
    """Returns entries of nearest neighbors in a periodic 3D array"""

    def idx(j):
        """Returns periodic index"""
        global L 
        return j % L 
    
    x,y,z = i

    if getIdx:
        return np.array([
            (x-1,y,z), (idx(x+1),y,z),
            (x,y-1,z), (x,idx(y+1),z),
            (x,y,z-1), (x,y,idx(z+1))
        ])

    return np.array([
        a[z][y][x-1], a[z][y][idx(x+1)],
        a[z][y-1][x], a[z][idx(y+1)][x],
        a[z-1][y][x], a[idx(z+1)][y][x]
    ])

print(K1Single((1,1,1)))