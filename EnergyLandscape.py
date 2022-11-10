import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

class EnergyLandscape:

    energies = None
    L = None
    chemPot = None

    def __init__(self, L, pattern="checker") -> None:
        self.L = L
        if pattern == "checker":
            self.createCheckerPattern()

    def createCheckerPattern(self, scale=1):
        if self.L%2 == 0:
            raise ValueError("Checker pattern needs odd length")
        
        self.energies = np.array([i%2 for i in range(self.L**3)])
        self.energies = np.reshape(self.energies, (self.L,self.L,self.L))

    def getJumpRates(self, model, attemptFreq, beta, minJumpBarrier=0):
        if model == "CSP":
            saddlePointEnergy = max(np.reshape(self.energies, self.L**3))+ minJumpBarrier
            return attemptFreq * np.exp(-beta*saddlePointEnergy) * np.exp(beta*self.energies)
        else:
            return 0

    def getEqOccupationNumbers(self, c, beta, nBins=1000):
        """Returns mean occupation numbers of sites in equilibrium"""
        return 1/(np.exp(beta*(self.energies - self.getChemPot(c, beta, nBins)))+1)

    def getChemPot(self, c, beta, nBins=1000):
        """Determines the chemical potential of the system,
        c.f. eq.11a in "Hopping conductivity in ideal Fermi lattice gases
        with site energy disorder" """

        if self.chemPot:
            return self.chemPot

        # calaculate site energy distribution g(E)
        sortedEnergies = np.reshape(self.energies, self.L**3)
        energyDist, bins = np.histogram(sortedEnergies, bins=nBins, density=True) # energyDist ^= g(E), bins ^= E

        # equation determining mu
        def approx(mu):
            dEBin = (bins[1]-bins[0])/2 # offset bin edge value - site energy value
            sum = 0
            for i in range(len(bins)-1):
                sum += energyDist[i] / (np.exp(beta*(bins[i]+dEBin - mu)) + 1)
            
            return abs(sum - c)

        # optimize result of approximation of c (|sum(mu) - c| = 0) varying mu
        solution = optimize.minimize_scalar(approx, method='brent')
        print(f'Chemical potential: {solution.x}')
        print(f'Root finding converged?: {solution.success}')
        print(f'Cause of termination of root finding: {solution.message}')

        self.chemPot = solution.x

        return solution.x

if __name__ == "__main__":
    el = EnergyLandscape(3, "checker")
    el.getChemPot(1, 1)
    #print(el.getJumpRates("CSP", 1, 1, minJumpBarrier=1))