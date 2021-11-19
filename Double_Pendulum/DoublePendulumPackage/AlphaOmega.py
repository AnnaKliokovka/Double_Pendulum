import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from sympy import diff, symbols, cos, sin
import math

class OmegaAlpha:
    alpha=0
    omegas = []
    counts=[]
    def __init__(self, alpha_in, omega_in):
       self.alpha = alpha_in
       for i in range(len(omega_in)):
           for j in range(len(omega_in[i])):
                if omega_in[i][j] in self.omegas:
                    index = list(self.omegas).index(omega_in[i][j])
                    self.counts[index] = self.counts[index]+1
                else:
                    self.omegas = np.append(self.omegas,omega_in[i][j])
                    self.counts = np.append(self.counts,1)

    def Omega():
        return omegas
    def Counts():
        return counts

