import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from sympy import diff, symbols, cos, sin
import math
# Глобальные переменные - угол вибрации, частота вибрации
# Global values - 
alpha=0 

def minor1(theta1_in, theta2_in, omega_in):
    theta1, theta2, omega = symbols('theta1 theta2 omega')
    DPE_theta1 =6*omega**2*sin(theta1) + (-648*sin(alpha - theta1)*cos(alpha - theta1) + 324*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - \
            324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))/(18*cos(theta1 - theta2)**2 - 32) + 36*(-324*cos(alpha - theta1)**2 + \
            324*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 144*cos(alpha - theta2)**2)*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2
        
    minor1= diff(DPE_theta1,theta1)
    return  minor1.subs([(theta1, theta1_in), (theta2, theta2_in), (omega, omega_in)])


def minor2(theta1, theta2, omega):
    return (6*omega**2*cos(theta1) + (-648*sin(alpha - theta1)**2 - 648*sin(alpha - theta1)*sin(theta1 - theta2)*cos(alpha - theta2) + 648*cos(alpha - theta1)**2 - \
        648*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2))/(18*cos(theta1 - theta2)**2 - 32) + 2*(-23328*sin(alpha - theta1)*cos(alpha - theta1) + \
        11664*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 -\
       (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 +\
      (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 +\
     (-839808*cos(alpha - theta1)**2 + 839808*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 373248*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3)*(2*omega**2*cos(theta2) +\
    (-288*sin(alpha - theta2)**2 + 648*sin(alpha - theta2)*sin(theta1 - theta2)*cos(alpha - theta1) - 648*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) + 288*cos(alpha - theta2)**2)/(18*cos(theta1 - theta2)**2 - 32) -\
   2*(11664*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) - 10368*sin(alpha - theta2)*cos(alpha - theta2) + 11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 - \
   (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 + (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) -\
  5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 + (-839808*cos(alpha - theta1)**2 + 839808*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 373248*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3) - \
  ((324*sin(alpha - theta1)*sin(alpha - theta2)*cos(theta1 - theta2) + 324*sin(alpha - theta1)*sin(theta1 - theta2)*cos(alpha - theta2) - 324*sin(alpha - theta2)*sin(theta1 - theta2)*cos(alpha - theta1) + 324*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2))/(18*cos(theta1 - theta2)**2 - 32) - (-23328*sin(alpha - theta1)*cos(alpha - theta1) + \
  11664*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 + (11664*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) - 10368*sin(alpha - theta2)*cos(alpha - theta2) + \
  11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 + (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 -\
 (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 - (-839808*cos(alpha - theta1)**2 + 839808*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 373248*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3)**2

def stability(alpha_in, theta1,theta2,omega):
    global alpha
    alpha = alpha_in
    minor11=(minor1(theta1,theta2,omega))
    minor22=round(minor2(theta1,theta2,omega))
    if minor11>0 and minor22>0 :
        return 1
    else:
        return 0