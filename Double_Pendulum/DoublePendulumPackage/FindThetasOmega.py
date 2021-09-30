import numpy as np
import matplotlib.pyplot as plt
from sympy import diff, symbols, cos, sin
import math
# Global values - 
alpha=0 
# Глобальные переменные 
Theta1 = np.array([1.7, 2, 1.9, 2.1, 1.1, 4.2])
Theta2 = np.array([-1, -0.8, 1.6, 3.3, 2.2, 1.6])
Omega = np.array([1.4, 1.5, 1 , 0.3, 0.2])
def DPE_theta1_theta1(theta1, theta2, omega):
    return 6*omega**2*cos(theta1) + (-648*sin(alpha - theta1)**2 - 648*sin(alpha - theta1)*sin(theta1 - theta2)*cos(alpha - theta2) + 648*cos(alpha - theta1)**2 - \
        648*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2))/(18*cos(theta1 - theta2)**2 - 32) + (-23328*sin(alpha - theta1)*cos(alpha - theta1) + \
        11664*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 + \
        36*(-648*sin(alpha - theta1)*cos(alpha - theta1) + 324*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - \
        324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 -\
       (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 +\
      (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 + \
      72*(-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3
#Вторая частная производная потенциальной энергии
def DPE_theta1_theta2(theta1, theta2, omega):
    return (324*sin(alpha - theta1)*sin(alpha - theta2)*cos(theta1 - theta2) + 324*sin(alpha - theta1)*sin(theta1 - theta2)*cos(alpha - theta2) - 324*sin(alpha - theta2)*sin(theta1 - theta2)*cos(alpha - theta1) + \
        324*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2))/(18*cos(theta1 - theta2)**2 - 32) - 36*(-648*sin(alpha - theta1)*cos(alpha - theta1) + 324*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) -\
       324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 + (11664*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) -\
      10368*sin(alpha - theta2)*cos(alpha - theta2) + 11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 +\
     (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 -\
    (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 - \
    72*(-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3
#Вторая частная производная потенциальной энергии
def DPE_theta2_theta1(theta1, theta2, omega):
    return (324*sin(alpha - theta1)*sin(alpha - theta2)*cos(theta1 - theta2) + 324*sin(alpha - theta1)*sin(theta1 - theta2)*cos(alpha - theta2) - 324*sin(alpha - theta2)*sin(theta1 - theta2)*cos(alpha - theta1) +\
       324*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2))/(18*cos(theta1 - theta2)**2 - 32) - (-23328*sin(alpha - theta1)*cos(alpha - theta1) + 11664*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) -\
      11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 + 36*(324*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) -\
     288*sin(alpha - theta2)*cos(alpha - theta2) + 324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 +\
    (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 -\
   (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 - \
   72*(-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3
#Вторая частная производная потенциальной энергии
def DPE_theta2_theta2(theta1, theta2, omega):
    return 2*omega**2*cos(theta2) + (-288*sin(alpha - theta2)**2 + 648*sin(alpha - theta2)*sin(theta1 - theta2)*cos(alpha - theta1) - 648*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) + \
        288*cos(alpha - theta2)**2)/(18*cos(theta1 - theta2)**2 - 32) - 36*(324*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) - 288*sin(alpha - theta2)*cos(alpha - theta2) + \
        324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 - (11664*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) -\
       10368*sin(alpha - theta2)*cos(alpha - theta2) + 11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 -\
      (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 +\
     (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 +\
    72*(-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3

def f(x):
        minor2 = DPE_theta1_theta1(x[0,0], x[1,0],x[2,0])*DPE_theta2_theta2(x[0,0], x[1,0],x[2,0]) - DPE_theta1_theta2(x[0,0], x[1,0],x[2,0])*DPE_theta2_theta1(x[0,0], x[1,0],x[2,0])
        DPE_theta1 =6*x[2,0]**2*sin(x[0,0]) + (-648*sin(alpha - x[0,0])*cos(alpha - x[0,0]) + 324*sin(alpha - x[0,0])*cos(alpha - x[1,0])*cos(x[0,0] - x[1,0]) - \
            324*sin(x[0,0] - x[1,0])*cos(alpha - x[0,0])*cos(alpha - x[1,0]))/(18*cos(x[0,0] - x[1,0])**2 - 32) + 36*(-324*cos(alpha - x[0,0])**2 + \
            324*cos(alpha - x[0,0])*cos(alpha - x[1,0])*cos(x[0,0] - x[1,0]) - 144*cos(alpha - x[1,0])**2)*sin(x[0,0] - x[1,0])*cos(x[0,0] - x[1,0])/(18*cos(x[0,0] - x[1,0])**2 - 32)**2
        DPE_theta2 = 2*x[2,0]**2*sin(x[1,0]) + (324*sin(alpha - x[1,0])*cos(alpha - x[0,0])*cos(x[0,0] - x[1,0]) - 288*sin(alpha - x[1,0])*cos(alpha - x[1,0]) + \
            324*sin(x[0,0] - x[1,0])*cos(alpha - x[0,0])*cos(alpha - x[1,0]))/(18*cos(x[0,0] - x[1,0])**2 - 32) - 36*(-324*cos(alpha - x[0,0])**2 + \
            324*cos(alpha - x[0,0])*cos(alpha - x[1,0])*cos(x[0,0] - x[1,0]) - 144*cos(alpha - x[1,0])**2)*sin(x[0,0] - x[1,0])*cos(x[0,0] - x[1,0])/(18*cos(x[0,0] - x[1,0])**2 - 32)**2
        return np.array([[DPE_theta1], [DPE_theta2], [minor2]],dtype='float64')


def minor2():
    theta1, theta2, omega = symbols('theta1 theta2 omega')
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

def DPE_theta1_omega(x):
    return 12*x[2,0]*sin(x[0,0])

def DPE_theta2_omega(x):
    return 4*x[2,0]*sin(x[1,0])

def minor2_theta1(x):
    theta1, theta2, omega = symbols('theta1, theta2, omega')
    DM_theta1 = diff(minor2(),theta1)
    return  DM_theta1.subs([(theta1, x[0,0]), (theta2, x[1,0]), (omega, x[2,0])])

def minor2_theta2(x):
     theta1, theta2, omega = symbols('theta1, theta2, omega')
     DM_theta2 = diff(minor2(),theta2)
     return  DM_theta2.subs([(theta1, x[0,0]), (theta2, x[1,0]), (omega, x[2,0])])

def minor2_omega(x):
     theta1, theta2, omega = symbols('theta1, theta2, omega')
     DM_theta_omega = diff(minor2(),omega)
     return  DM_theta_omega.subs([(theta1, x[0,0]), (theta2, x[1,0]), (omega, x[2,0])])

def df(x):
    return np.array([[DPE_theta1_theta1(x[0,0],x[1,0], x[2,0]), DPE_theta1_theta2(x[0,0],x[1,0], x[2,0]), DPE_theta1_omega(x)],
                     [DPE_theta2_theta1(x[0,0],x[1,0], x[2,0]), DPE_theta2_theta2(x[0,0],x[1,0], x[2,0]), DPE_theta2_omega(x)],
                     [minor2_theta1(x), minor2_theta2(x), minor2_omega(x)]],dtype='float64')



"""
usage:
    f: system of non-linear equations
    df: derivative of f
    x0: starting guess
    e: desired tolerance of error
    realSolution: the exact solution of the equation
"""
def newtons_method(f, df, x0, e):
    print(df(x0))
    print(f(x0))
    delta = np.linalg.norm(np.dot(np.linalg.inv(df(x0)),f(x0)))
    while np.linalg.norm(delta) > e:
        x0 = x0 - np.dot(np.linalg.inv(df(x0)),f(x0))
        delta = np.linalg.norm(np.dot(np.linalg.inv(df(x0)),f(x0)))
    return x0

def calculate(alpha_in):
    global  alpha, Theta1, Theta2, Omega
    alpha = alpha_in
    x0 = np.array([[Theta1[0]],[Theta2[0]],[Omega[0]]])
    x0_now = newtons_method(f, df, x0 ,10**-5)
    print(x0_now)

