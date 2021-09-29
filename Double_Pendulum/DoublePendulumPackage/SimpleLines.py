import numpy as np
import matplotlib.pyplot as plt
from sympy import diff, symbols, cos, sin
import math
# Глобальные переменные - угол вибрации, частота вибрации
# Global values - 
alpha=0 
omega = 10
# Глобальные переменные 
X = np.array([0])
Y = np.array([0])

X1 = np.array([0])
Y1 = np.array([math.pi])

X2 = np.array([math.pi])
Y2 = np.array([math.pi])

X3 = np.array([math.pi])
Y3 = np.array([0])
def f(theta):
        DPE_theta1 =6*omega**2*sin(theta[0,0]) + (-648*sin(alpha - theta[0,0])*cos(alpha - theta[0,0]) + 324*sin(alpha - theta[0,0])*cos(alpha - theta[1,0])*cos(theta[0,0] - theta[1,0]) - \
            324*sin(theta[0,0] - theta[1,0])*cos(alpha - theta[0,0])*cos(alpha - theta[1,0]))/(18*cos(theta[0,0] - theta[1,0])**2 - 32) + 36*(-324*cos(alpha - theta[0,0])**2 + \
            324*cos(alpha - theta[0,0])*cos(alpha - theta[1,0])*cos(theta[0,0] - theta[1,0]) - 144*cos(alpha - theta[1,0])**2)*sin(theta[0,0] - theta[1,0])*cos(theta[0,0] - theta[1,0])/(18*cos(theta[0,0] - theta[1,0])**2 - 32)**2
        DPE_theta2 = 2*omega**2*sin(theta[1,0]) + (324*sin(alpha - theta[1,0])*cos(alpha - theta[0,0])*cos(theta[0,0] - theta[1,0]) - 288*sin(alpha - theta[1,0])*cos(alpha - theta[1,0]) + \
            324*sin(theta[0,0] - theta[1,0])*cos(alpha - theta[0,0])*cos(alpha - theta[1,0]))/(18*cos(theta[0,0] - theta[1,0])**2 - 32) - 36*(-324*cos(alpha - theta[0,0])**2 + \
            324*cos(alpha - theta[0,0])*cos(alpha - theta[1,0])*cos(theta[0,0] - theta[1,0]) - 144*cos(alpha - theta[1,0])**2)*sin(theta[0,0] - theta[1,0])*cos(theta[0,0] - theta[1,0])/(18*cos(theta[0,0] - theta[1,0])**2 - 32)**2
        return np.array([[DPE_theta1], [DPE_theta2]],dtype='float64') 
 
#Вторая частная производная потенциальной энергии
def DPE_theta1_theta1(theta1, theta2):
    return 6*omega**2*cos(theta1) + (-648*sin(alpha - theta1)**2 - 648*sin(alpha - theta1)*sin(theta1 - theta2)*cos(alpha - theta2) + 648*cos(alpha - theta1)**2 - \
        648*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2))/(18*cos(theta1 - theta2)**2 - 32) + (-23328*sin(alpha - theta1)*cos(alpha - theta1) + \
        11664*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 + \
        36*(-648*sin(alpha - theta1)*cos(alpha - theta1) + 324*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - \
        324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 -\
       (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 +\
      (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 + \
      72*(-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3
#Вторая частная производная потенциальной энергии
def DPE_theta1_theta2(theta1, theta2):
    return (324*sin(alpha - theta1)*sin(alpha - theta2)*cos(theta1 - theta2) + 324*sin(alpha - theta1)*sin(theta1 - theta2)*cos(alpha - theta2) - 324*sin(alpha - theta2)*sin(theta1 - theta2)*cos(alpha - theta1) + \
        324*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2))/(18*cos(theta1 - theta2)**2 - 32) - 36*(-648*sin(alpha - theta1)*cos(alpha - theta1) + 324*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) -\
       324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 + (11664*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) -\
      10368*sin(alpha - theta2)*cos(alpha - theta2) + 11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 +\
     (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 -\
    (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 - \
    72*(-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3
#Вторая частная производная потенциальной энергии
def DPE_theta2_theta1(theta1, theta2):
    return (324*sin(alpha - theta1)*sin(alpha - theta2)*cos(theta1 - theta2) + 324*sin(alpha - theta1)*sin(theta1 - theta2)*cos(alpha - theta2) - 324*sin(alpha - theta2)*sin(theta1 - theta2)*cos(alpha - theta1) +\
       324*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2))/(18*cos(theta1 - theta2)**2 - 32) - (-23328*sin(alpha - theta1)*cos(alpha - theta1) + 11664*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) -\
      11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 + 36*(324*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) -\
     288*sin(alpha - theta2)*cos(alpha - theta2) + 324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 +\
    (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 -\
   (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 - \
   72*(-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3
#Вторая частная производная потенциальной энергии
def DPE_theta2_theta2(theta1, theta2):
    return 2*omega**2*cos(theta2) + (-288*sin(alpha - theta2)**2 + 648*sin(alpha - theta2)*sin(theta1 - theta2)*cos(alpha - theta1) - 648*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) + \
        288*cos(alpha - theta2)**2)/(18*cos(theta1 - theta2)**2 - 32) - 36*(324*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) - 288*sin(alpha - theta2)*cos(alpha - theta2) + \
        324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 - (11664*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) -\
       10368*sin(alpha - theta2)*cos(alpha - theta2) + 11664*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2 -\
      (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 +\
     (-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**2 +\
    72*(-11664*cos(alpha - theta1)**2 + 11664*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 5184*cos(alpha - theta2)**2)*sin(theta1 - theta2)**2*cos(theta1 - theta2)**2/(18*cos(theta1 - theta2)**2 - 32)**3
def df(theta):
    return np.array([[DPE_theta1_theta1(theta[0,0],theta[1,0]), DPE_theta1_theta2(theta[0,0],theta[1,0])],
     [DPE_theta2_theta1(theta[0,0],theta[1,0]), DPE_theta2_theta2(theta[0,0],theta[1,0])]],dtype='float64')

def g(x):
    return np.array([[2*x[0,0] + 3*x[1,0] + 4*x[1,0]**2 - 5*x[1,0]**3 + 16],[-x[0,0] - 2*x[1,0] + 3*x[1,0]**2 - 7*x[1,0]**3 + 49]])

def dg(x):
    return np.array([[2, 3 - 8*x[1,0] - 15*x[1,0]**2],[-1, -2 + 6*x[1,0] - 21*x[1,0]**2]])


"""
usage:
    f: system of non-linear equations
    df: derivative of f
    x0: starting guess
    e: desired tolerance of error
    realSolution: the exact solution of the equation
"""
def newtons_method(f, df, x0, e):
    delta = np.linalg.norm(np.dot(np.linalg.inv(df(x0)),f(x0)))
    while np.linalg.norm(delta) > e:
        deltax = np.dot(np.linalg.inv(df(x0)),f(x0))
        x0 = x0 - np.dot(np.linalg.inv(df(x0)),f(x0))
        delta = np.linalg.norm(np.dot(np.linalg.inv(df(x0)),f(x0)))
    return x0
 
def calculate(alpha_in):
    global X, Y, X1, Y1, X2, Y2, X3, Y3, alpha, omega
    alpha = alpha_in
    omega = 10;
    x0 = np.array([[0],[0]])
#Первая линия от (0,0)
    while omega >= 0:
          x0_now = newtons_method(f, df, x0, 10**-10)  
          X = np.append(X, x0_now[0])
          Y = np.append(Y, x0_now[1])
          x0 = x0_now
          omega = omega - 0.1
#Вторая линия из (0, Pi)
    x0 = np.array([[0],[math.pi]])
    omega = 10
    while omega >= 0:
          x0_now = newtons_method(f, df, x0, 10**-10)  
          X1 =np.append(X1, x0_now[0])
          Y1 = np.append(Y1, x0_now[1])
          x0 = x0_now
          omega = omega - 0.1
#Третья линия из (Pi, Pi)
    x0 = np.array([[math.pi],[math.pi]])
    omega = 10
    while omega >= 0:
          x0_now = newtons_method(f, df, x0, 10**-10)  
          X2 =np.append(X2, x0_now[0])
          Y2 = np.append(Y2, x0_now[1])
          x0 = x0_now
          omega = omega - 0.1
#Четвертая линия из (Pi, 0)
    x0 = np.array([[math.pi],[0]])
    omega = 10
    while omega >= 0:
          x0_now = newtons_method(f, df, x0, 10**-10)  
          X3 =np.append(X3, x0_now[0])
          Y3 = np.append(Y3, x0_now[1])
          x0 = x0_now
          omega = omega - 0.1
def plot_result():
    plt.axis([-1,4,-1,4.2])
    plt.plot(X,Y,"salmon")
    plt.plot(X1,Y1,"lime")
    plt.plot(X2,Y2,"pink")
    plt.plot(X3,Y3,"blue")
    plt.grid(True)
    plt.xlabel("Theta 1"); plt.ylabel("Theta 2")
    plt.legend(['First series','Second series','3 series','4 series'], loc=2)
    plt.show()


