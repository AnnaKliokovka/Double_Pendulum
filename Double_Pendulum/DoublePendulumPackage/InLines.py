import DoublePendulumPackage.FindThetasOmega as FTO
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from sympy import diff, symbols, cos, sin
import math
# Глобальные переменные - угол вибрации, частота вибрации
# Global values - 
alpha=0 
# Глобальные переменные 
X = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
Y = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
def f(theta):
       DPE_theta1 =6*omega**2*sin(theta[0]) + (-648*sin(alpha - theta[0])*cos(alpha - theta[0]) + 324*sin(alpha - theta[0])*cos(alpha - theta[1])*cos(theta[0] - theta[1]) - \
            324*sin(theta[0] - theta[1])*cos(alpha - theta[0])*cos(alpha - theta[1]))/(18*cos(theta[0] - theta[1])**2 - 32) + 36*(-324*cos(alpha - theta[0])**2 + \
            324*cos(alpha - theta[0])*cos(alpha - theta[1])*cos(theta[0] - theta[1]) - 144*cos(alpha - theta[1])**2)*sin(theta[0] - theta[1])*cos(theta[0] - theta[1])/(18*cos(theta[0] - theta[1])**2 - 32)**2
       DPE_theta2 = 2*omega**2*sin(theta[1]) + (324*sin(alpha - theta[1])*cos(alpha - theta[0])*cos(theta[0] - theta[1]) - 288*sin(alpha - theta[1])*cos(alpha - theta[1]) + \
            324*sin(theta[0] - theta[1])*cos(alpha - theta[0])*cos(alpha - theta[1]))/(18*cos(theta[0] - theta[1])**2 - 32) - 36*(-324*cos(alpha - theta[0])**2 + \
            324*cos(alpha - theta[0])*cos(alpha - theta[1])*cos(theta[0] - theta[1]) - 144*cos(alpha - theta[1])**2)*sin(theta[0] - theta[1])*cos(theta[0] - theta[1])/(18*cos(theta[0] - theta[1])**2 - 32)**2
       return np.array([DPE_theta1, DPE_theta2],dtype='float64') 
 
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


"""
usage:
    f: system of non-linear equations
    df: derivative of f
    x0: starting guess
    e: desired tolerance of error
    realSolution: the exact solution of the equation
"""
def newtons_method(f, df, x0, e):
    k=0
    delta = np.linalg.norm(np.dot(np.linalg.inv(df(x0)),f(x0)))
    while np.linalg.norm(delta) > e and k<10:
        deltax = np.dot(np.linalg.inv(df(x0)),f(x0))
        x0 = x0 - np.dot(np.linalg.inv(df(x0)),f(x0))
        delta = np.linalg.norm(np.dot(np.linalg.inv(df(x0)),f(x0)))
        k = k+1
    if (k<10):
        return x0
    else:
        return None
 
def calculate(alpha_in):
    global X, Y,alpha, omega
    X = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    Y = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    alpha = alpha_in
    deltaOmega = 0.05
    i_now =0
    for i in [-1,0,1,2]:
        for j in [-1,0,1,2]:
            omega = 0.001
            x0 = np.array([alpha + i*0.5*math.pi,alpha + j*0.5*math.pi])
            while omega <= 10:
                  root = opt.fsolve(f, x0, full_output = 1) 
                  x0_now =np.array([root[0][0],root[0][1]])
                  x_check = np.array([[x0_now[0] - x0[0]],[x0_now[1] - x0[1]]])
                  if max(x_check) < 0.4 and root[2] == 1 :
                         X[i_now] = np.append(X[i_now], x0_now[0])
                         Y[i_now] = np.append(Y[i_now], x0_now[1])
                         x0 = x0_now
                         omega = omega + deltaOmega
                  else:
                      x_star = FTO.calculate(alpha, X[i_now][len(X[i_now])-1], Y[i_now][len(Y[i_now])-1],omega)
                      if (x_star[0] == X[i_now][len(X[i_now])-1] and x_star[1] == Y[i_now][len(Y[i_now])-1] and max(x_check) <= 0.55):
                          X[i_now] = np.append(X[i_now], x0_now[0])
                          Y[i_now] = np.append(Y[i_now], x0_now[1])
                      else:
                          X[i_now] = np.append(X[i_now], x_star[0])
                          Y[i_now] = np.append(Y[i_now], x_star[1])
                          break        
            i_now = i_now+1
    plot_result(i_now)

def plot_result(CountLines):
    plt.figure(figsize=(8,6))
    plt.axis([-math.pi/2 - 0.05 ,3*math.pi/2 + 0.2,-math.pi/2 - 0.2 ,3*math.pi/2 + 0.2])
    for i in range(CountLines):
        j = i+1
        plt.text(X[i][2] - 0.2,Y[i][2]+0.1,'%i' % j)
        plt.plot(X[i],Y[i],'black',linewidth='0.8')        
        plt.plot(X[i][0],Y[i][0],'ko')
    plt.grid(True,linestyle = '--')
    plt.xlabel('theta 1')
    plt.ylabel('theta 2')
    plt.xticks([alpha - 0.5*math.pi,0,alpha, alpha+ 0.5*math.pi, math.pi, alpha+math.pi]) 
    plt.yticks([alpha - 0.5*math.pi,0,alpha, alpha+ 0.5*math.pi, math.pi, alpha+math.pi])
    a = alpha*180/math.pi
    plt.title('Положения равновесия при альфа = %f ' % a)    
    plt.savefig('%i.png' % a)



