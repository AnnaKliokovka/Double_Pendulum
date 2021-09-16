import math
from sympy import diff, symbols, cos, sin
import matplotlib.pyplot as plt
theta1, theta2, omega, alpha = symbols('theta1 theta2 omega alpha')
#Получение частной производной потенциальной энергии от угла тета1
def DPE_theta1():
    PotentialEnergy = -2*((3*omega*omega)*cos(theta1)+(omega*omega)*cos(theta2))+\
                      (4*(81*cos(theta1-theta2)*cos(alpha-theta1)*cos(alpha-theta2)-81*cos(alpha-theta1)*cos(alpha-theta1)- \
                      36*cos(alpha-theta2)*cos(alpha-theta2)))/(18*cos(theta1-theta2)*cos(theta1-theta2)-32)
    return diff(PotentialEnergy, theta1)
#Получение частной производной потенциальной энергии от угла тета2
def DPE_theta2():
    PotentialEnergy = -2*((3*omega*omega)*cos(theta1)+(omega*omega)*cos(theta2))+\
                      (4*(81*cos(theta1-theta2)*cos(alpha-theta1)*cos(alpha-theta2)-81*cos(alpha-theta1)*cos(alpha-theta1)- \
                      36*cos(alpha-theta2)*cos(alpha-theta2)))/(18*cos(theta1-theta2)*cos(theta1-theta2)-32)
    return diff(PotentialEnergy, theta2)
def Phi1(theta1, theta2):
    return math.asin(((-648*sin(alpha - theta1)*cos(alpha - theta1) + 324*sin(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - \
        324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))/(18*cos(theta1 - theta2)**2 - 32) + 36*(-324*cos(alpha - theta1)**2 + \
        324*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 144*cos(alpha - theta2)**2)*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2)/(-6*omega**2))
def Phi2(theta1, theta2):
    return math.asin(((324*sin(alpha - theta2)*cos(alpha - theta1)*cos(theta1 - theta2) - 288*sin(alpha - theta2)*cos(alpha - theta2) + \
       324*sin(theta1 - theta2)*cos(alpha - theta1)*cos(alpha - theta2))/(18*cos(theta1 - theta2)**2 - 32) - 36*(-324*cos(alpha - theta1)**2 + \
      324*cos(alpha - theta1)*cos(alpha - theta2)*cos(theta1 - theta2) - 144*cos(alpha - theta2)**2)*sin(theta1 - theta2)*cos(theta1 - theta2)/(18*cos(theta1 - theta2)**2 - 32)**2)/(-2*omega**2))

def equations(p):
    theta1, theta2 = p
    return (DPE_theta1(), DPE_theta2())

def delta_max(theta1_now , theta1_last, theta2_now, theta2_last):
    deltamax = 1
    delta1 = math.fabs(theta1_now - theta1_last)
    delta2 = math.fabs(theta2_now - theta2_last)
    if k>0:
        deltamax = max(delta1, delta2)
    return deltamax

theta1_now, theta2_now = (math.pi, math.pi)
theta1_last, theta2_last = (0, 0)
k = 0
delta = 0.000001
omega = 2.5
alpha = math.pi/4
Theta1 =[]
Theta2 = []
while omega > 0:
    while delta_max(theta1_now , theta1_last, theta2_now, theta2_last)>= delta:
         theta1_last = theta1_now
         theta2_last = theta2_now
         theta1_now = Phi1(theta1_last, theta2_last)
         theta2_now = Phi2(theta1_now,theta2_last)
         k=k+1
    omega = omega - 0.1
    if omega == 2:
        print("eee")
    k=0
    Theta1.append(theta1_now)
    Theta2.append(theta2_now)
plt.axis([-math.pi/2,2*math.pi/3,-math.pi/2,2*math.pi/3])
plt.plot(Theta1,Theta2,'crimson')
plt.grid(True)
plt.show()
   