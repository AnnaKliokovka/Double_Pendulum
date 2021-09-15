import math
from sympy import diff, symbols, cos, sin
theta1, theta2, omega, alpha = symbols('theta1 theta2 omega alpha')
PotentialEnergy = -2*((3*omega*omega)*cos(theta1)+(omega*omega)*cos(theta2))+\
   (4*(81*cos(theta1-theta2)*cos(alpha-theta1)*cos(alpha-theta2)-81*cos(alpha-theta1)*cos(alpha-theta1)- \
   36*cos(alpha-theta2)*cos(alpha-theta2)))/(18*cos(theta1-theta2)*cos(theta1-theta2)-32)
DPE_theta1 = diff(PotentialEnergy, theta1)
DPE_theta2 = diff(PotentialEnergy, theta2)
print(DPE_theta1)
print(DPE_theta2)




