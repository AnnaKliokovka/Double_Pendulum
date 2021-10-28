# -*- coding: utf-8 -*-
"""
author: chris manucredo, chris.manucredo@gmail.com
about:
    this is an implementaion of newton's method for solving systems of 
    non-linear equations. as an example we use two different non-linear
    systems. at the end you'll see a work-accuracy plot for both systems.
"""
import DoublePendulumPackage as DP
import DoublePendulumPackage.SimpleLines as SL
import DoublePendulumPackage.FindThetasOmega as FTO
import DoublePendulumPackage.InLines as IL
import numpy as np
import random
from tkinter import *
import math
#A = int(input("Введите угол альфа: "))
#alpha = A*math.pi/180
i=0
while i < 5:
    A = random.randint(0,90)
    alpha = A*math.pi/180
    print('Вычисляется угол %i' % A)
#FTO.calculate(alpha, alpha - 0.5 *math.pi, alpha + 0.5 *math.pi, 1.6 )
#FTO.calculate(alpha, -0.5, 2.27, 1.7 )
#SL.calculate(alpha )
    IL.calculate(alpha)
    i=i+1
