# -*- coding: utf-8 -*-
"""
author: chris manucredo, chris.manucredo@gmail.com
about:
    this is an implementaion of newton's method for solving systems of 
    non-linear equations. as an example we use two different non-linear
    systems. at the end you'll see a work-accuracy plot for both systems.
"""
import DoublePendulumPackage.InLines as IL
import numpy as np
import random
from tkinter import *
import math
#A = int(input("Введите угол альфа: "))
A = 32
step = 1
while A < 90 :
    A = A + step
    alpha = A*math.pi/180
    IL.calculate(alpha)
    print('Посчитан угол %i' % A)
print('Вычисления выполнены!')
