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
from tkinter import *
import math
#SL.calculate(math.pi/4)
#SL.plot_result()
IL.calculate(math.pi/4)