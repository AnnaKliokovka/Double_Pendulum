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
import matplotlib.pyplot as plt
#A = int(input("Введите угол альфа: "))
A = 89
Classs = []
step = 100
while A < 90 :    
    alpha = A*math.pi/180   
    Classs=np.append(Classs, IL.calculate(alpha) )
    print('Посчитан угол %i' % A)
    A = A + step
print("Строим в плоскости альфа омега")
plt.figure(figsize=(10,6))
plt.axis([0 ,4,0,90])
plt.yticks([5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90]) 
plt.xticks([0,0.5,1,1.5,2,2.5,3,3.5,4])
for i in range(len(Classs)):
    alpha = Classs[i].alpha*180/math.pi
    omegas = Classs[i].omegas
    counts = Classs[i].counts
    for j in range(len(omegas)):
        if counts[j]==4:
            plt.plot(omegas[j],alpha,'o', color = "orange", ms ='1') 
        elif counts[j]==6:
            plt.plot(omegas[j],alpha,'ro',ms ='1') 
        elif counts[j]==8:
            plt.plot(omegas[j],alpha,'go',ms ='1') 
        elif counts[j]==10:
            plt.plot(omegas[j],alpha,'bo',ms ='1')
        elif counts[j]==12:
            plt.plot(omegas[j],alpha,'mo',ms ='1') 
        elif counts[j]==14:
            plt.plot(omegas[j],alpha,'o', color = "aquamarine",ms ='1') 
        elif counts[j]==16:
            plt.plot(omegas[j],alpha,'co',ms ='1') 
        elif counts[j]==18:
            plt.plot(omegas[j],alpha,'ko',ms ='1') 
        else: 
            print('Точек ', counts[j], ', при омеге равной ',omegas[j])
            plt.plot(omegas[j],alpha,'^', color = "yellow",ms ='1') 
plt.title("Зависимость альфа омега")    
plt.savefig("Dependence.png")


           
print('Вычисления выполнены!')
