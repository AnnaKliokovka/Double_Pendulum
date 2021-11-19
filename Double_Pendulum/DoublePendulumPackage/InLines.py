import DoublePendulumPackage.FindThetasOmega as FTO
import DoublePendulumPackage.CheckStability as CS
import DoublePendulumPackage.AlphaOmega as AO
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from sympy import diff, symbols, cos, sin
import math
# Глобальные переменные - угол вибрации, частота вибрации
# Global values - 
alpha=0 
# Глобальные переменные 
OMEGA = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
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

def calculate_30(alpha_in):
    global X, Y, OMEGA, alpha, omega
    X = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    Y = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    OMEGA = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    alpha = alpha_in
    deltaOmega = 0.005
    i_now =0
    for i in [-1,0,1,2]:
        for j in [2,1,0,-1]: 
            omega = 0.001
            x0 = np.array([alpha + i*0.5*math.pi,alpha + j*0.5*math.pi])
            while omega <= 10:
                    root = opt.fsolve(f, x0, full_output = 1) 
                    x0_now =np.array([root[0][0],root[0][1]])
                    x_check = np.array([[x0_now[0] - x0[0]],[x0_now[1] - x0[1]]])
                    if max(x_check) < 0.4 and root[2] == 1 :
                            X[i_now] = np.append(X[i_now], x0_now[0])
                            Y[i_now] = np.append(Y[i_now], x0_now[1])
                            OMEGA[i_now] = np.append(OMEGA[i_now],omega)
                            x0 = x0_now
                            omega = round(omega + deltaOmega, 4)
                    else:
                        x_star = FTO.calculate(alpha, X[i_now][len(X[i_now])-1], Y[i_now][len(Y[i_now])-1],omega)
                        if (x_star[0] == X[i_now][len(X[i_now])-1] and x_star[1] == Y[i_now][len(Y[i_now])-1] and max(x_check) <= 0.55):
                            X[i_now] = np.append(X[i_now], x0_now[0])
                            Y[i_now] = np.append(Y[i_now], x0_now[1])
                        else:
                            X[i_now] = np.append(X[i_now], x_star[0])
                            Y[i_now] = np.append(Y[i_now], x_star[1])
                            OMEGA[i_now] = np.append(OMEGA[i_now],omega)
                            break   
            if (i_now ==4 ):
                x0_now =np.array([X[i_now][len(X[i_now])-1],Y[i_now][len(Y[i_now])-1]])
                #X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                #Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                #OMEGA[i_now+1] = np.append(OMEGA[i_now+1],omega)
                x_not =np.array([X[i_now][len(X[i_now])-2],Y[i_now][len(Y[i_now])-2]])
                delta_x = 0.001
                delta_y = 0
                x0 =np.array([X[i_now][len(X[i_now])-1] +delta_x,Y[i_now][len(Y[i_now])-1]-delta_y])
                omega = round(omega - deltaOmega, 4)
                b = 1
                while omega > 0:
                    root = opt.fsolve(f, x0, full_output = 1) 
                    x0_now =np.array([root[0][0],root[0][1]])
                    x_check = np.array([[abs(x0_now[0] - x_not[0])],[abs(x0_now[1] - x_not[1])]])
                    k=0
                    while max(x_check)<0.0001 and k< 500 and b==1:
                        x0 =np.array([x0[0] +delta_x,x0[1]-delta_y])
                        root = opt.fsolve(f, x0, full_output = 1) 
                        x0_now =np.array([root[0][0],root[0][1]])
                        x_check = np.array([[abs(x0_now[0] - x_not[0])],[abs(x0_now[1] - x_not[1])]])
                        k=k+1
                    if k==500:
                        break
                    else:
                        b=0
                        x_check = np.array([[x0_now[0] - x0[0]],[x0_now[1] - x0[1]]])
                        if max(x_check) < 0.4  and root[2] == 1 :
                                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1],omega)
                                x0 = x0_now
                                omega = round(omega - deltaOmega, 4)
                        else:
                            x_star = FTO.calculate(alpha, X[i_now][len(X[i_now])-1], Y[i_now][len(Y[i_now])-1],omega)
                            if (x_star[0] == X[i_now][len(X[i_now])-1] and x_star[1] == Y[i_now][len(Y[i_now])-1] and max(x_check) <= 0.55):
                                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                            else:
                                X[i_now+1] = np.append(X[i_now+1], x_star[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x_star[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1],omega)
                                break  
                i_now = i_now+1
            #if (i_now ==4 and alpha == 33*math.pi/180):
                #i_now = i_now+1
            if (i_now ==9):
                x0_now =np.array([X[i_now][len(X[i_now])-1],Y[i_now][len(Y[i_now])-1]])
                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                OMEGA[i_now+1] = np.append(OMEGA[i_now+1], omega)
                x_not =np.array([X[i_now][len(X[i_now])-2],Y[i_now][len(Y[i_now])-2]])
                delta_x = 0.0001
                delta_y = 0
                x0 =np.array([X[i_now][len(X[i_now])-1] -delta_x,Y[i_now][len(Y[i_now])-1]-delta_y])
                omega = round(omega - deltaOmega, 4)
                b = 1
                while omega > 0:
                    root = opt.fsolve(f, x0, full_output = 1) 
                    x0_now =np.array([root[0][0],root[0][1]])
                    x_check = np.array([[abs(x0_now[0] - x_not[0])],[abs(x0_now[1] - x_not[1])]])
                    k=0
                    while max(x_check)<0.0001 and k< 500 and b==1:
                        x0 =np.array([x0[0] -delta_x,x0[1]-delta_y])
                        root = opt.fsolve(f, x0, full_output = 1) 
                        x0_now =np.array([root[0][0],root[0][1]])
                        x_check = np.array([[abs(x0_now[0] - x_not[0])],[abs(x0_now[1] - x_not[1])]])
                        k=k+1
                    if k==500:
                        break
                    else:
                        b = 0
                        x_check = np.array([[abs(x0_now[0] - x0[0])],[abs(x0_now[1] - x0[1])]])
                        if max(x_check) < 0.4  and root[2] == 1 :
                                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1], omega)
                                x0 = x0_now
                                omega = round(omega - deltaOmega, 4)
                        else:
                            x_star = FTO.calculate(alpha, X[i_now+1][len(X[i_now+1])-1], Y[i_now+1][len(Y[i_now+1])-1],omega)
                            if (x_star[0] == X[i_now+1][len(X[i_now+1])-1] and x_star[1] == Y[i_now+1][len(Y[i_now+1])-1] and max(x_check) <= 0.55):
                                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1], omega)
                            else:
                                X[i_now+1] = np.append(X[i_now+1], x_star[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x_star[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1], omega)
                                break  
                i_now = i_now+1
            i_now = i_now+1  
    if alpha==33*math.pi/180:
       plot_result(i_now )
   # else:
    #plot_result(i_now)
    return AO.OmegaAlpha(alpha, OMEGA)

def calculate_70(alpha_in):
    if alpha_in >= 71*math.pi/180 and alpha_in <73*math.pi/180:
        index_line_1=9
        index_line_2 = 11
    elif alpha_in >= 69*math.pi/180 and alpha_in <= 70*math.pi/180:
        index_line_1 = 10
        index_line_2=12
    else:
        index_line_1 = 4
        index_line_2= 9
    global X, Y, OMEGA, alpha, omega
    X = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    Y = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    OMEGA = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    alpha = alpha_in
    deltaOmega = 0.005
    i_now =0
    for i in [-1,0,1,2]:
        for j in [2,1,0,-1]:
            omega = 0.001
            #if omega == 1.036:
               #print("catchit")
            x0 = np.array([alpha + i*0.5*math.pi,alpha + j*0.5*math.pi])
            while omega <= 10:
                    root = opt.fsolve(f, x0, full_output = 1) 
                    x0_now =np.array([root[0][0],root[0][1]])
                    x_check = np.array([[x0_now[0] - x0[0]],[x0_now[1] - x0[1]]])
                    if max(x_check) < 0.4 and root[2] == 1 :
                            X[i_now] = np.append(X[i_now], x0_now[0])
                            Y[i_now] = np.append(Y[i_now], x0_now[1])
                            OMEGA[i_now] = np.append(OMEGA[i_now],omega)
                            x0 = x0_now
                            omega = round(omega + deltaOmega, 4)
                    else:
                        x_star = FTO.calculate(alpha, X[i_now][len(X[i_now])-1], Y[i_now][len(Y[i_now])-1],omega)
                        if (x_star[0] == X[i_now][len(X[i_now])-1] and x_star[1] == Y[i_now][len(Y[i_now])-1] and max(x_check) <= 0.55):
                            X[i_now] = np.append(X[i_now], x0_now[0])
                            Y[i_now] = np.append(Y[i_now], x0_now[1])
                            OMEGA[i_now] = np.append(OMEGA[i_now],omega)
                        else:
                            X[i_now] = np.append(X[i_now], x_star[0])
                            Y[i_now] = np.append(Y[i_now], x_star[1])
                            OMEGA[i_now] = np.append(OMEGA[i_now],omega)
                            break   
            if (i_now ==index_line_1):
                x0_now =np.array([X[i_now][len(X[i_now])-1],Y[i_now][len(Y[i_now])-1]])
                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                OMEGA[i_now+1] = np.append(OMEGA[i_now+1],omega)
                x_not =np.array([X[i_now][len(X[i_now])-2],Y[i_now][len(Y[i_now])-2]])
                delta_x = 0.001
                delta_y = 0
                x0 =np.array([X[i_now][len(X[i_now])-1] +delta_x,Y[i_now][len(Y[i_now])-1]-delta_y])
                omega = round(omega - deltaOmega, 4)
                b = 1
                while omega > 0:
                    root = opt.fsolve(f, x0, full_output = 1) 
                    x0_now =np.array([root[0][0],root[0][1]])
                    x_check = np.array([[abs(x0_now[0] - x_not[0])],[abs(x0_now[1] - x_not[1])]])
                    k=0
                    while max(x_check)<0.0001 and k< 500 and b==1:
                        x0 =np.array([x0[0] +delta_x,x0[1]-delta_y])
                        root = opt.fsolve(f, x0, full_output = 1) 
                        x0_now =np.array([root[0][0],root[0][1]])
                        x_check = np.array([[abs(x0_now[0] - x_not[0])],[abs(x0_now[1] - x_not[1])]])
                        k=k+1
                    if k==500:
                        break
                    else:
                        b=0
                        x_check = np.array([[x0_now[0] - x0[0]],[x0_now[1] - x0[1]]])
                        if max(x_check) < 0.4  and root[2] == 1 :
                                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1],omega)
                                x0 = x0_now
                                omega = round(omega - deltaOmega, 4)
                        else:
                            x_star = FTO.calculate(alpha, X[i_now+1][len(X[i_now+1])-1], Y[i_now+1][len(Y[i_now+1])-1],omega)
                            if (x_star[0] == X[i_now+1][len(X[i_now+1])-1] and x_star[1] == Y[i_now+1][len(Y[i_now+1])-1] and max(x_check) <= 0.55):
                                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1],omega)
                            else:
                                X[i_now+1] = np.append(X[i_now+1], x_star[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x_star[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1],omega)
                                break  
                i_now = i_now+1
            if (i_now ==index_line_2):
                x0_now =np.array([X[i_now][len(X[i_now])-1],Y[i_now][len(Y[i_now])-1]])
                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                OMEGA[i_now+1] = np.append(OMEGA[i_now+1], omega)
                x_not =np.array([X[i_now][len(X[i_now])-2],Y[i_now][len(Y[i_now])-2]])
                delta_x = 0.0001
                delta_y = 0
                x0 =np.array([X[i_now][len(X[i_now])-1] -delta_x,Y[i_now][len(Y[i_now])-1]-delta_y])
                omega = round(omega - deltaOmega, 4)
                b = 1
                while omega > 0:
                    root = opt.fsolve(f, x0, full_output = 1) 
                    x0_now =np.array([root[0][0],root[0][1]])
                    x_check = np.array([[abs(x0_now[0] - x_not[0])],[abs(x0_now[1] - x_not[1])]])
                    k=0
                    while max(x_check)<0.0001 and k< 500 and b==1:
                        x0 =np.array([x0[0] -delta_x,x0[1]-delta_y])
                        root = opt.fsolve(f, x0, full_output = 1) 
                        x0_now =np.array([root[0][0],root[0][1]])
                        x_check = np.array([[abs(x0_now[0] - x_not[0])],[abs(x0_now[1] - x_not[1])]])
                        k=k+1
                    if k==500:
                        break
                    else:
                        b = 0
                        x_check = np.array([[abs(x0_now[0] - x0[0])],[abs(x0_now[1] - x0[1])]])
                        if max(x_check) < 0.4  and root[2] == 1 :
                                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1], omega)
                                x0 = x0_now
                                omega = round(omega - deltaOmega, 4)
                        else:
                            x_star = FTO.calculate(alpha, X[i_now+1][len(X[i_now+1])-1], Y[i_now+1][len(Y[i_now+1])-1],omega)
                            if (x_star[0] == X[i_now+1][len(X[i_now+1])-1] and x_star[1] == Y[i_now+1][len(Y[i_now+1])-1] and max(x_check) <= 0.55):
                                X[i_now+1] = np.append(X[i_now+1], x0_now[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x0_now[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1], omega)
                            else:
                                X[i_now+1] = np.append(X[i_now+1], x_star[0])
                                Y[i_now+1] = np.append(Y[i_now+1], x_star[1])
                                OMEGA[i_now+1] = np.append(OMEGA[i_now+1], omega)
                                break  
                i_now = i_now+1
            i_now = i_now+1    
    plot_result(i_now)
    return AO.OmegaAlpha(alpha, OMEGA)

def calculate(alpha_in):
    if alpha_in < 33*math.pi/180 and alpha_in > 23*math.pi/180 or alpha_in  > 62*math.pi/180  and alpha_in  < 67*math.pi/180  :
        return calculate_30(alpha_in)
    elif alpha_in>= 69*math.pi/180  and alpha_in< 73*math.pi/180:
        return calculate_70(alpha_in)
    else:
        global X, Y, OMEGA,alpha, omega
        X = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        Y = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        OMEGA = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        alpha = alpha_in
        deltaOmega = 0.005
        delta_x = 0.0001
        delta_y=0.0001
        i_now =0
        count_i=0
        for i in [-1,0,1,2]:
            for j in [2,1,0,-1]:
                omega = 0.001
                x0 = np.array([alpha + i*0.5*math.pi,alpha + j*0.5*math.pi])
                while omega <= 10:
                      if len(X[0]) ==138:
                          print("")
                      root = opt.fsolve(f, x0, full_output = 1) 
                      x0_now =np.array([root[0][0],root[0][1]])
                      x_check = np.array([[abs(x0_now[0] - x0[0])],[abs(x0_now[1] - x0[1])]])
                      if max(x_check) < 0.4 and root[2] == 1:
                             X[i_now] = np.append(X[i_now], x0_now[0])
                             Y[i_now] = np.append(Y[i_now], x0_now[1])
                             OMEGA[i_now] = np.append(OMEGA[i_now], omega)
                             x0 = x0_now
                             omega = round(omega + deltaOmega, 4)
                      else:
                          x_star = FTO.calculate(alpha, X[i_now][len(X[i_now])-1], Y[i_now][len(Y[i_now])-1],omega)
                          if (x_star[0] == X[i_now][len(X[i_now])-1] and x_star[1] == Y[i_now][len(Y[i_now])-1] and max(x_check) <= 0.55):
                              X[i_now] = np.append(X[i_now], x0_now[0])
                              Y[i_now] = np.append(Y[i_now], x0_now[1])
                          else:
                              X[i_now] = np.append(X[i_now], x_star[0])
                              Y[i_now] = np.append(Y[i_now], x_star[1])
                              OMEGA[i_now] = np.append(OMEGA[i_now], omega)
                              break        
                i_now = i_now+1        
        plot_result(i_now)
        return AO.OmegaAlpha(alpha, OMEGA)

def plot_result(CountLines):
    plt.figure(figsize=(8,6))
    plt.axis([-math.pi/2 - 0.05 ,3*math.pi/2 + 0.2,-math.pi/2 - 0.2 ,3*math.pi/2 + 0.2])
    for i in range(CountLines):
        count_lines = i+1
        plt.text(X[i][0] - 0.2,Y[i][0]+0.1,'%i' % count_lines)
        plt.plot(X[i],Y[i],'k',linewidth='0.8')  
        #for j in range(len(X[i])-1):
            #if CS.stability(alpha, X[i][j], Y[i][j], OMEGA[i][j]) == 1:
               # plt.plot([X[i][j],X[i][j+1]],[Y[i][j],Y[i][j+1]],'r',linewidth='0.8') 
            #else:
               # plt.plot([X[i][j],X[i][j+1]],[Y[i][j],Y[i][j+1]],'k',linewidth='0.8') 
        plt.plot(X[i][0],Y[i][0],'ko',markerfacecolor='w')
        plt.plot(X[i][len(X[i])-1],Y[i][len(Y[i])-1],'kD',markerfacecolor='w')
    plt.grid(True,linestyle = '--')
    plt.xlabel(r"$\theta_1$")
    plt.ylabel(r"$\theta_2$")
    plt.xticks([alpha - 0.5*math.pi,0,alpha, alpha+ 0.5*math.pi, math.pi, alpha+math.pi]) 
    plt.yticks([alpha - 0.5*math.pi,0,alpha, alpha+ 0.5*math.pi, math.pi, alpha+math.pi])
    a = alpha*180/math.pi
    plt.title(r'Положения равновесия при $\alpha$ = %f ' % a)    

    plt.savefig('%f.png' % a)



