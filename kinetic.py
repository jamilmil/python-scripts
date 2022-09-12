#Complex kinetic system
from operator import delitem
import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
def plotting(filename):
    plt.figure()
    data = np.loadtxt(filename, delimiter=',', unpack=True)

    x = data[0]
    y1 = data[1] # the first graph
    y2 = data[2] # the second graph
    y3 = data[3] # the third graph
    plt.plot(x,y1, label='Compound 1')
    plt.plot(x,y2, label='Compound 2')
    plt.plot(x,y3, label='Compound 3')
    plt.xlabel('Time(s)')
    plt.ylabel('Concentration(mol/L)')
    plt.legend()
    plt.title(filename)
    plt.show()
#Euler's method
dt=float(input("Enter the step size: "))
k1,k2,k3,k4=0.2,0.1,0.1,0.02
conc_1=[1]
conc_2=[0]
conc_3=[0]
conc_1_prime=[1-(k1*dt)/2]
conc_2_prime=[k1]
conc_3_prime=[0]
conc_1_rk2=[]
conc_2_rk2=[]
conc_3_rk2=[]
conc_1_rk4=[]
conc_2_rk4=[]
conc_3_rk4=[]
time=[t for t in np.arange(0, 50+dt, dt)]
f_conc_1=[-k1]
f_conc_2=[k1]
f_conc_3=[0]
f_conc_1_prime=[-k1*conc_1_prime[0]+k2*conc_2_prime[0]]
f_conc_2_prime=[k1*conc_1_prime[0]-k2*conc_2_prime[0]-k3*conc_3_prime[0]+k4*conc_3_prime[0]]
f_conc_3_prime=[k3*conc_2_prime[0]-k4*conc_3_prime[0]]
conc_1_prime2=[conc_1[0]+f_conc_1_prime[0]*(dt/2)]
conc_2_prime2=[conc_2[0]+f_conc_2_prime[0]*(dt/2)]
conc_3_prime2=[conc_3[0]+f_conc_3_prime[0]*(dt/2)]
f_conc_1_prime2=[-k1*conc_1_prime2[0]+k2*conc_2_prime2[0]]
f_conc_2_prime2=[k1*conc_1_prime2[0]-k2*conc_2_prime2[0]-k3*conc_3_prime2[0]+k4*conc_3_prime2[0]]
f_conc_3_prime2=[k3*conc_2_prime2[0]-k4*conc_3_prime2[0]]
conc_1_prime3=[conc_1[0]+f_conc_1_prime2[0]*dt]
conc_2_prime3=[conc_2[0]+f_conc_2_prime2[0]*dt]
conc_3_prime3=[conc_3[0]+f_conc_3_prime2[0]*dt]
f_conc_1_prime3=[-k1*conc_1_prime3[0]+k2*conc_2_prime3[0]]
f_conc_2_prime3=[k1*conc_1_prime3[0]-k2*conc_2_prime3[0]-k3*conc_3_prime3[0]+k4*conc_3_prime3[0]]
f_conc_3_prime3=[k3*conc_2_prime3[0]-k4*conc_3_prime3[0]]
for i in range(0,len(time)):
    conc_1.append(conc_1[i]+f_conc_1[i]*dt)
    f_conc_1.append(-k1*conc_1[i]+k2*conc_2[i])
    conc_2.append(conc_2[i]+f_conc_2[i]*dt)
    f_conc_2.append(k1*conc_1[i]-k2*conc_2[i]-k3*conc_2[i]+k4*conc_3[i])
    conc_3.append(conc_3[i]+f_conc_3[i]*dt)
    f_conc_3.append(k3*conc_2[i]-k4*conc_3[i])
with open("Euler.dat", "w") as f:
    for i in range(len(time)):
        data=[time[i], conc_1[i], conc_2[i], conc_3[i]]
        writer=csv.writer(f)
        writer.writerow(data)
#Runge-Kutta 2nd order
for i in range(0,len(time)):
    conc_1_prime.append(conc_1[i]+f_conc_1[i]*(dt/2))
    f_conc_1_prime.append(-k1*conc_1_prime[i]+k2*conc_2_prime[i])
    conc_1_rk2.append(conc_1[i]+f_conc_1_prime[i]*dt)
    conc_2_prime.append(conc_2[i]+f_conc_2[i]*(dt/2))
    f_conc_2_prime.append(k1*conc_1_prime[i]-k2*conc_2_prime[i]-k3*conc_2_prime[i]+k4*conc_3_prime[i])
    conc_2_rk2.append(conc_2[i]+f_conc_2_prime[i]*(dt/2))
    conc_3_prime.append(conc_3[i]+f_conc_3[i]*(dt/2))
    f_conc_3_prime.append(k3*conc_2_prime[i]-k4*conc_3_prime[i])
    conc_3_rk2.append(conc_3[i]+f_conc_3_prime[i]*(dt/2))
with open("RK2.dat", "w") as f:
    for i in range(len(time)):
        data=[time[i], conc_1_rk2[i], conc_2_rk2[i], conc_3_rk2[i]]
        writer=csv.writer(f)
        writer.writerow(data)
#Runge-Kutta 4th order
for i in range(0,len(time)):
    conc_1_prime2.append(conc_1[i]+f_conc_1_prime[i]*(dt/2))
    f_conc_1_prime2.append(k1*conc_1_prime2[i]+k2*conc_2_prime2[i])
    conc_1_prime3.append(conc_1[i]+f_conc_1_prime2[i]*dt)
    f_conc_1_prime3.append(k1*conc_1_prime3[i]+k2*conc_2_prime3[i])
    conc_1_rk4.append(conc_1[i]+(f_conc_1[i]+2*f_conc_1_prime[i]+2*f_conc_1_prime2[i]+f_conc_1_prime3[i])*(dt/6))
    conc_2_prime2.append(conc_2[i]+f_conc_2_prime[i]*(dt/2))
    f_conc_2_prime2.append(k1*conc_1_prime2[i]-k2*conc_2_prime2[i]-k3*conc_2_prime2[i]+k4*conc_3_prime2[i])
    conc_2_prime3.append(conc_2[i]+f_conc_2_prime2[i]*dt)
    f_conc_2_prime3.append(k1*conc_1_prime3[i]-k2*conc_2_prime3[i]-k3*conc_2_prime3[i]+k4*conc_3_prime3[i])
    conc_2_rk4.append(conc_2[i]+(f_conc_2[i]+2*f_conc_2_prime[i]+2*f_conc_2_prime2[i]+f_conc_2_prime3[i])*(dt/6))
    conc_3_prime2.append(conc_3[i]+f_conc_3_prime[i]*(dt/2))
    f_conc_3_prime2.append(k3*conc_2_prime2[i]-k4*conc_3_prime2[i])
    conc_3_prime3.append(conc_3[i]+f_conc_3_prime2[i]*dt)
    f_conc_3_prime3.append(k3*conc_2_prime3[i]-k4*conc_3_prime3[i])
    conc_3_rk4.append(conc_3[i]+(f_conc_3[i]+2*f_conc_3_prime[i]+2*f_conc_3_prime2[i]+f_conc_3_prime3[i])*(dt/6))
with open("RK4.dat", "w") as f:
    for i in range(len(time)):
        data=[time[i], conc_1_rk4[i], conc_2_rk4[i], conc_3_rk4[i]]
        writer=csv.writer(f)
        writer.writerow(data)
for x in ["Euler.dat", "RK2.dat", 'RK4.dat']:
    plotting(x)
