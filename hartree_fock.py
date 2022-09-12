#Hartree-Fock method
import math
import numpy as np
from scipy.linalg import eig, inv
z=int(input('Enter the Atomic NUmber of your atom: '))
list1=[]
if z==2:
    file = open('basis_He.txt')
    for line in file:
        fields = line.strip().split()
        list1.append(float(fields[0]))
    #print(list1)
elif z==4:
    file = open('basis_Be.txt')
    for line in file:
        fields = line.strip().split()
        list1.append(float(fields[0]))
    #print(list1)
else: print('For the time being i only process Helium and Berylium.')
R, C=int(len(list1)), int(len(list1))
# Initialize the Overlap matrix 
t = [] 
# For user input
S_matrix=[] 
for i in list1[0:]:        # A for loop for row entries 
    a =[] 
    for j in list1[0:]:     # A for loop for column entries 
        b=math.sqrt(math.pow(2*(math.sqrt(i*j)/(i+j)),3))
        a.append(b)
    S_matrix.append(a)
# For printing the T matrix 
#for i in S_matrix:
    #print(i)
#Initialize the T matrix
T_matrix=[] 
for i in range(0,int(len(list1)),1):
    a=[]
    for j in range(0,int(len(list1)),1):
        b=((3*list1[i]*list1[j])/(list1[i]+list1[j]))*S_matrix[i][j]
        a.append(b)
    T_matrix.append(a)
print("-------------------Kinetic Energy Matrix-------------------")
for i in T_matrix:
  print(i)
#Initialize the Vmatrix
V_matrix=[] 
for i in range(0,int(len(list1)),1):
    a=[]
    for j in range(0,int(len(list1)),1):
        b=-(2/(math.sqrt(math.pi)))*z*(math.sqrt(list1[i]+list1[j]))*S_matrix[i][j]
        a.append(b)
    V_matrix.append(a)
print("-------------------Potential Energy Matrix-------------------")
for i in V_matrix:
  print(i)
#Initialize the H matrix
H_matrix=[] 
for i in range(0,int(len(list1)),1):
    a=[]
    for j in range(0,int(len(list1)),1):
        b=T_matrix[i][j]+V_matrix[i][j]
        a.append(b)
    H_matrix.append(a)
print("-------------------Hamiltonian Matrix-------------------")
for i in H_matrix:
  print(i)
eVals,eVecs = eig(S_matrix)    
#print(eVals)
#print(eVecs)
D = np.zeros((len(list1),len(list1)))
for i in range(0,len(eVals)):
    D[i,i] = eVals[i].real
Sinv = inv(eVecs)
print("-------------------EigenVectors of the S matrix-------------------")
print(eVecs)
print("-------------------EigenValues of the S matrix--------------------------")
print(D)
print("-------------------Inverse of the S matrix--------------------------")
print(Sinv)
fao=H_matrix    #Starting with the bar Hamiltonian matrix
for i in range(11):
    print('Iteration number {}'.format(i))
    