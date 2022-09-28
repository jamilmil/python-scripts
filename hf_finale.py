from math import *;import numpy as np;import numpy.linalg as linalg
exponents=np.loadtxt('data.dat');z=2;nbas=exponents.size

def overlap():
    S=np.zeros((nbas,nbas))
    for i in range(nbas):
        for j in range(nbas):
            S[i][j]=pow((2*(sqrt(exponents[i]*exponents[j]))/(exponents[i]+exponents[j])),1.5)
    return S

def kinetic():
    T=np.zeros((nbas,nbas))
    for i in range(nbas):
        for j in range(nbas):
            T[i][j]=3*(exponents[i]*exponents[j]/(exponents[i]+exponents[j]))*S[i][j]
    return T

def potential():
    V=np.zeros((nbas,nbas))
    for i in range(nbas):
        for j in range(nbas):
            V[i][j]=(-2/sqrt(pi))*z*(sqrt(exponents[i]+exponents[j]))*S[i][j]
    return V

S=overlap();T=kinetic();V=potential();H=T+V

def twoint(i,j,k,l):
    try:
        twoint=(2/sqrt(pi))*sqrt(((exponents[i]+exponents[j])*(exponents[k]+exponents[l]))/(exponents[i]+exponents[j]+exponents[k]+exponents[l]))*S[i][j]*S[k][l]
        return twoint
    except: return 0

def diagonalization(matrix):     
    eVals,eVecs=linalg.eigh(matrix)
    idx = eVals.argsort()[::1] #There is a routine that orders the eVals increasingly
    eVals = eVals[idx]          #
    eVecs = eVecs[:,idx]        #
    matrix_diag = np.zeros((nbas,nbas))
    for i in range(0,len(eVals)): matrix_diag[i,i] = eVals[i].real
    return matrix_diag, eVecs, eVals
S_diag=diagonalization(S)[0];eVecs=diagonalization(S)[1]

def inverse_sqrt(matrix,eVecs):
    a=np.zeros((nbas,nbas))
    a=np.linalg.inv(np.sqrt(matrix))
    matrix_inv_sqrt=np.matmul(np.matmul(eVecs,a),np.transpose(eVecs))
    return matrix_inv_sqrt
S_inv_sqrt=inverse_sqrt(S_diag,eVecs)

def density(C,D):
    D_old=np.zeros((nbas,nbas))
    for i in range(nbas):
        for j in range(nbas):
            D_old[i,j]=D[i,j]
            D[i,j]=0.0e0
            for k in range(int(z/2)):
                D[i][j]=D[i][j]+2*(C[i][k]*C[j][k])
    return D,D_old

def make_fock_matrix(H,D):
    F=np.zeros((nbas,nbas))
    for i in range(nbas):
        for j in range(nbas):
            F[i,j]=H[i,j]
            for k in range(nbas):
                for l in range(nbas):
                    F[i][j]=F[i][j]+(twoint(i,j,k,l)-0.5*twoint(i,l,k,j))*D[k][l]
    return F

def fock_prime_matrix(S_inv_sqrt,F):
    F_prime=np.matmul(np.transpose(S_inv_sqrt),np.matmul(F,S_inv_sqrt))
    return F_prime

def energy(H,F,D):
    E0=0.0e0
    for i in range(nbas):
        for j in range(nbas):
            E0=E0+(H[i,j]+F[i,j])*D[i,j]*0.5
    return E0

def delta_D(D,D_old):
    delta=0.0e0
    for i in range(nbas):
        for j in range(nbas):
            delta=delta+((D[i,j]-D_old[i,j])**2)
    delta=(delta/4)**0.5
    return delta
D=np.zeros((nbas,nbas));delta=1.0
while delta > 0.000000000001:
    F=make_fock_matrix(H,D)
    #   Transforming the Fock Matrix to the Orthonormal Basis
    F_prime=fock_prime_matrix(S_inv_sqrt,F)
    #   Diagonalizing the Fock matrix
    C_prime=diagonalization(F_prime)[1]
    #   Constructing the new SCF eigenvector matrix
    C=np.matmul(S_inv_sqrt,C_prime)
    #   Forming the new Density Matrix
    D,D_old=density(C,D)
    #   Calculate the electronic energy
    E=energy(H,F,D)
    #   Test for convergence of the Density
    delta=delta_D(D,D_old)
print("SCF Energy:",E)
#######################################################################################
#                                                                                     #
#                                     MP2 Energy                                      #
#                                                                                     #
#######################################################################################
#   Initializing the orbital energies epsilon and the two electron integrals
epsilon=list(diagonalization(F_prime)[2])
#   Transforming the two electrons itegrals from AO to MO basis
def twoint_mo():
    pqrs=np.zeros((nbas,nbas,nbas,nbas))
    for p in range(0,nbas):  
        for q in range(0,nbas):  
            for r in range(0,nbas):  
                for s in range(0,nbas):  
                    for mu in range(0,nbas):  
                        for nu in range(0,nbas):  
                            for lam in range(0,nbas):  
                                for sig in range(0,nbas):  
                                    pqrs[p,q,r,s] += C[mu,p]*C[nu,q]*C[lam,r]*C[sig,s]*twoint(mu,nu,lam,sig)
    return pqrs
pqrs=twoint_mo()
#   Transforming spatial MO to spin MO
nbas=nbas*2
integrals=np.zeros((nbas,nbas,nbas,nbas))
try:
    for p in range(1,nbas+1):
        for q in range(1,nbas+1):
            for r in range(1,nbas+1):
                for s in range(1,nbas+1):
                    prqs=pqrs[(p+1)//2][(r+1)//2][(q+1)//2][(s+1)//2]*(p%2==r%2)*(q%2==s%2)
                    psqr=pqrs[(p+1)//2][(s+1)//2][(q+1)//2][(r+1)//2]*(p%2==s%2)*(q%2==r%2)
                    integrals[p-1,q-1,r-1,s-1]=prqs-psqr
except: IndexError: 0
#   Spin basis fock matrix eigenvalues
F_spin=np.zeros((nbas))
for i in range(nbas):
    F_spin[i]=epsilon[i//2]
#   MP2 iteration
EMP2=0.0
for i in range(z):
    for j in range(z):
        for a in range(z,nbas):
            for b in range(z,nbas):
                EMP2+=0.25*pow(integrals[i,b,j,a],2)/(F_spin[i]+F_spin[j]-F_spin[a]-F_spin[b])
print("EMP2 correlation Energy:",EMP2)