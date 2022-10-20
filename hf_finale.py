from __future__ import division
from inspect import CO_NESTED
from math import *
from socket import CAN_ISOTP;import numpy as np
from timeit import default_timer as timer
import multiprocessing
###############################################################################################################
#                                                                                                             #
#                                           Hartree-Fock method                                               #
#                                                                                                             #
###############################################################################################################
start_time=timer()
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

S=overlap();T=kinetic();V=potential()

def twoint(i,j,k,l):
    try:
        twoint=(2/sqrt(pi))*sqrt(((exponents[i]+exponents[j])*(exponents[k]+exponents[l]))/(exponents[i]+exponents[j]+exponents[k]+exponents[l]))*S[i][j]*S[k][l]
        return twoint
    except: return 0

def build_new_fock(D,nbas,H):
    F=np.zeros((nbas,nbas))
    for i in range(nbas):
        for j in range(nbas):
            F[i,j]=H[i,j]
            for k in range(nbas):
                for l in range(nbas):
                    F[i,j]+=D[k,l]*(2*twoint(i,j,k,l)-twoint(i,k,j,l))
    return F

def diagonalization(F,s_inv_sqrt,H):
    Fao=np.matmul(s_inv_sqrt,np.matmul(F,s_inv_sqrt))
    epsilon,C_prime=np.linalg.eig(Fao)
    idx = epsilon.argsort()[::1] #There is a routine that orders the eVals increasingly
    epsilon = epsilon[idx]          #
    C_prime = C_prime[:,idx]        #
    C=np.matmul(s_inv_sqrt,C_prime)
    C_occ=C[:,0:(z//2)]
    D=np.matmul(C_occ,C_occ.transpose())
    Eelec=np.trace(np.matmul(D,(H+F)))
    return D,Eelec,Fao,C,epsilon,C_occ

#   Building a new error matrix
def build_erros_vector(F,D,S):
    error_vec=(np.matmul(F,np.matmul(D,S))-np.matmul(S,np.matmul(D,F))).reshape(nbas*nbas,1)
    norm_error_vec=sqrt(np.matmul(error_vec.transpose(),error_vec)[0,0])
    return error_vec, norm_error_vec

#   Building the B matrix and solve for the DIIS coefficient
def get_diis_coeff(errvec_subspace):
    b_mat=np.zeros((len(errvec_subspace)+1,len(errvec_subspace)+1))
    b_mat[-1,:]=-1; b_mat[:,-1]=-1; b_mat[-1,-1]=0
    rhs=np.zeros((len(errvec_subspace)+1,1))
    rhs[-1,-1]=-1
    for i in range(len(errvec_subspace)):
        for j in range(i+1):
            b_mat[i,j]=np.matmul(errvec_subspace[i].transpose(),errvec_subspace[j])
            b_mat[j,i]=b_mat[i,j]
    *diis_coeff, _ =np.linalg.solve(b_mat,rhs)
    return diis_coeff

#   Building an extrapolated fock matrix from a linear combination
def extrapolated_fock(F_subspace,diis_coeff):
    extrapolated_Fmat=np.zeros((nbas,nbas))
    for i in range(len(F_subspace)):
        extrapolated_Fmat+=F_subspace[i]*diis_coeff[i]
    return extrapolated_Fmat

#   A class that stores error vectors, fock and density matrices
class subspace(list):
    def append(self, item):
        list.append(self,item)
        if len(self) > dimsubspace: del self[0]

print("="*58)
print(f"* Output for the SCF energy calculation *")
print("="*58,"\n\n")

#   Global variables
H=T+V
iterations=12
scf_tol=1e-12
density_tol=1e-12
error_tol=1e-12
print("Total number of basis function: {}".format(nbas))
print("SCF energy convergence criterion: {}".format(scf_tol))
print("Density matrix convergence criterion: {}".format(density_tol))
print("DIIS error vector convergence criterion: {}".format(error_tol))

#   S^(-1/2) matrix calculation
s_eigval, s_eigvec=np.linalg.eig(S)
s_halfeig=np.zeros((nbas,nbas))
for i in range(nbas):
    s_halfeig[i,i]=s_eigval[i]**(-0.5)
a=np.matmul(s_eigvec,s_halfeig)
b=s_eigvec.transpose()
s_inv_sqrt=np.matmul(a,b)
#print(s_inv_sqrt)

##########################
########## DIIS ##########
##########################
dimsubspace=8

#   creating instances of the class subspace
errvec_subspace= subspace()
F_subspace= subspace()
D_subspace= subspace()

#   Initializing variables for convergence check
old_energy, old_density=0, np.zeros((nbas,nbas))

#   Initial guess for the Fock matrix
F=np.zeros((nbas,nbas))
for i in range(nbas):
    for j in range(nbas):
        F[i,j]=H[i,j]

print("\nHF SCF ITERATIONS:\n------------------")

#   The DIIS-SCF iterations
for i in range(iterations):
    print(f"Iteration {i+1}:")
    if i <=1:
        D=diagonalization(F,s_inv_sqrt,H)[0]
        D_subspace.append(D)
        energy=diagonalization(F,s_inv_sqrt,H)[1]
        print(f"     * HF electronic energy: {energy} Ha.")
        F=build_new_fock(D,nbas,H)
        F_subspace.append(F)
        errvec, norm_errvec=build_erros_vector(F,D,S)
        errvec_subspace.append(errvec)
        print(f"    * Norm of DIIS error vector = {norm_errvec}")
    else:
        #   build and solve the B matrix
        diis_coeff=get_diis_coeff(errvec_subspace)
        #   build extrapolated fock matrices
        extrapolated_F=extrapolated_fock(F_subspace,diis_coeff)
        #   Density matrix from extrapolated Fock
        D=diagonalization(extrapolated_F,s_inv_sqrt,H)[0]
        D_subspace.append(D)
        #   next entry in the fock subspace obtained from the recent density matrix obtained
        F=build_new_fock(D,nbas,H)
        F_subspace.append(F)
        #  electronic energy
        energy=diagonalization(F,s_inv_sqrt,H)[1]
        print(f"    * HF electronic energy: {energy} Ha.")
        #   error vector calculation
        errvec,norm_errvec=build_erros_vector(F,D,S)
        errvec_subspace.append(errvec)
        print(f"    * Norm of DIIS error vector = {norm_errvec}")
        #   DIIS ends here
    #   Checking for various convergence criteria
    delta_E=energy-old_energy
    print(f"    * Energy convergence = {delta_E} Eh")
    D_vector, old_D_vector=np.reshape(D,(nbas*nbas,1)),np.reshape(old_density,(nbas*nbas,1))
    norm_D=D_vector-old_D_vector
    D_conv=(np.matmul(norm_D.transpose(),norm_D))**0.5
    print(f"    * Norm of density matrix convergence = {D_conv[0,0]}\n")
    #   Breaking the iterations when convergence is reached
    if abs(delta_E) < scf_tol and abs(D_conv[0,0]) < density_tol:
        print(f"SCF ENERGY AND DENSITY CONVERGED WITHIN {i+1} ITERATIONS!\n\n")
        MO_energies=diagonalization(F,s_inv_sqrt,H)[4]
        print(f"ORBITAL ENERGIES:\n-----------------\n{MO_energies}\n\n")
        break
    else:
        old_density=D
        old_energy=energy
MO_energies=diagonalization(F,s_inv_sqrt,H)[4]
print(f"FINAL RESULT:\n-------------")
print(f"     * DIIS-HF electronic energy: {energy} Ha.")
print(f"     * HOMO-LUMO gap = {MO_energies[z//2]-MO_energies[(z-1)//2]} Ha\n")
###############################################################################################################
#                                                                                                             #
#                                           MP2                                                               #
#                                                                                                             #
###############################################################################################################
C=diagonalization(F,s_inv_sqrt,H)[3]
#   Transforming the two electrons itegrals from AO to MO basis
def twoint_transform():
    pqrs=np.zeros((nbas,nbas,nbas,nbas))
    for p in range(nbas):  
        for q in range(nbas):  
            for r in range(nbas):  
                for s in range(nbas):  
                    for mu in range(nbas):  
                        for nu in range(nbas):  
                            for lam in range(nbas):  
                                for sig in range(nbas):  
                                    pqrs[p,q,r,s] += C[mu,p]*C[nu,q]*C[lam,r]*C[sig,s]*twoint(mu,nu,lam,sig)
    return pqrs

b=np.zeros((nbas,nbas,nbas,nbas))
b=twoint_transform()
def twointMO(i,j,k,l):
    a=np.zeros((nbas,nbas,nbas,nbas))
    a[i,j,k,l]=b[i,j,k,l]
    return b[i,j,k,l]
#   Transforming spatial MO to spin MO
ndim=nbas*2     #Each spatial orbital is now a two spin orbitals
spin_int=np.zeros((ndim,ndim,ndim,ndim))
for p in range(ndim):
    for q in range(ndim):
        for r in range(ndim):
            for s in range(ndim):
                prqs=twointMO(p//2,r//2,q//2,s//2)*(p%2==r%2)*(q%2==s%2)
                psqr=twointMO(p//2,s//2,q//2,r//2)*(p%2==s%2)*(q%2==r%2)
                spin_int[p,q,r,s]=prqs-psqr
F_spin=np.zeros((ndim))
for i in range(ndim):
    F_spin[i]=MO_energies[i//2]
#   MP2 iteration
EMP2=0.0
for i in range(z):
    for j in range(z):
        for a in range(z,ndim):
            for b in range(z,ndim):
                EMP2+=0.25*pow(spin_int[i,j,a,b],2)/(F_spin[i]+F_spin[j]-F_spin[a]-F_spin[b])
print(f"     * MP2 correlation energy: {EMP2} Ha.")
print(f"     * MP2 electronic energy: {energy+EMP2} Ha\n.")
###############################################################################################################
#                                                                                                             #
#                         Coupled Cluster Doubles : doi/10.1002/qua.560140503                                 #
#                                                                                                             #
###############################################################################################################
#   Spin basis Fock matrix Eigen valus, Putting MO energies in a diagonal matrix
F_spin=np.diag(F_spin)
#   Initialize the u and v vectors
u2=np.zeros((ndim,ndim,ndim,ndim))
v2=np.zeros((ndim,ndim,ndim,ndim))
delta_abij=np.zeros((ndim,ndim,ndim,ndim))
for a in range(z,ndim):
    for b in range(z,ndim):
        for i in range(z):
            for j in range(z):
                delta_abij[a,b,i,j]=F_spin[a,a]+F_spin[b,b]-F_spin[i,i]-F_spin[j,j]

#   Initialize the "A" vector
a2=np.zeros((ndim,ndim,ndim,ndim))
def make_a2(u2,v2):
    a2=np.zeros((ndim,ndim,ndim,ndim))
    for a in range(z,ndim):
        for b in range(z,ndim):
            for i in range(z):
                for j in range(z):
                    a2[a,b,i,j]=-(1/delta_abij[a,b,i,j])*(spin_int[a,b,i,j]+u2[a,b,i,j]+v2[a,b,i,j])
    return a2
#   Making the linear array
def make_u2(a1,a2):
    for a in range(z,ndim):
        for b in range(z,ndim):
            for i in range(z):
                for j in range(z):
                    for k in range(z):
                        u2[a,b,i,j]=-spin_int[k,b,i,j]*a1[a,k]+spin_int[k,a,i,j]*a1[b,k]
                        for l in range(z):
                            u2[a,b,i,j]+=0.5*spin_int[k,l,i,j]*a2[a,b,k,l]
                        for c in range(z,ndim):
                            u2[a,b,i,j]+=-spin_int[k,b,j,c]*a2[a,c,i,k]-spin_int[k,a,j,c]*a2[c,b,i,k]\
                            -spin_int[k,a,i,c]*a2[a,b,k,j]-spin_int[k,b,i,c]*a2[a,c,k,j]
                    for c in range(z,ndim):
                        u2[a,b,i,j]+=spin_int[a,b,c,j]*a1[c,i]-spin_int[a,b,c,i]*a1[c,j]
                        for d in range(z,ndim):
                            u2[a,b,i,j]+=0.5*spin_int[a,b,c,d]*a2[c,d,i,j]
    return u2
#   Making the quadratic array
def make_v2(a2):
    for i in range(z):
        for j in range(z):
            for a in range(z,ndim):
                for b in range(z,ndim):
                    for k in range(z):
                        for l in range(z):
                            for c in range(z,ndim):
                                for d in range(z,ndim):
                                    v2[i,j,a,b]=0.25*spin_int[k,l,c,d]*(a2[c,d,i,j]*a2[a,b,k,l]-2*(a2[a,c,i,j]*a2[b,d,k,l] + a2[b,d,i,j]*a2[a,c,k,l])\
                                        -2*(a2[a,b,i,k]*a2[c,d,j,l] + a2[c,d,i,k]*a2[a,b,j,l])+4*(a2[a,c,i,k]*a2[b,d,j,l] + a2[b,d,i,k]*a2[a,c,j,l]))
    return v2

#   ECCD energy
def ccdenergy(a2):
    ECCD=0.0
    for i in range(z):
        for j in range(z):
            for a in range(z,ndim):
                for b in range(z,ndim):
                    ECCD+=0.25*spin_int[i,j,a,b]*a2[a,b,i,j]
    return ECCD

#   Iterative process
ECCD = 0
DECCD = 1.0
while DECCD > 1e-12: # arbitrary convergence criteria
    OLDECCD = ECCD
    a2=make_a2(u2,v2)
    ECCD=ccdenergy(a2)
    #print(ECCD)
    DECCD = abs(ECCD - OLDECCD)
    u2=make_u2(np.zeros((ndim,ndim)),a2)
    v2=make_v2(a2)
print(f"     * CCD correlation energy: {ECCD} Ha.")
print(f"     * CCD electronic energy: {energy+ECCD} Ha\n.")
###############################################################################################################
#                                                                                                             #
#                                           Coupled Cluster Single Double                                     #
#                                                                                                             #
###############################################################################################################
#   Initializing and T1 and T2 arrays
t1=np.zeros((ndim,ndim))
t2=np.zeros((ndim,ndim,ndim,ndim))
#   Initial guess for the t2 from previous MP2 calculations
for a in range(z,ndim):
    for b in range(z,ndim):
        for i in range(z):
            for j in range(z):
                t2[a,b,i,j]+=spin_int[i,j,a,b]/(F_spin[i,i]+F_spin[j,j]-F_spin[a,a]-F_spin[b,b])

#   Denominator arrays Dai Daibij
Dai=np.zeros((ndim,ndim))
for a in range(z,ndim):
    for i in range(z):
        Dai[a,i]=F_spin[i,i]-F_spin[a,a]
Dabij=np.zeros((ndim,ndim,ndim,ndim))
for a in range(z,ndim):
    for b in range(z,ndim):
        for i in range(z):
            for j in range(z):
                Dabij[a,b,i,j]=F_spin[i,i]+F_spin[j,j]-F_spin[a,a]-F_spin[b,b]
#   Make the effective two particle excitation amplitudes
def tau_tilde(a,b,i,j):
    tau_tilde=t2[a,b,i,j]+0.5*(t1[a,i]*t1[b,j]-t1[b,i]*t1[a,j])
    return tau_tilde

def Tau(a,b,i,j):
    Tau=t2[a,b,i,j]+t1[a,i]*t1[b,j]-t1[b,i]*t1[a,j]
    return Tau

#   Define the updated F and W intermediates, the Kronecker delta will be used
def updateintermediates():
        Fae=np.zeros((ndim,ndim))
        for a in range(z,ndim):
            for e in range(z,ndim):
                Fae[a,e]=(1-(a==e))*F_spin[a,e]
                for m in range(z):
                    Fae[a,e]+=-0.5*F_spin[m,e]*t1[a,m]
                    for f in range(z,ndim):
                        Fae[a,e]+=t1[f,m]*spin_int[m,a,f,e]
                    for n in range(z):
                        for f in range(z,ndim):
                            Fae[a,e]+=-0.5*tau_tilde(a,f,m,n)*spin_int[m,n,e,f]
        
        Fmi=np.zeros((ndim,ndim))
        for m in range(z):
            for i in range(z):
                Fmi[m,i]=(1-(m==i))*F_spin[m,i]
                for e in range(z,ndim):
                    Fmi[m,i]+=0.5*t1[e,i]*F_spin[m,e]
                    for n in range(z):
                        Fmi[m,i]+=t1[e,n]*spin_int[m,n,i,e]
                for n in range(z):
                    for e in range(z,ndim):
                        for f in range(z,ndim):
                            Fmi[m,i]+=0.5*tau_tilde(e,f,i,n)*spin_int[m,n,e,f]
        
        Fme=np.zeros((ndim,ndim))
        for m in range(z):
            for e in range(z,ndim):
                Fme[m,e]=F_spin[m,e]
                for n in range(z):
                    for f in range(z,ndim):
                        Fme[m,e]+=t1[f,n]*spin_int[m,n,e,f]

        Wmnij=np.zeros((ndim,ndim,ndim,ndim))
        for m in range(z):
            for n in range(z):
                for i in range(z):
                    for j in range(z):
                        Wmnij[m,n,i,j]=spin_int[m,n,i,j]
                        for e in range(z,ndim):
                            Wmnij[m,n,i,j]+=t1[e,j]*spin_int[m,n,i,e]-t1[e,i]*spin_int[m,n,j,e]
                            for f in range(z,ndim):
                                Wmnij[m,n,i,j]+=0.25*Tau(e,f,i,j)*spin_int[m,n,e,f]

        Wabef=np.zeros((ndim,ndim,ndim,ndim))
        for a in range(z,ndim):
            for b in range(z,ndim):
                for e in range(z,ndim):
                    for f in range(z,ndim):
                        Wabef[a,b,e,f]=spin_int[a,b,e,f]
                        for m in range(z):
                            Wabef[a,b,e,f]+=-t1[b,m]*spin_int[a,m,e,f]+t1[a,m]*spin_int[b,m,e,f]
                            for n in range(z):
                                Wabef[a,b,e,f]+=0.25*Tau(a,b,m,n)*spin_int[m,n,e,f]

        Wmbej=np.zeros((ndim,ndim,ndim,ndim))
        for m in range(z):
            for b in range(z,ndim):
                for e in range(z,ndim):
                    for j in range(z):
                        Wmbej[m,b,e,j]=spin_int[m,b,e,j]
                        for f in range(z,ndim):
                            Wmbej[m,b,e,j]+=t1[f,j]*spin_int[m,b,e,f]
                        for n in range(z):
                            Wmbej[m,b,e,j]+=-t1[b,n]*spin_int[m,n,e,j]
                            for f in range(z,ndim):
                                Wmbej[m,b,e,j]+=-(0.5*t2[f,b,j,n] + t1[f,j]*t1[b,n])*spin_int[m,n,e,f]
        return Fae, Fmi, Fme, Wmnij, Wabef, Wmbej

#   Make T1 and T2 amplitudes
def makeT1(t1,t2):
        t1new=np.zeros((ndim,ndim))
        for a in range(z,ndim):
            for i in range(z):
                t1new[a,i]=F_spin[i,a]
                for e in range(z,ndim):
                    t1new[a,i]+=t1[e,i]*Fae[a,e]
                for m in range(z):
                    t1new[a,i]+=-t1[a,m]*Fmi[m,i]
                    for e in range(z,ndim):
                        t1new[a,i]+=t2[a,e,i,m]*Fme[m,e]
                        for f in range(z,ndim):
                            t1new[a,i]+=-0.5*t2[e,f,i,m]*spin_int[m,a,e,f]
                        for n in range(z):
                            t1new[a,i]+=-0.5*t2[a,e,m,n]*spin_int[n,m,e,i]
                for n in range(z):
                    for f in range(z,ndim):
                        t1new[a,i]+=-t1[f,n]*spin_int[n,a,i,f]
                t1new[a,i]=t1new[a,i]/Dai[a,i]
        return t1new

def makeT2(t1,t2):
        t2new = np.zeros((ndim,ndim,ndim,ndim))
        for a in range(z,ndim):
            for b in range(z,ndim):
                for i in range(z):
                    for j in range(z):
                        t2new[a,b,i,j] = spin_int[i,j,a,b]
                        for e in range(z,ndim):
                            t2new[a,b,i,j] += t2[a,e,i,j]*Fae[b,e] - t2[b,e,i,j]*Fae[a,e]
                            for m in range(z):
                                t2new[a,b,i,j] += -0.5*t2[a,e,i,j]*t1[b,m]*Fme[m,e] + 0.5*t2[b,e,i,j]*t1[a,m]*Fme[m,e]
                        for m in range(z):
                            t2new[a,b,i,j] += -t2[a,b,i,m]*Fmi[m,j] + t2[a,b,j,m]*Fmi[m,i]
                            for e in range(z,ndim):
                                t2new[a,b,i,j] += -0.5*t2[a,b,i,m]*t1[e,j]*Fme[m,e] + 0.5*t2[a,b,j,m]*t1[e,i]*Fme[m,e]
                        for e in range(z,ndim):
                            t2new[a,b,i,j] += t1[e,i]*spin_int[a,b,e,j] - t1[e,j]*spin_int[a,b,e,i]
                            for f in range(z,ndim):
                                t2new[a,b,i,j] += 0.5*Tau(e,f,i,j)*Wabef[a,b,e,f]
                        for m in range(z):
                            t2new[a,b,i,j] += -t1[a,m]*spin_int[m,b,i,j] + t1[b,m]*spin_int[m,a,i,j]  
                            for e in range(z,ndim):
                                t2new[a,b,i,j] +=  t2[a,e,i,m]*Wmbej[m,b,e,j] - t1[e,i]*t1[a,m]*spin_int[m,b,e,j]
                                t2new[a,b,i,j] += -t2[a,e,j,m]*Wmbej[m,b,e,i] + t1[e,j]*t1[a,m]*spin_int[m,b,e,i]
                                t2new[a,b,i,j] += -t2[b,e,i,m]*Wmbej[m,a,e,j] + t1[e,i]*t1[b,m]*spin_int[m,a,e,j]
                                t2new[a,b,i,j] +=  t2[b,e,j,m]*Wmbej[m,a,e,i] - t1[e,j]*t1[b,m]*spin_int[m,a,e,i]
                            for n in range(z):
                                t2new[a,b,i,j] += 0.5*Tau(a,b,m,n)*Wmnij[m,n,i,j]
                        t2new[a,b,i,j] = t2new[a,b,i,j]/Dabij[a,b,i,j] 
        return t2new
#   CCSD Energy
def ccsdenergy():
    ECCSD = 0.0
    for i in range(z):
        for a in range(z,ndim):
            ECCSD+=F_spin[i,a]*t1[a,i]
        for j in range(z):
            for a in range(z,ndim):
                for b in range(z,ndim):
                    ECCSD+=0.25*spin_int[i,j,a,b]*t2[a,b,i,j] + 0.5*spin_int[i,j,a,b]*(t1[a,i]*t1[b,j])
    return ECCSD

#   CCSD iterations
ECCSD = 0
DECCSD = 1.0
while DECCSD > 1e-12: # arbitrary convergence criteria
    OLDECCSD = ECCSD
    Fae,Fmi,Fme,Wmnij,Wabef,Wmbej = updateintermediates()
    ECCSD = ccsdenergy()
    DECCSD = abs(ECCSD - OLDECCSD)
    t1 = makeT1(t1,t2)
    t2 = makeT2(t1,t2)
print(f"     * CCSD correlation energy: {ECCSD} Ha.")
print(f"     * CCSD electronic energy: {energy+ECCSD} Ha\n.")
###############################################################################################################
#                                                                                                             #
#             Quadratic configuration interaction Single Double, doi: 10.1063/1.453520                        #
#                                                                                                             #
###############################################################################################################
a1=np.zeros((ndim,ndim))
a2=np.zeros((ndim,ndim,ndim,ndim))
w1=np.zeros((ndim,ndim))
v1=np.zeros((ndim,ndim))
w2=np.zeros((ndim,ndim,ndim,ndim))
v2=np.zeros((ndim,ndim,ndim,ndim))
#   The equation providing the amplitude a1
def make_a1(w1,v1):
    for a in range(z,ndim):
        for i in range(z):
            a1[a,i]=(w1[i,a]+v1[i,a])/Dai[a,i]
    return a1
#   The equation providing the amplitude a2
def make_a2(w2,v2):
    for a in range(z,ndim):
        for b in range(z,ndim):
            for i in range(z):
                for j in range(z):
                    a2[a,b,i,j]=(spin_int[a,b,i,j]+w2[i,j,a,b]+v2[i,j,a,b])/Dabij[a,b,i,j]
    return a2
#   Computing the w1 term
def make_w1(a1,a2):
    for i in range(z):
        for a in range(z,ndim):
            for j in range(z):
                for b in range(z,ndim):
                    w1[i,a]=-spin_int[j,a,i,b]*a1[b,j]
                    for c in range(z,ndim):
                        w1[i,a]+=-0.5*spin_int[j,a,b,c]*a2[b,c,i,j]
                for k in range(z):
                    for b in range(z,ndim):
                        w1[i,a]+=-0.5*spin_int[j,k,i,b]*a2[a,b,j,k]
    return w1
#   Compute the w2 term
def make_w2(a1,a2):
    for i in range(z):
        for j in range(z):
            for a in range(z,ndim):
                for b in range(z,ndim):
                    for c in range(z,ndim):
                        w2[i,j,a,b]=spin_int[a,b,c,j]*a1[c,i]-spin_int[a,b,c,i]*a1[c,j]
                        for d in range(z,ndim):
                            w2[i,j,a,b]+=0.5*spin_int[a,b,c,d]*a2[c,d,i,j]
                    for k in range(z):
                        w2[i,j,a,b]+=-spin_int[k,b,i,j]*a1[a,k]+spin_int[k,a,i,j]*a1[b,k]
                        for l in range(z):
                            w2[i,j,a,b]+=0.5*spin_int[k,l,i,j]*a2[a,b,k,l]
                        for c in range(z,ndim):
                            w2[i,j,a,b]+=-spin_int[k,b,j,c]*a2[a,c,i,k]-spin_int[k,a,j,c]*a2[c,b,i,k]-spin_int[k,b,i,c]*a2[a,c,k,j]-spin_int[k,a,i,c]*a2[c,b,k,j]
    return w2
#   Compute the quadratic array v1
def make_v1(a1,a2):
    for i in range(z):
        for a in range(z,ndim):
            for j in range(z):
                for k in range(z):
                    for b in range(z,ndim):
                        for c in range(z,ndim):
                            v1[i,a]=0.5*spin_int[j,k,b,c]*(a1[b,i]*a2[c,a,j,k] + a1[a,j]*a2[c,b,i,k] + 2*a1[b,j]*a2[a,c,i,k])
    return v1
#   The quadratic v2 array is fetched from the CCD block
#   The QCISD   Energy
def qcisd_energy(a2):
    qcisd=0.0
    for i in range(z):
        for j in range(z):
            for a in range(z,ndim):
                for b in range(z,ndim):
                    qcisd+=0.25*spin_int[a,b,i,j]*a2[a,b,i,j]
    return qcisd
#   QCISD Iterations
EQCISD = 0
DEQCISD = 1.0
while DEQCISD > 1e-12: # arbitrary convergence criteria
    OLDEQCISD = EQCISD
    a1=make_a1(w1,v1)
    a2=make_a2(w2,v2)
    EQCISD=qcisd_energy(a2)
    v1=make_v1(a1,a2)
    v2=make_v2(a2)
    w1=make_w1(a1,a2)
    w2=make_w2(a1,a2)
    DEQCISD= abs(EQCISD - OLDEQCISD)
print(f"     * QCISD correlation energy: {EQCISD} Ha.")
print(f"     * QCISD electronic energy: {energy+EQCISD} Ha\n.")
end_time = timer()
print(f"Wall time = {end_time-start_time} seconds")
#print("NUmber of cpu: ", multiprocessing.cpu_count())