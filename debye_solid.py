#Integration of functions/ Debye Solid
import math
n=int(input('Enter an even maximum number for sampling: '))
a1=float(input('Enter the lowest boundary of the Integral: '))
b1=float(input('Enter the greates boudary of the Integral (1/\u03BD): '))
r=8.31446261815324  #gas constant
#composite simpson's 1/3 rule-----------------------------------------------------------------------------
def simpson(a,b,c): 
    list1=[]    #initializing a list for the abscissa
    list2=[]    #initializing a list for the ordinates
    h=(b-a)/c   #step size
    sum1=0  #initializing the sum for the odd ordinates
    sum2=0  #initializing the sum for even ordinates
    for i in range(0,c+1,1):
        list1.append(a+(i*h)) #Calculating and storing all calculated abscissa
        x=a+(i*h)
        list2.append((math.pow(x,4)*math.exp(x))/(math.pow((math.exp(x)-1),2))) #Calculating and storing all abscissa
    for i in range(1,int(c/2)+1,1):sum1=sum1+list2[(2*i)-1] #looping over and summing the odd ordinates
    for i in range(1,int((c/2)-1)+1,1):sum2=sum2+list2[2*i] #looping over and summing the even ordinates
    integral=(h/3)*(list2[0]+(4*sum1)+(2*sum2)+list2[-1])   #integral calculation
    cv=(9*r*math.pow(1/b,3))*integral   #heat capacity calculation
    return integral, cv
#trapezoidal rule/ Uniform grid--------------------------------------------------------------------------
def trapezium(a,b,c):
    list1=[]
    list2=[]
    h=(b-a)/c   #grid spacing
    #print('h= ', h)
    sum1=0
    list1.append(a)
    list1.append(b)
    list2.append((math.pow(a,4)*math.exp(a))/(math.pow((math.exp(a)-1),2)))
    list2.append((math.pow(b,4)*math.exp(b))/(math.pow((math.exp(b)-1),2)))
    x=a
    for i in range(int(a),int(n)-1):
        x+=h
        list1.append(x)
        list2.append((math.pow(x,4)*math.exp(x))/(math.pow((math.exp(x)-1),2)))
    for i in list2[2:]: sum1=sum1+i
    integral=h*(sum1+(list2[0]+list2[1])/2)
    cv=(9*r*math.pow(1/b,3))*integral   #heat capacity calculation
    return integral, cv
#Gauss-Legendre Quadrature-------------------------------------------------------------------------------
def gaussL(a,b):
    list1=[]    #xi
    list2=[]    #wi
    list3=[]
    list4=[]
    sum1=0
    file = open('values.txt')
    for line in file:
        fields = line.strip().split()
        list1.append(float(fields[0]))
        list2.append(float(fields[1]))
    h=(b-a)/2
    d=(a+b)/2
    for i in list1:list3.append((h*i)+d)
    for i in list3:list4.append((math.pow(i,4)*math.exp(i))/(math.pow((math.exp(i)-1),2)))
    for i in range(4):sum1=sum1+(list2[i]*list4[i])
    integral=h*sum1
    cv=(9*r*math.pow(1/b,3))*integral   #heat capacity calculation
    return integral, cv
print("-The Trapezium rule gives:               I={0} and Cv={1} J.k\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}.mol\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE} for \u03BD={2}.".format(trapezium(a1,b1,n)[0],trapezium(a1,b1,n)[1],float(1/b1)))
print("-The composite Simpson's 1/3 rule gives: I={0} and Cv={1} J.k\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}.mol\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE} for \u03BD={2}.".format(simpson(a1,b1,n)[0],simpson(a1,b1,n)[1],float(1/b1)))
print("-The Gauss-Legendre Quadrature gives:    I={0} and Cv={1} J.k\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}.mol\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE} for \u03BD={2}.".format(gaussL(a1,b1)[0],gaussL(a1,b1)[1],float(1/b1)))