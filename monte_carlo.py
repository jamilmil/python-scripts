#Monte Carlo Integration
n1, n2, n3, n4=input('Enter the 1st number of random points you want: ').split()
list_n=[int(n1), int(n2), int(n3), int(n4)]
import math
import random
from tabulate import tabulate
def monte_carlo(b):
    print('1) <I> and \u0394<I> Calculations for k=10:')
    data =[]
    pin_list=[]
    for a in b:
        list= []
        int_list= []
        for k in range(1,11,1):
            list.clear()
            for i in range(1,a+1,1):
                x=random.uniform(0,1)
                y=random.uniform(0,1)
                y_func=math.exp(-math.pow(x,2))
                if (y <= y_func):
                    list.append(y)
            number_of_elements=int(len(list))
            integral=number_of_elements/a
            int_list.append(integral)
        sum1=sum(int_list)
        average_int=sum1/10
        squared_numbers = [number ** 2 for number in int_list]
        sum2=sum(squared_numbers)
        delta_I=((sum2)-(math.pow(sum1,2))/10)/90
        #create data
        data.append([a*10, delta_I])
        pin_list.append([a,(2/math.sqrt(math.pi))*average_int])
        print('For n={0}: <I> ={1} and \u0394<I> ={2}.'.format(a, average_int, delta_I))
    print()
    print("2) Convergence Investigation:")
    #define headers name
    col_name = ["n.k", "\u0394<I>"]
    print(tabulate(data, headers=col_name))
    print()
    print('3) Pin value:')
    col_name = ["n", "Pin"]
    print(tabulate(pin_list, headers=col_name))
    print()
    print('4) More efficient estimiations using an infinite sum:')
    sum3=0
    for i in range(0,170,1):
        I2=(math.pow(-1,i))/((2*i+1)*math.factorial(i))
        sum3=sum3+I2
    print('<I> ={0}'.format(sum3))
    return
print(monte_carlo(list_n))