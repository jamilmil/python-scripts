#Finding the zero point of a function
#Interval method
a=float(input('Enter the lowest boundary: '))
b=float(input('Enter the greatest boundary: '))
import math
def interval_method(a,b):
    list1=[]
    y1=math.sin(a)
    y2=math.sin(b)
    list1.append(y1)
    list1.append(y2)
    c=0
    while c != math.pi:
        c=(a+b)/2
        y3=math.sin(c)
        list1.append(y3)
        if (list1[0]*list1[-1] < 0): b,c=c,b
        else: a,c=c,a
    return c
print('The Interval Method gives {0}={1}.'.format('\u03C0',interval_method(a,b)))
#Regula Falsi
def regula_falsi(a,b):
    list2=[]
    c=0
    y1=math.sin(a)
    y2=math.sin(b)
    list2.append(y1)
    list2.append(y2)
    while (list2[1]-list2[0]) !=0:
        if abs(list2[0]) > abs(list2[1]):
            a,b=b,a
            list2[0],list2[1]=list2[1],list2[0]
        c=(a*list2[1]-(b*list2[0]))/(list2[1]-list2[0])
        y3=math.sin(c)
        list2.append(y3)
        if abs(list2[-1]) < abs(list2[1]):
            b=c
            list2[1],list2[-1]=list2[-1],list2[1]
            #print('a-b=', a-b)
        else:
            print('fail')
    return c
print('The Regula Falsi Method gives {0}={1}.'.format('\u03C0',regula_falsi(a,b)))
#Newton method
def newton(a,b):
    list3=[]
    c=0
    while b != math.pi:
        list3.clear()
        y1=math.sin(a)
        y11=math.cos(a)
        list3.insert(0,y1)
        list3.append(y11)
        b=a-(list3[0]/list3[1])
        y2=math.sin(b)
        list3.append(y2)
        if abs(list3[-1]) < abs(y1):
            a=b
    return b
print("Newton's Method gives {0}={1}.".format('\u03C0',newton(a,b)))