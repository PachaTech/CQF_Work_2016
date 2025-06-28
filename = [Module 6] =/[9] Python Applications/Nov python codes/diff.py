# Numerical differentiation and root finding
from math import *
from matplotlib import *
from pylab import *
    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Q1: function of x^{x/3} 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def func(x):
    y = x ** (x / 3)
    return y
    
a = func(1)
b = func(2)
print 'f(a) = ', func(a)
print 'f(b) = ', func(b)


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Q2: derivative of function x^{x/3} 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def func_diff(x, diff_step, CASE):
    if CASE == 'forward':
        y = (func(x + diff_step) - func(x)) / diff_step
    elif CASE == 'backward':
        y = (func(x) - func(x - diff_step)) / diff_step
    elif CASE == 'central':
        y = (func(x + diff_step) - func(x - diff_step)) \
                / diff_step / 2
    else:
        print('No numerical differentiation rule found')
    return y
    
def real_func_diff(x):
    y = (log(x) + 1) * x**(x/3) / 3
    return y

d = 0.0000001
a_diff = func_diff(1, d, 'forward')
b_diff = func_diff(2, d, 'central')
a_real_diff = real_func_diff(1)
b_real_diff = real_func_diff(2)
print 'forward derivative ', (a_diff)
print 'central derivative ', (b_diff)
print(a_real_diff)
print(b_real_diff)


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Aux: Graph of x^{x/3}
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
import matplotlib.pyplot as plt
import numpy as np
plt.figure(figsize = (4, 3))
x = np.linspace(0, 5, 100)
plt.plot(x, func(x) - 2, 'r-')
plt.plot(x, np.zeros_like(x), 'k-')
show()

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Q3: bisection solve x^{x/3} = 2
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def bisec_func(value, low_bnd, up_bnd, xtol = 0.00001):
    low = low_bnd * 1
    up = up_bnd * 1
    if (func(low)-value) * (func(up)-value) > 0:
        print('Input interval with different output signs')
    else:
        middle = (low + up) / 2
        while((up - low) / 2 > xtol):
            equ_mid = func(middle) - value
            equ_low = func(low) - value
            if equ_mid == 0:
                break
            if equ_low * equ_mid < 0:
                up = middle
            else:
                low = middle
            middle = (up + low) / 2
        return middle
    
bisec2 = bisec_func(2, 1, 3)
print(bisec2)
print(func(bisec2)-2)
            

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Q4: Newton-Raphson solve x^{x/3} = 2
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def newton_func(value, x0, xtol = 0.01, Nmax = 100):
    xn = x0
    prev = 0
    N = 0
    while(xn - prev > xtol and N < Nmax):
        N += 1
        xn = xn - (func(xn) - 2)/func_diff(xn, \
                0.0000001,'forward')
    return xn
    
NR2 = newton_func(2, 2)
print(NR2)
print(func(NR2)-2)
