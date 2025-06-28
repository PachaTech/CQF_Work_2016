import numpy as np
from scipy import integrate
 
def fun(x):
    return np.log(x)
     
value, error = integrate.quad(fun,0,1)
print(value)
print(error)