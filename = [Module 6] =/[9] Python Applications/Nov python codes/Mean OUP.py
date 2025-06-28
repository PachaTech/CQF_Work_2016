import numpy as np
import numpy.random as npr
 
def Process(t, Z, x0 = 1, theta = 1, \
            mu = 1, sigma = 0.5):
    y1 = x0 * np.exp(-theta * t)  
    y2 = mu * (1 - np.exp(-theta * t))
    y3 = sigma * np.sqrt((1 - np.exp(-theta * t)) \
            / 2 / theta) * Z
    return y1 + y2 + y3
 
# code begin:
Z1 = npr.standard_normal((10000))
myProc = Process(10, Z1)
myMean = np.mean(myProc)
print('The mean is %f' %myMean)