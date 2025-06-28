import numpy.random as npr
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
X = npr.standard_normal((5000)) # N(0,1)
Y = npr.normal(-1, 1, (5000))   # N(mu, var)
Z = npr.uniform(-3, 3, (5000))  # U(a,b)
W = npr.lognormal(0, 1, (5000)) 

plt.hist(W, bins = 50, cumulative=False, color='g')
plt.title('lognormal LN(0,1)', fontsize = 18)
plt.xlabel('RV', color='k', fontsize=12)
#hist(Y, bins = 50, cumulative=False, color='k')