import numpy as np
import numpy.random as npr
from scipy.stats import norm
 
S = 100; K = 100; T = 1
r = 0.05; vol = 0.5
 
I = 10000 # MC paths
Z = npr.standard_normal(I)
ST = S * np.exp((r - 0.5 * vol**2) \
        * T + vol * np.sqrt(T) * Z)
V = np.mean(np.exp(-r * T) \
    * np.maximum(ST - K, 0))
print(V)
 
 
std = np.std(np.exp(-r * T) \
    * np.maximum(ST - K, 0)) / np.sqrt(I)
print(std)
 
d1 = (np.log(S/K) + (r + 0.5 * vol**2) \
        * T) / vol / np.sqrt(T)
d2 = (np.log(S/K) + (r - 0.5 * vol**2) \
        * T) / vol / np.sqrt(T)
BS = S * norm.cdf(d1) - np.exp(-r * T) \
        * K * norm.cdf(d2)
print(BS)