import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
 
S = 100
T = 1.0
r = 0.05
vol = 0.2
 
I = 100 # MC paths
N = 252
dt=T/N
Z = npr.standard_normal((N+1, I))
St = S * np.ones((N+1, I))
rnd = np.exp((r - 0.5 * vol**2) * dt + vol * np.sqrt(dt) * Z)
for i in range(1, N):
    St[i] = St[i-1] * rnd[i-1]
 
plt.figure(figsize = (8, 6))
for k in range(I):
    plt.plot(np.arange(N+1), St[:, k])
plt.xlim(0,252)
plt.title('Monte Carlo path for a stock', fontsize = 18)