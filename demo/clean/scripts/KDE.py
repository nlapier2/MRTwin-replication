import numpy as np
import math
from scipy.interpolate import interp1d
from scipy.stats import norm

def BinDist(x, weights, lo, up, n):
    y = np.array([0.] * (n * 2))
    x_min = 0
    x_max = n - 2
    x_del = (up - lo) / (n-1)
    for i in range(len(x)):
        if np.isfinite(x[i]):
            xpos = (x[i] - lo)/x_del
            if not np.isfinite(xpos):
                continue
            ix = int(math.floor(xpos))
            fx = xpos - ix
            wi = weights[i]
            if x_min <= ix & ix <= x_max:
                y[ix] += (1 - fx) * wi
                y[ix + 1] += fx * wi
            elif ix == -1:
                y[0] += fx*wi
            elif ix == (x_max + 1):
                y[ix] += (1 - fx) * wi
    return y


def KDE(x, bw, weights): 
    n = max(512, len(x))
    if n > 512:
        n = 2**(math.ceil(math.log(n, 2)))
    from_ = min(x) - 3*bw/4
    to_ = max(x) + 3*bw/4
    lo = from_ - 4*bw/4
    up = to_ + 4*bw/4
    y = BinDist(x, weights, lo, up, n)
    kords = np.linspace(0, 2 * (up - lo), int(2 * n))
    kords[(n+1):(2*n)] = -1 * (kords[1:n])[::-1]
    kords = norm.pdf(kords, scale=bw/4)
    kords = np.fft.ifft(np.fft.fft(y) * np.conjugate(np.fft.fft(kords)))/ len(x)
    kords = np.maximum(kords.real[0:n]/len(y), 0 * len(kords))
    #print(kords[0])
    xords = np.linspace(lo, up, n)
    #print(xords[0])
    x_a = np.linspace(from_, to_, n)
    f = interp1d(xords, kords)
    y_o = f(x_a)
    return x_a, y_o

