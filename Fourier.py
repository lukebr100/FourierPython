import math
import time
import matplotlib.pyplot as plt 
import sympy as sym
import numpy as np
# an O(n^2) algorithm for Fourier Coefficients
# a method called Fast Fouier Transform can do it 
# in O(n*log(n))


# inputs: 
# f piecewise continuous function
# n number of terms
# length of half interval L
#
#ouput:
# finite Fouier series of f on [-L, L]
def fourier(f, n, L):
    x = sym.Symbol('x') # declare symolic variable x using sympy
    pi = math.pi # store math.pi into variable named pi
    LI = 1 / L # the inverse of L
    an = [] # initialize list x2
    bn = []
    a0 = 1 / (2 * L) * sym.integrate(f, (x, -L, L)) # calculate a0 term
    an = an + [a0] # concatenate a0 to an
    bn = bn + [0] # b0 term is always 0, do this so that a loop later one is easier
    # next we calculate the desired number of an terms
    for i in range(1, n):
        c = sym.cos(i * pi * LI * x ) # store symbolic function cosine into c
        t = sym.integrate(LI* f * c, (x, -L, L)) # integrate LI * c * f wrt x over [-L, L]
        s = t.evalf() # evaluate t so we have a float instead of symbolic expression
        an.append(s) # append this float s to lits of an's
    for i in range(1, n):
      # do the same thing as before, but with sine, bn instead of cosine, an
        s = sym.sin(i * pi * x / L)
        t = sym.integrate(LI* f * s, (x, -L, L))
        s = t.evalf()
        bn.append(s)
    F = 0 # initialize F
    for i in range(1, n - 1):
        F = F + an[i] * sym.cos(i * pi * LI * x) + bn[i] * sym.sin(i * pi * LI * x)
    return a0 + F
