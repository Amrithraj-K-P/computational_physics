import numpy as np

def midpoint_integ(f,lim,n): #lim : [a,b]
    a,b=lim
    h=abs((b-a)/n)
    x_i=a #dummy variable
    out=0
    while x_i<=b:
        out += f((x_i + x_i + h)/2)
        x_i += h
    return out*h


def trapezoidal_integ(f,lim,n): #lim : [a,b]
    a,b=lim
    h=abs((b-a)/n)
    out = f(a) + f(b)
    x_i = a+h
    while x_i<b: #summation is done over terms with weight = 2
        out += 2*f(x_i)
        x_i+=h
    return out*h/2