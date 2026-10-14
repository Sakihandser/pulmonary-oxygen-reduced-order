import math
from scipy.integrate import quad
from scipy.optimize import brentq
def sat(p): return p**2.7/(p**2.7+26.8**2.7)
def content(p,h): return 1.34*h*sat(p)+.0031*p
def slope(p,h): return 1.34*h*2.7*sat(p)*(1-sat(p))/p+.0031
def invert(c,h): return brentq(lambda p:content(p,h)-c,1e-8,1000,xtol=1e-11)
def independent_end(mode,q,share,cv,h,d):
    pv=invert(cv,h)
    if mode=='saturation':
        se=sat(100)-(sat(100)-sat(pv))*math.exp(-60*86*share/(1000*q*.3))
        return content(26.8*(se/(1-se))**(1/2.7),h)
    if mode=='content': return content(100,h)-(content(100,h)-cv)*math.exp(-60*86*share/(1000*q*.3))
    # Integrate the separation-of-variables solution, not the production ODE.
    u=d*share/(10*q)
    fun=lambda z:slope(100-(100-pv)*math.exp(-z),h)
    top=max(1, u/.0031+1)
    z=brentq(lambda z:quad(fun,0,z,epsabs=1e-11,epsrel=1e-10)[0]-u,0,top,xtol=1e-10)
    return content(100-(100-pv)*math.exp(-z),h)

