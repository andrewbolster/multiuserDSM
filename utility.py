# This is all for experiments at the minute and not being called by anything... yet

__author__="bolster"
__date__ ="$02-Dec-2010 18:48:38$"

import math, cmath

material={
        "r_0c"	:179,
        "a_c"	:35.89e-3,
        "r_0s"	:0,
        "a_s"	:0,

        "l_0"	:0.695e-3,
        "l_inf"	:585e-6,

        "b"		:1.2,
        "f_m"	:1000,

        "c_inf"	:55e-9,
        "c_0"	:1e-9,
        "c_e"	:0.1,

        "g_0"	:0.5e-9,
        "g_e"	:1.033
        }

def do_transfer_function(length,freq,type=2):#TODO Where on earth to Zs/Zl come from? They do nothing!
    """
    Z=impedance/l, Y=admittance/l
    Z=R+jwL, Y=G+jwC
    Z0=Characteristic Impedence, gamma=propagation constant
    Z0=sqrt(Z/Y), gamma=sqrt(Z*Y)
    Should Use PyGSL, but start off with cmath
    """

    w = 2 * math.pi * freq

    Z = complex(_R(freq),w*_L(freq))
    Y = complex(_G(freq),w*_C(freq))

    Z0 = cmath.sqrt(Z/Y)
    gamma = cmath.sqrt(Z*Y)

    gammad=gamma*length

    Zs = complex(100,0)
    Zl = complex(100,0)
    print("w",w)
    print("R",_R(freq),"L",_L(freq),"G",_G(freq),"C",_C(freq))
    print("Z",Z,"Y",Y,"Gamma",gamma,"Gammad",gammad)
    upper = Z0 * ( 1/cmath.cosh(gammad)) #sech=1/cosh
    lower = Zs * ( (Z0/Zl) + cmath.tanh(gammad) ) + Z0 * ( 1+ (Z0/Zl)*cmath.tanh(gammad) )

    H = upper/lower
    return cmath.polar(H)[0]        #polar returns (abs(h),real(h))

def _R(freq):
    """
    Return R Parameter for transfer function
    """
    c_partial = math.pow(material["a_c"]*math.pow(freq,2)+math.pow(material["r_0c"],4),(0.25))

    if material["r_0s"] > 0:
        s_partial = math.pow(material["a_s"]*math.pow(freq,2)+math.pow(material["r_0s"],4),(0.25))
        return (c_partial*s_partial)/(c_partial+s_partial)
    else:
        return c_partial

    return

def _L(freq):
    """
    Return L Parameter for transfer function
    """
    return (material["l_0"]+material["l_inf"]*math.pow(_L_partial(freq),material["b"]))/(1+_L_partial(freq))

def _L_partial(freq):
    return freq*1e-3/material["f_m"]

def _C(freq):
    return material["c_inf"]+material["c_0"]*math.pow(freq,-material["c_e"])

def _G(freq):
    return material["g_0"]*math.pow(freq,material["g_e"])

if __name__ == "__main__":

    print do_transfer_function(1000,8.5e4)
