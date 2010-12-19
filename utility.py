# This is all for experiments at the minute and not being called by anything... yet

__author__="bolster"
__date__ ="$02-Dec-2010 18:48:38$"

import math, cmath, numpy

import logging

# Log everything, and send it to stderr.
log = logging.getLogger('multiuserdsm')
log.setLevel(logging.DEBUG)

h = logging.StreamHandler()
f = logging.Formatter('%(levelname)-7s %(module)s %(lineno)d %(message)s')
h.setFormatter(f)
log.addHandler(h)

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

"""
Let the transfer function default to type 2;
    Allows for easy default change later
    Allows for easy 'case-based' changes
"""

def do_transfer_function(length,freq,type=2, measure="m"):#TODO Where on earth to Zs/Zl come from? They do nothing!
    """
    Z=impedance/l, Y=admittance/l
    Z=R+jwL, Y=G+jwC
    Z0=Characteristic Impedence, gamma=propagation constant
    Z0=sqrt(Z/Y), gamma=sqrt(Z*Y)
    Should Use PyGSL, but start off with cmath
    """
    
    if measure == "m":
        length /= 1000
    
    w = 2 * math.pi * freq

    log.debug("Freq=%f,Len=%d",freq,length)
    Z = complex(_R(freq),w*_L(freq))
    Y = complex(_G(freq),w*_C(freq))

    Z0 = cmath.sqrt(Z/Y)
    gamma = cmath.sqrt(Z*Y)

    gammad=gamma*length

    Zs = complex(100,0)
    Zl = complex(100,0)
    log.debug("w:%.3f",w)
    log.debug("R:%.3f L:%.3f G:%.3f C:%.3f",_R(freq),_L(freq),_G(freq),_C(freq))
    log.debug("Z:%s Y:%s Z0:%s Gamma:%s GammaD:%s",complex2str(Z),complex2str(Y),complex2str(Z0),complex2str(gamma),complex2str(gammad))
    upper = Z0 * ( 1/cmath.cosh(gammad)) #sech=1/cosh
    lower = Zs * ( (Z0/Zl) + cmath.tanh(gammad) ) + Z0 * ( 1+ (Z0/Zl)*cmath.tanh(gammad) )

    H = upper/lower
    log.debug("Returned Transfer Function: %g",cmath.polar(H)[0])
    return cmath.polar(H)[0]        #polar returns (abs(h),real(h))

def _R(freq):
    """
    Return R Parameter for transfer function
    """
    c_partial = math.pow(material["a_c"]*freq*freq+math.pow(material["r_0c"],4),(0.25))

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
    return (material["l_0"]+material["l_inf"]*math.pow(freq*1e-3/material["f_m"],material["b"]))/(1+freq*1e-3/material["f_m"])

def _C(freq):
    return material["c_inf"]+material["c_0"]*math.pow(freq,-material["c_e"])

def _G(freq):
    return material["g_0"]*math.pow(freq,material["g_e"])

"""
Mathematical Utilities
"""
def dbmhz_to_watts(psd):
    return UndB(psd)*1e-3

def UndB(input):
    return math.pow(10,input)

def watts_to_dbmhz(psd):
    return TodB(psd*1e3)

def TodB(input):
    return 10*math.log10(input)

def freq_on_tone(K): #TODO
        """
        Assume ADSL downstream for now
        """
        return K * 4312.5 + 140156.25;

def complex2str(complex):
    return '{0:.3f}{1:+.3f}i'.format(complex.real,complex.imag)

if __name__ == "__main__":

    print do_transfer_function(2,8.5e4)
