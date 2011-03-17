# This is all for experiments at the minute and not being called by anything... yet

__author__="bolster"
__date__ ="$02-Dec-2010 18:48:38$"

import math, cmath, numpy, scipy.special as sps
import functools,cPickle

import logging

# Log everything, and send it to stderr.
log = logging.getLogger('multiuserdsm')
log.setLevel(logging.DEBUG)

h = logging.StreamHandler()
f = logging.Formatter('%(levelname)-7s %(module)s %(lineno)d %(message)s')
h.setFormatter(f)
log.addHandler(h)


material=[{ # awg 26
            "r_0c":286.17578,     # ohms/km
            "a_c":0.14769620,

            "l_0":675.36888e-6,   # H/km
            "l_inf":488.95186e-6, # H/km

            "b":0.92930728,
            "f_m":806.33863,      # kHz  

            "c_inf":49e-9,        # F/km
            "c_0":0,
            "c_e":0,

            "g_0":43e-9,          # S/km
            "g_e":0.70,
            },{
            # BT_DWUG
            "r_0c":179,           # ohms/km
            "a_c":35.89e-3,

            "l_0":0.695e-3,       # H/km
            "l_inf":585e-6,       # H/km

            "b":1.2,
            "f_m":1000,           # kHz  

            "c_inf":55e-9,        # F/km
            "c_0":1e-9,
            "c_e":0.1,

            "g_0":0.5e-9,         # S/km
            "g_e":1.033
            },{
            # awg 24
            "r_0c":174.55888,     # ohms/km
            "a_c":0.053073,

            "l_0":617.29e-6,      # H/km
            "l_inf":478.97e-6,    # H/km

            "b":1.1529,
            "f_m":553.760,        # kHz  

            "c_inf":50e-9,        # F/km
            "c_0":0.0,
            "c_e":0.0,

            "g_0":234.87476e-15,  # S/km
            "g_e":1.38,
    }]

t=2 #Type 3 in transferfn.c


CHANNEL_BANDWIDTH = 4312.5 #from include/multiuser_load.h

"""
Fancy memoise decorator to make my life easy
"""
def memoize(fctn):
        memory = {}
        @functools.wraps(fctn)
        def memo(*args,**kwargs):
                haxh = cPickle.dumps((args, sorted(kwargs.iteritems())))

                if haxh not in memory:
                        memory[haxh] = fctn(*args,**kwargs)

                return memory[haxh]
        if memo.__doc__:
            memo.__doc__ = "\n".join([memo.__doc__,"This function is memoized."])
        return memo   

"""
Let the transfer function default to type 2;
    Allows for easy default change later
    Allows for easy 'case-based' changes
"""
def do_transfer_function(length,freq,type=3, measure="m"):
    """
    Z=impedance/l, Y=admittance/l
    Z=R+jwL, Y=G+jwC
    Z0=Characteristic Impedence, gamma=propagation constant
    Z0=sqrt(Z/Y), gamma=sqrt(Z*Y)
    Should Use PyGSL, but start off with cmath
    """
    
    #Length must be in KM
    if measure == "m":
        length /= 1000
    
    w = 2 * math.pi * freq

    #log.debug("Freq=%f,Len=%d",freq,length)
    Z = complex(_R(freq),w*_L(freq))
    Y = complex(_G(freq),w*_C(freq))

    Z0 = cmath.sqrt(Z/Y)
    gamma = cmath.sqrt(Z*Y)

    gammad=gamma*length

    Zs = complex(100,0)
    Zl = complex(100,0)
    #log.debug("w:%.3f",w)
    #log.debug("R:%e L:%e G:%e C:%e",_R(freq),_L(freq),_G(freq),_C(freq))
    #log.debug("Z:%s Y:%s Z0:%s Gamma:%s GammaD:%s",complex2str(Z),complex2str(Y),complex2str(Z0),complex2str(gamma),complex2str(gammad))
    upper = Z0 * ( 1/cmath.cosh(gammad)) #sech=1/cosh
    lower = Zs * ( (Z0/Zl) + cmath.tanh(gammad) ) + Z0 * ( 1+ (Z0/Zl)*cmath.tanh(gammad) )

    H = upper/lower
    H*=H
    #log.debug("Returned Transfer Function: %g",cmath.polar(H)[0])
    return cmath.polar(H)[0]        #polar returns (abs(h),real(h))

def _R(freq):
    """
    Return R Parameter for transfer function
    """
    c_partial = math.pow(material[t]["a_c"]*freq*freq+math.pow(material[t]["r_0c"],4),(0.25))

    try:
        if material[t]["r_0s"] > 0:
            s_partial = math.pow(material[t]["a_s"]*math.pow(freq,2)+math.pow(material[t]["r_0s"],4),(0.25))
            return (c_partial*s_partial)/(c_partial+s_partial)
    except KeyError:
        return c_partial
    
def _L(freq):
    """
    Return L Parameter for transfer function
    """
    upper=(material[t]["l_0"]+material[t]["l_inf"]*math.pow(freq*1e-3/material[t]["f_m"],material[t]["b"]))
    lower=(1+math.pow(freq*1e-3/material[t]["f_m"],material[t]["b"]))
    return (upper/lower)

def _C(freq):
    return material[t]["c_inf"]+material[t]["c_0"]*math.pow(freq,-material[t]["c_e"])

def _G(freq):
    return material[t]["g_0"]*math.pow(freq,material[t]["g_e"])

"""
Mathematical Utilities
"""
def dbmhz_to_watts(psd):
    return UndB(psd)*1e-3*CHANNEL_BANDWIDTH

def watts_to_dbmhz(psd):
    return TodB((psd*1e3)/CHANNEL_BANDWIDTH)

def UndB(input):
    try:
        return math.pow(10,input/10)
    except ValueError:
        log.error("Caught Exception on UndB(%f)"%input)
        raise ValueError

def TodB(input):
    try:
        return 10*math.log10(input)
    except ValueError:
        log.error("Caught Exception on TodB(%f)"%input)
        return 0

def freq_on_tone(K): #TODO
        """
        Assume ADSL downstream for now
        """
        return K * 4312.5 + 140156.25;

def complex2str(complex):
    return '{0:.3f}{1:+.3f}i'.format(complex.real,complex.imag)

def _Q_1(value):
    return math.sqrt(2)*sps.erfc(1-2*value)

def _Ne(M): #from http://www.docstoc.com/docs/21599433/Basics-of-Digital-Modulation
    rtM=math.sqrt(M)
    return ((2*rtM)-1)/rtM

#Uncoded SNR Gap (Probability-bit-error:Pe,N Nearest Neighbours:Ne)
def get_GAMMA(Pe,M):
    return (pow(_Q_1(Pe/_Ne(M)),2)/3)

#Combination Generator with replacements
def combinations(iterable, r,type=int):
    # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
    pool = tuple(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield tuple(pool[i] for i in indices)
        
#since numpy matrixes are a pain to convert to arrays
def mat2arr(matrix):
    return numpy.squeeze(numpy.asarray(matrix))

#makes things human
def psd2str(psd):
    assert(isinstance(psd,numpy.ndarray))
    return str(map(watts_to_dbmhz,psd))

if __name__ == "__main__":

    test_len=3000
    test_f=freq_on_tone(0)
    test_result=0.000853
    
    mine=do_transfer_function(test_len,test_f)
    assert mine == test_result, "Transfer function isn't right"
