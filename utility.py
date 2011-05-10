# This is all for experiments at the minute and not being called by anything... yet

__author__="bolster"
__date__ ="$02-Dec-2010 18:48:38$"

import cmath, numpy as np, sys#, scipy.special as sps
from math import *
import functools,cPickle
import logging

# Log everything, and send it to stderr.
log = logging.getLogger('multiuserdsm')
log.setLevel(logging.INFO)

h = logging.StreamHandler()
f = logging.Formatter('%(levelname)-7s %(module)s %(lineno)d %(message)s')
h.setFormatter(f)
log.addHandler(h)
h = logging.FileHandler('multiuserdsm.log')
h.setFormatter(f)
log.addHandler(h)

#Directories
rawdir="raw_results/"
profdir="profiles/"
graphdir="graphs/"

mp=True

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

tiny= 1e-300
half=  5.00000000000000000000e-01
one =  1.00000000000000000000e+00
two =  2.00000000000000000000e+00
erx =  8.45062911510467529297e-01

## Coefficients for approximation to  erf  in [0.84375,1.25] 
pa0  = -2.36211856075265944077e-03
pa1  =  4.14856118683748331666e-01
pa2  = -3.72207876035701323847e-01
pa3  =  3.18346619901161753674e-01
pa4  = -1.10894694282396677476e-01
pa5  =  3.54783043256182359371e-02
pa6  = -2.16637559486879084300e-03
qa1  =  1.06420880400844228286e-01
qa2  =  5.40397917702171048937e-01
qa3  =  7.18286544141962662868e-02
qa4  =  1.26171219808761642112e-01
qa5  =  1.36370839120290507362e-02
qa6  =  1.19844998467991074170e-02

def erf1(x):
    '''erf(x) for x in [0,0.84375]'''
    e, i = frexp(x)
    if abs(i)>28:
        if abs(i)>57:
            return 0.125*(8.0*x+efx8*x)
        return x + efx*x
    z = x*x
    r = pp0+z*(pp1+z*(pp2+z*(pp3+z*pp4)))
    s = one+z*(qq1+z*(qq2+z*(qq3+z*(qq4+z*qq5))))
    y = r/s
    return x + x*y

def erfc1(x):
    '''erfc(x)for x in [0,0.84375]'''
    e,i = frexp(x)
    if abs(i)>56:
        return one-x
    z = x*x
    r = pp0+z*(pp1+z*(pp2+z*(pp3+z*pp4)))
    s = one+z*(qq1+z*(qq2+z*(qq3+z*(qq4+z*qq5))))
    y = r/s
    if (x<0.25):
        return one-(x+x*y)
    else:
        r = x*y
        r += (x-half)
        return half - r

## Coefficients for approximation to  erf  in [0.84375,1.25] 

pa0  = -2.36211856075265944077e-03
pa1  =  4.14856118683748331666e-01
pa2  = -3.72207876035701323847e-01
pa3  =  3.18346619901161753674e-01
pa4  = -1.10894694282396677476e-01
pa5  =  3.54783043256182359371e-02
pa6  = -2.16637559486879084300e-03
qa1  =  1.06420880400844228286e-01
qa2  =  5.40397917702171048937e-01
qa3  =  7.18286544141962662868e-02
qa4  =  1.26171219808761642112e-01
qa5  =  1.36370839120290507362e-02
qa6  =  1.19844998467991074170e-02

def erf2(x):
    '''erf(x) for x in [0.84375,1.25]'''
    s = fabs(x)-one
    P = pa0+s*(pa1+s*(pa2+s*(pa3+s*(pa4+s*(pa5+s*pa6)))))
    Q = one+s*(qa1+s*(qa2+s*(qa3+s*(qa4+s*(qa5+s*qa6)))))
    if x>=0:
        return erx + P/Q
    return -erx - P/Q

def erfc2(x):
    '''erfc(x) for x in [0.84375, 1.25]'''
    return one-erf2(x)

## Coefficients for approximation to  erfc in [1.25,1/0.35]

ra0  = -9.86494403484714822705e-03
ra1  = -6.93858572707181764372e-01
ra2  = -1.05586262253232909814e+01
ra3  = -6.23753324503260060396e+01
ra4  = -1.62396669462573470355e+02
ra5  = -1.84605092906711035994e+02
ra6  = -8.12874355063065934246e+01
ra7  = -9.81432934416914548592e+00
sa1  =  1.96512716674392571292e+01
sa2  =  1.37657754143519042600e+02
sa3  =  4.34565877475229228821e+02
sa4  =  6.45387271733267880336e+02
sa5  =  4.29008140027567833386e+02
sa6  =  1.08635005541779435134e+02
sa7  =  6.57024977031928170135e+00
sa8  = -6.04244152148580987438e-02

def erf3(x):
    '''erf(x) for x in [1.25,2.857142]'''
    x0=x
    x = fabs(x)
    s = one/(x*x)
    R=ra0+s*(ra1+s*(ra2+s*(ra3+s*(ra4+s*(ra5+s*(ra6+s*ra7))))))
    S=one+s*(sa1+s*(sa2+s*(sa3+s*(sa4+s*(sa5+s*(sa6+s*(sa7+s*sa8)))))))
    z = ldexp(x0,0)
    r = exp(-z*z-0.5625)*exp((z-x)*(z+x)+R/S)
    if(x0>=0):
        return one-r/x
    else:
        return  r/x-one;

def erfc3(x):
    '''erfc(x) for x in [1.25,1/0.35]'''
    return one-erf3(x)

## Coefficients for approximation to  erfc in [1/.35,28]

rb0  = -9.86494292470009928597e-03
rb1  = -7.99283237680523006574e-01
rb2  = -1.77579549177547519889e+01
rb3  = -1.60636384855821916062e+02
rb4  = -6.37566443368389627722e+02
rb5  = -1.02509513161107724954e+03
rb6  = -4.83519191608651397019e+02
sb1  =  3.03380607434824582924e+01
sb2  =  3.25792512996573918826e+02
sb3  =  1.53672958608443695994e+03
sb4  =  3.19985821950859553908e+03
sb5  =  2.55305040643316442583e+03
sb6  =  4.74528541206955367215e+02
sb7  = -2.24409524465858183362e+01

def erf4(x):
    '''erf(x) for x in [1/.35,6]'''
    x0=x
    x = fabs(x)
    s = one/(x*x)
    R=rb0+s*(rb1+s*(rb2+s*(rb3+s*(rb4+s*(rb5+s*rb6)))))
    S=one+s*(sb1+s*(sb2+s*(sb3+s*(sb4+s*(sb5+s*(sb6+s*sb7))))))
    z  = ldexp(x0,0)
    r  =  exp(-z*z-0.5625)*exp((z-x)*(z+x)+R/S)
    if(z>=0):
        return one-r/x
    else:
        return  r/x-one;

def erfc4(x):
    '''erfc(x) for x in [2.857142,6]'''
    return one-erf4(x)

def erf5(x):
    '''erf(x) for |x| in [6,inf)'''
    if x>0:
        return one-tiny
    return tiny-one

def erfc5(x):
    '''erfc(x) for |x| in [6,inf)'''
    if (x>0):
        return tiny*tiny
    return two-tiny

#############
##inf = float('inf')
##nan = float('nan')
###########
inf = float(9e999)

def erf(x):
    '''return the error function of x'''
    f = float(x)
    if (f == inf):
        return 1.0
    elif (f == -inf):
        return -1.0
##    elif (f is nan):
##        return nan
    else:
        if (abs(x)<0.84375):
            return erf1(x)
        elif (0.84375<=abs(x)<1.25):
            return erf2(x)
        elif (1.25<=abs(x)<2.857142):
            return erf3(x)
        elif (2.857142<=abs(x)<6):
            return erf4(x)
        elif (abs(x)>=6):
            return erf5(x)
    
def erfc(x):
    '''return the complementary of error function of x'''
    f = float(x)
    if (f == inf):
        return 0.0
    elif (f is -inf):
        return 2.0
##    elif (f == nan):
##        return nan
    else:
        if (abs(x)<0.84375):
            return erfc1(x)
        elif (0.84375<=abs(x)<1.25):
            return erfc2(x)
        elif (1.25<=abs(x)<2.857142):
            return erfc3(x)
        elif (2.857142<=abs(x)<6):
            return erfc4(x)
        elif (abs(x)>=6):
            return erfc5(x)
CHANNEL_BANDWIDTH = 4312.5 #from include/multiuser_load.h

'''
Fancy memoise decorator to make my life easy
'''
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

'''
Let the transfer function default to type 2;
    Allows for easy default change later
    Allows for easy 'case-based' changes
'''
def do_transfer_function(length,freq,type=3, measure="m"):
    '''
    Z=impedance/l, Y=admittance/l
    Z=R+jwL, Y=G+jwC
    Z0=Characteristic Impedence, gamma=propagation constant
    Z0=sqrt(Z/Y), gamma=sqrt(Z*Y)
    Should Use PyGSL, but start off with cmath
    '''
    
    #Length must be in KM
    if measure == "m":
        length /= 1000
    elif measure == "km":
        pass
    else: raise TypeError('Improper measurement scale used') 
    
    w = 2 * pi * freq

    #log.debug("Freq=%f,Len=%d",freq,length)
    Z = complex(_R(freq),w*_L(freq))
    Y = complex(_G(freq),w*_C(freq))

    Z0 = cmath.sqrt(Z/Y)
    gamma = cmath.sqrt(Z*Y)

    gammad=gamma*length

    Zs = complex(100,0)
    Zl = complex(100,0)
    upper = Z0 * ( 1/cmath.cosh(gammad)) #sech=1/cosh
    lower = Zs * ( (Z0/Zl) + cmath.tanh(gammad) ) + Z0 * ( 1+ (Z0/Zl)*cmath.tanh(gammad) )

    H = upper/lower
    H*=H
    return cmath.polar(H)[0]        #polar returns (abs(h),real(h))

def _R(freq):
    '''
    Return R Parameter for transfer function
    '''
    c_partial = pow(material[t]["a_c"]*freq*freq+pow(material[t]["r_0c"],4),(0.25))

    try:
        if material[t]["r_0s"] > 0:
            s_partial = pow(material[t]["a_s"]*pow(freq,2)+pow(material[t]["r_0s"],4),(0.25))
            return (c_partial*s_partial)/(c_partial+s_partial)
    except KeyError:
        return c_partial
    
def _L(freq):
    '''
    Return L Parameter for transfer function
    '''
    upper=(material[t]["l_0"]+material[t]["l_inf"]*pow(freq*1e-3/material[t]["f_m"],material[t]["b"]))
    lower=(1+pow(freq*1e-3/material[t]["f_m"],material[t]["b"]))
    return (upper/lower)

def _C(freq):
    return material[t]["c_inf"]+material[t]["c_0"]*pow(freq,-material[t]["c_e"])

def _G(freq):
    return material[t]["g_0"]*pow(freq,material[t]["g_e"])

'''
Mathematical Utilities
'''
def dbmhz_to_watts(psd):
    return UndB(psd)*1e-3*CHANNEL_BANDWIDTH

def watts_to_dbmhz(e):
    return TodB((e*1e3)/CHANNEL_BANDWIDTH)

def UndB(input):
    try:
        return pow(10,input/10)
    except ValueError:
        #log.debug("Caught Exception on UndB(%f)"%input)
        raise ValueError
    except OverflowError:
        log.info("Overflowed on UndB(%f)"%input)
        raise OverflowError

def TodB(input):
    try:
        return 10*log10(input)
    except ValueError:
        #log.debug("Caught Exception on TodB(%f)"%input)
        return -np.inf #FIXME Either make it dynamic across the package or use something like -inf; 

def freq_on_tone(K): #TODO Memoize
        '''
        Assume ADSL downstream for now
        '''
        #return K * 4312.5 + 140156.25;
        assert False==True, "Someone Tried to use freq_on_tone"

def complex2str(complex):
    return '{0:.3f}{1:+.3f}i'.format(complex.real,complex.imag)

def _Q_1(value):
    return sqrt(2)*erfc(1-2*value)

def _Ne(M): #from http://www.docstoc.com/docs/21599433/Basics-of-Digital-Modulation
    rtM=sqrt(M)
    return ((2*rtM)-1)/rtM

#Uncoded SNR Gap (Probability-bit-error:Pe,N Nearest Neighbours:Ne)
def get_GAMMA(Pe,M):
    return (pow(_Q_1(Pe/_Ne(M)),2)/3)

#Combination Generator with replacements
def combinations(iterable, r,type=int):
    #FIXME Creation of a iterable for every tone for every loop is pointless
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
    return np.squeeze(np.asarray(np.copy(matrix)))

#makes things human
def psd2str(psd):
    assert(isinstance(psd,np.ndarray))
    return str(map(watts_to_dbmhz,psd))

#Class overload for cost values NOT USED ATM
class CostValue(object):
    def __set__(self,obj,val):
        self.val= val if val > 0 else sys.maxint

#bitload from int(used for lk max)
def bitload_from_id(id,N,mbpt):
    bitload=np.zeros(N)
    for i in range(N):
        bitload[i]=id%mbpt;
        id/=mbpt;
    return bitload

if __name__ == "__main__":

    test_len=3000
    test_f= 140156.25 #first channel for ADSL Downstream
    test_result=0.000853
    
    mine=do_transfer_function(test_len,test_f)
    assert mine == test_result, "Transfer function isn't right"
