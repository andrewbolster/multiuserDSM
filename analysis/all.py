"""
Algorithm Parent module 
"""

import sys
import utility
import math
import numpy


class Algorithm(object):
    def __init__(self,bundle):
        self.name="Default Algorithm Name; You should never see this!"
        self.bundle=bundle
        
        """
        Algorithm Constants; usually set within child classes
        """
        self.MAXRATE = None #Need to set this for load_fm and ISB
        self.MAXPOWER = None
        self.MAXBITSPERTONE = 15


        
    
    def preamble(self):
        #Stuff that needs to be done for any algorithm        
        utility.log.info("Lets get this party started!")
    
    def postscript(self):
        utility.log.info("All Done Here")
        
    #FIXME I'm fairly sure nothing below here is layed out properly yet; am_load_ra is used in SOME algos...
    def am_load_ra(self,line):  #am_load.c used in IWF, Multiuser_greedy, RR_IWF, RR_OSB
        line.b = numpy.zeros(self.bundle.K)
        tone_full = numpy.tile(False,self.bundle.K)
        p = numpy.zeros(self.bundle.K)
        p_total=0
        b_total=0
        delta_p = []
        
        self._calc_delta_p(line)
        
        while (0 <= min(delta_p)): #TODO I know this doesn't work, Don't understand the breaking state see am_load.c find_min()
            tone = min(delta_p)
            line.b[tone]+=1
            p_total += delta_p[tone]
            p[tone] += delta_p[tone]
            
            if (line.b[tone] >= self.MAXBITSPERTONE):
                tone_full[tone] = True
            
            if (p_total > self.MAXPOWER):
                line.b[tone]-=1
                p_total -=delta_p[tone]
                break
            elif (p[tone] > utility.dbmhz_to_watts(-20)): #if this tone is greater than the masking power
                line.b[tone]-=1
                p[tone]-=delta_p[tone]
                p_total -= delta_p[tone]
                tone_full[tone] = True
            
            self._calc_delta_p(line)
        else:
            utility.log.info("All Tones are full!") #Breaking statement where min(delta_p) < 0
        
        self.update_psds(line)
        
        
        line.service = numpy.zeros(self.bundle.K) #Still have no idea what this does, but reset it to zero anyway
        
        return line.rate()
    
    def am_load_fm(self,line,half=0):  #am_load.c used in IWF, Multiuser_greedy, RR_IWF, Single_IWF
        #half:(-1:bottom,0:both,1:top)
        """
        Initialise Everything
        """
        line.b = numpy.zeros(self.bundle.K)
        tone_full = numpy.tile(False,self.bundle.K)     #set all tones to full
        p_total=0
        b_total=0
               
        """
        Calculate Initial bit/tone costs
        """
        delta_p=self._calc_delta_p(line,tone_full)
        
        """
        Deal with Halves
        """ #Why is there a check on K==512?
        if ( half == 1 ): #Top Half
            for tone in range(self.bundle.K/2):
                tone_full[tone] = True
                utility.log.debug("Top Half only")
        elif (half == -1): #Bottom Half
            for tone in range(self.bundle.K/2,self.bundle.K):
                tone_full[tone] = True        
                utility.log.debug("Bottom Half only")
        
        while (0 < min(delta_p)):
            tone = numpy.argmin(delta_p) #return the tone with the lowest bit-adding cost
            line.b[tone]+=1
            b_total += 1            #Keep track of yer bits!
            p_total += delta_p[tone]
            
            """Rate/Power/Bit Checking"""
            if (b_total == self.MAXRATE):
                utility.log.debug("n:%d,k:%d - rate-satisfied"%(line.id,tone))
                break #I'm done          
            if (line.b[tone] >= self.MAXBITSPERTONE):
                tone_full[tone] = True
                utility.log.info("n:%d,k:%d - tone full"%(line.id,tone))
            if (p_total > self.MAXPOWER):
                utility.log.info("n:%d,k:%d - exceeded power budget, rolling back"%(line.id,tone))
                line.b[tone]-=1
                b_total -=1
                p_total -=delta_p[tone]
                break
            
            #Recalculate bit/tone costs
            delta_p=self._calc_delta_p(line,tone_full)
        else:
            utility.log.info("All Tones are full!") #Breaking statement where min(delta_p) <= 0
        
        #Update powers
        self.update_psds(line)
                        
        line.service = numpy.zeros(self.bundle.K) #Still have no idea what this does, but reset it to zero anyway
        
        if b_total != self.MAXRATE: #This should be right as b_total is incrementally added
            utility.log.error("Could not reach target data rate. Desired:",self.MAXRATE," Achieved:",b_total)                
        return b_total
    
    """
    Calculate cost of additional bit on all the tones in a line
    :from am_load.c (implicit)
    Optimised from original version (Thank you wolframalpha)
    """
    def _calc_delta_p(self,line,tone_full):
        delta_p=numpy.zeros(self.bundle.K)
        for tone in self.bundle.K:
            if not tone_full[tone]:
                delta_p[tone] = (pow(2,(line.b[tone]-1)) * 3 - 2 )* self.bundle.gamma_hat/line.cnr[tone]
            else:
                delta_p[tone] = 100
"""
Class that describes the DSL bundle object
"""

#Global imports
import sys
import cmath
import numpy
import pydot
import pyparsing
import scipy.special as ss
import pylab
from math import pow,log10,sqrt

#Local Imports
from utility import *
from line import Line

class Bundle(object):
    def __init__(self,network_file,K=512,rates_file=None, weights_file=None):
        self.K = int(K)                # the number of DMT channels
        self.N = 0                  #Filled in on loading file, here to remind me
        self.lines = []            # the DSL line objects
        self.xtalk_gain = []    # XT Matrix (NB Must be FULLY instantiated before operating)

        #Assuming that each bundle is going to be one medium
        """
        Material Parameters
        """
        self.material={
            "r_0c"    :0,
            "a_c"    :0,
            "r_0s"    :0,
            "a_s"    :0,
        
            "l_0"    :0,
            "l_inf"    :0,
        
            "b"        :0,
            "f_m"    :0,
        
            "c_inf"    :0,
            "c_0"    :0,
            "c_e"    :0,
        
            "g_0"    :0,
            "g_e"    :0
            }
        self.MARGIN = 3.0   #Performance Margin db
        self.C_G = 0        #Coding Gain db
        self.GAMMA= get_GAMMA(1e-7, 4)       #Uncoded SNR Gap (setup elsewere)
        
        self.gamma_hat = pow(10,(self.GAMMA+self.MARGIN-self.C_G)/10)


        """
        Try to open and parse the network configuration file
        """
        try:
            with open(network_file,"r") as nf:
                for n,line in enumerate(nf):
                    # nt and lt are network termination (CO or RT) and line termination (CPE)
                    nt,lt = line.split(",")
                    self.lines.append(Line(nt,lt,n,self))

            self.N = len(self.lines)
            log.info("Read %d Lines",self.N)
            

        except IOError:
            log.error("Cannot open the network file: %s",network_file)
            sys.exit(1)


        log.info("Calculating the channel matrix")
        #Initialise the gains
        self.xtalk_gain = numpy.zeros((self.N,self.N,self.K))
        #The real work begins
        self.calc_channel_matrix()
        
        log.info("Printing xtalk_gain")
        print self.xtalk_gain
        #self.graph_channel_matrix()
        log.info("Running self check:")
        self.check_xtalk_gains() #This is only ever used once; could be sent into calc_channel_matrix?
        


    """
    Calculates the bundle channel gain matrix, generating xtalk_gain[][][]
    :from channel_matrix.c
    """
    def calc_channel_matrix(self): #TODO
        for k in range(self.K):                         #For Each Tone
            for i,li in enumerate(self.lines):          #Between Every Victim
                for j,lj in enumerate(self.lines):      # and every xtalker
                    if i == j:                          #If you're talking to yourself, do lazy transfer_fn
                        li.gain[k] = self.xtalk_gain[i,j,k] = li.transfer_fn(freq_on_tone(k))
                    else:                               #Otherwise look at XT
                        self.xtalk_gain[i,j,k] = self.calc_fext_xtalk_gain(li,lj,freq_on_tone(k),"DOWNSTREAM") #This makes more sense in passing line objects instead of id's
                    #log.debug("Set channel %d,%d,%d to %g",i,j,k,self.xtalk_gain[i,j,k])
    
    """
    Check Normalised XT Gains
    :from channel_matrix.c
    """
    def check_xtalk_gains(self):
        yeses=0
        """ Original Way to do it
        for k in self.K:            #tone
            for x,lineX in enumerate(self.lines):    #xtalker
                for v,lineV in enumerate(self.lines):#victim
                    if x==v: continue
                    if self.xtalk_gain[x,v,k]/lineV.gain[k] > 0.5:
                        yeses+=1
        """
        
        """Lets turn this on its head"""
        for v,victim in enumerate(self.lines):
            #listcomprehension for x,v,k on xtalk_gains and gain[k]
            gainratio=[self.xtalk_gain[x,v,k]/victim.gain[k] for x in range(self.N) for k in range(self.K)]
            yeses+=len([1 for i in gainratio if i>0.5])
        
        log.info("Total:%d,%%Yes:%f%%"%(len(gainratio),yeses/(1.0*len(gainratio))))
    
    """
    Calculate Far End XT Gain
    :from channel_matrix.c 
    Uses:
        transfer_func
        insertion_loss
        fext
    (Upstream/Downstream) Cases:
        victim and xtalker == lt
            victim length < xtalker length (4/8)
            victim length > xtalker length (6/2)
            victim and xtalker == length (5)
        victim lt < xtalker lt
            victim nt < xtalker nt (1/9)
            victim nt < xtalker nt (3)
        victim lt > xtalker lt
            victim nt > xtalker nt (9/1)
            victim nt < xtalker nt (7)
        victim and xtalker == nt
            victim length < xtalker length (8/4)
            victim length > xtalker length (2/6)
        victim nt <= xtalker lt (0)
        victim lt >= xtalker nt (0)
    
    What do the cases mean?
    http://www.nordstrom.nu/tomas/publications/BalEtAl_2003_ETSITM6_033w07.pdf
    "Proposed Method on Crosstalk Calculations in a Distributed Environment"
    In Section 2: Top line is 'victim', A->B is nt->lt
    
    
    4 bundle segment types for 2 lines (xtalker and victim);
    H1:xtalker before victim (1,4,7)
    H2:xtalker and victim (all)
    H3:victim after xtalker(1,2,3)
    H4:victim before xtalker(3,6)
    
    """
    def calc_fext_xtalk_gain(self,victim,xtalker,freq,dir): #TODO Currently does UPSTREAM only, swap lt/nt for victim and xtalker for downstream
               
        if dir == "DOWNSTREAM":                 #If DOWNSTREAM, swap nt/lt's
            (dnt,dlt)=(xtalker.lt,xtalker.nt)
            (vnt,vlt)=(victim.lt,victim.nt)
        else:
            (dnt,dlt)=(xtalker.nt,xtalker.lt)
            (vnt,vlt)=(victim.nt,victim.lt)
            
        
        h1 = self.insertion_loss(abs(dnt-vnt), freq)
        h2 = self.insertion_loss(max(vnt,dnt)-dlt, freq)
        h3 = self.insertion_loss(dlt-vlt, freq)
        """
        This _H takes away the biggest pain of the fext gain cases; 
        H1 > 0 in cases 1,3,4,6,7 (also acting as H4 using abs())
        H2 active in all cases, but using max(v.nt,x.lt)
        H3 > 0 in cases 1,2,3
        """
        H = h1*h2*h3
        
        gain = H * self.fext(freq, max(vnt,dnt)-max(vlt,dlt))
        #log.debug("Returned xtalk_gain: %g",gain)

        return gain  
    
    """
    Calculate FEXT for given Freq and Length
    Uses:
        UndB
    Model from BT's simulation parameters document cp38-2 http://www.niccstandards.org.uk/files/current/NICC%20ND%201513%20(2010-01).pdf
    :from channel_matrix.c
    """
    def fext(self,freq,length):
       
        x=-55   #This is different from the NICC Model  #FUDGE
        x+= 20*log10(freq/90e3)         #f in Hz
        x+= 6*log10(len(self.lines)-1)/(49) #n disturbers 
        x+= 10*log10(length)            #shared length in km
        
        try:
            return UndB(x)
        except ZeroDivisionError:
            return 0
        
    """
    Calculate Insertion Loss
    :from transfer_fn.c
    """    
    def insertion_loss(self,length,freq): #TODO: Reintroduce factor removed from C version?
        """
        insertion loss makes more sense in the bundle as most of the time,
        it is looking at a common length between two lines
        """
        if length > 0:
            return do_transfer_function(length,freq )
        else: 
            return 1 #Leads straight into multiplication; so ends up x*0=0 otherwise
    
    """
    Calculate snr statistics for all lines in bundle
    :from snr.c
    """
    def calculate_snr(self):
        for line in self.lines:
            line.sanity()
            noise = numpy.zeros(self.K)
            
            for tone in xrange(self.K):
                
                noise = line.calc_fext_noise(tone) + line.alien_xtalk(tone) + dbmhz_to_watts(line.noise)
                line.cnr[tone] = line.gain[tone]/noise
                #log.debug("b%d,b%e,g%e"%(line.b[tone],(pow(2,line.b[tone])-1),(self.gamma_hat/line.cnr[tone])))
                line.p[tone] = watts_to_dbmhz((pow(2,line.b[tone])-1)*(self.gamma_hat/line.cnr[tone])) #TODO gamma_hat/cnr[tone] functionalise
                if (line.p[tone] < line.MINPSD) and (line.b[tone] != 0):
                    line.p[tone] = line.MINPSD
                line.snr[tone] = dbmhz_to_watts(line.p[tone])*line.cnr[tone]
                line.gamma_m[tone] = 10*log10(line.snr[tone]/pow(2,line.b[tone]-1))
                
            #line.symerr = [ self._calc_sym_err(line,xtalker) for xtalker in xrange(self.K)] #TODO _calc_sym_err
            line.p_total = sum(map(dbmhz_to_watts,line.p))
            line.b_total = sum(line.b)
            
            """
            With This whole Fractional thing; 
            Are b/_b ever different? If its a fractional line is b[k] ever needed?
            """
    """
    Calculate Symbol Error Rate on line
    :from symerr.c
    """
    def _calc_sym_err(self,line,k): #TODO
        
        M=pow(2,line.b[k])
               
        return 1 - (1 - (1 - 1/sqrt(M))*ss.erf(sqrt((3*line.snr[k])/(2*(M-1)))))
    
    """
    Checks inter-service margins
    :from snr.c
    #I suspect this can be ignored
    """
    def check_all_margins(self,gap):
        pass
            
    
    """
    Generate PSD vector matrix between two lines and return element
    :from psd_vector.c
    #TODO Memoise
    """
    def calc_psd(self,b1,b2,k):
        #Generate A Matrix (See Intro.pdf 2.23)
        A=numpy.asmatrix(numpy.zeros((self.N,self.N)))
        for i,linea in enumerate(self.lines):
            for j, lineb in enumerate(self.lines):
                if i==j: 
                    A[i,j]=1
                else:
                    A[i,j]=-self._f(linea,k)*self._h2(lineb,linea,k)/self._h2(lineb,linea,k)

        #Generate B Matrix (same reference)
        B=numpy.asmatrix(numpy.zeros(self.N))
        for i,line in enumerate(self.lines):
            B[0,i]=-self._f(line,k)*(dbmhz_to_watts(-140)+line.alien_xtalk(k))/self._h2(lineb,linea,k)
        
        #Yeah, lets twist again!
        B=B.T
                
        #Everyone loves linear algebra...dont they?
        if (False): #QR or regular way?
            q,r = numpy.linalg.qr(A)
            assert numpy.allclose(A,numpy.dot(q,r))
            p= numpy.dot(q.T,B)
            P=numpy.dot(numpy.linalg.inv(r),p)
        else:
            P=numpy.linalg.solve(A,B)
        
        #Because I'm paranoid
        assert numpy.allclose(B,(numpy.dot(A,P))), log.error("Turns out linear algebra is hard")
        #log.info("Calc_psd(b1:%d,b2:%d,k:%d)=P %s"%(b1,b2,k,P))
        return P

        
        
                    
    """
    Transfer function lookup - |h_{i,j}|^2
    This more or less does nothing but return xtalk_gain but is closer to the math
    """
    def _h2(self,line1,line2,channel):
        return self.xtalk_gain[line1.id][line2.id][channel]
    
    """
    F- Gamma function - 
    f(line,k)=\Gamma(2^{b_n(k)} -1)
    :from psd_vector.c
    #TODO I don't understand this
    """
    def _f(self,line,k):
        if hasattr(line,'g'):
            g=line.g[k] #TODO gamma values from symerr.c ?
        else:
            g=(12.8) #I'm lazy, using default value from intro.pdf
            
        b=line.b[k]
        return pow(10,(g+3)/10)*(pow(2,b)-1) #TODO initially, b=0, so this doesnt work
    """    
    Pretty Print channel Matrix
    #TODO I've got no idea how to display this....
    """
    def graph_channel_matrix(self):
        
        pylab.contourf(self.xtalk_gain[...,0])
        pylab.figure()
        pylab.show#!/usr/bin/env python

'''
generates call graph of given python code file
in dot format input for graphviz.

limitations:
* statically tried to figure out functions calls
* does not understand classes
* algorithm is naive and may not statically find
  all cases
'''

import sys
import parser
import symbol, token
import pprint
import optparse

try: s = set()
except: import sets; set = sets.Set

def annotate_ast_list(ast_list):
    code = ast_list[0]
    if code in symbol.sym_name: code = symbol.sym_name[code]
    else: code = token.tok_name[code]
    ast_list[0] = code

    for index, item in enumerate(ast_list):
        if index == 0: continue
        if isinstance(item, list):
            ast_list[index] = annotate_ast_list(item) 
    return ast_list
 
def get_atom_name(atom):
    first_child = atom[1]
    first_child_code = first_child[0]
    if first_child_code != token.NAME: return None
    return first_child[1]

def get_fn_call_data(ast_list):
    if len(ast_list) < 3: return None
    first_child, second_child = ast_list[1:3]
    first_child_code = first_child[0]
    if first_child_code != symbol.atom: return None
    fn_name = get_atom_name(first_child)

    second_child_code = second_child[0]
    if second_child_code != symbol.trailer: return None
    
    if len(second_child) < 3: return None
    if second_child[1][0] == token.LPAR and second_child[-1][0] == token.RPAR:
        return fn_name
    else: return None

def find_fn_call(ast_list, calls):
    code = ast_list[0]
    if code == symbol.power:
        fn_name = get_fn_call_data(ast_list)
        if fn_name != None and getattr(__builtins__, fn_name, None) == None: calls.add(fn_name) 
   
    for item in ast_list[1:]:
        if isinstance(item, list):
            find_fn_call(item, calls)

def process_fn(fn_ast_list, call_graph):
    dummy, dummy, func_name = fn_ast_list[:3]
    dummy, func_name = func_name

    calls = set()
    find_fn_call(fn_ast_list, calls)

    call_graph[func_name] = list(calls)

def construct_call_graph(ast_list, call_graph):
    code = ast_list[0]
    if code == symbol.funcdef:
        process_fn(ast_list, call_graph)

    for item in ast_list[1:]:
        if isinstance(item, list):
            construct_call_graph(item, call_graph)

    return call_graph

def generate_dot_code(python_code):
    ast = parser.suite(python_code)
    ast_list = parser.ast2list(ast)
    #annotated_ast_list = annotate_ast_list(ast_list)
    #pprint.pprint(annotated_ast_list)

    call_graph = {}
    construct_call_graph(ast_list, call_graph)
    #pprint.pprint(call_graph)

    dot = []

    dot.append("digraph G {")
    dot.append("rankdir=LR")
    for from_fn, to_fns in call_graph.iteritems():
        if not to_fns:
            dot.append('%s;' % from_fn)

        for to_fn in to_fns:
            if to_fn not in call_graph: continue
            dot.append('%s -> %s;' % (from_fn, to_fn))
    dot.append("}")

    return '\n'.join(dot)

if __name__ == '__main__':
    oparser = optparse.OptionParser()

    oparser.add_option('-i', '--input-file', default=None, metavar='FILE', help='python code file to process')

    options, args = oparser.parse_args()

    if options.input_file:
        python_code = open(options.input_file).read()
    else:
        python_code = sys.stdin.read()

    dot_code = generate_dot_code(python_code)
    print dot_code
"""
Algorithm Modules
"""

#Global Imports
import numpy
import math

#Local Imports
from bundle import Bundle
from algorithm import Algorithm

class ISB(Algorithm):
    """
    Iterative Spectrum Balancing
    """
    def __init__(self):
        self.name = "ISB"
        self.preamble
        
        self.current_rate = numpy.zeros(self.bundle.K)
        """
Algorithm Modules
"""

#Global Imports
import numpy
import math

#Local Imports
from bundle import Bundle
from algorithm import Algorithm
import utility

class IWF(Algorithm):
    """
    Iterative Water Filling
    """
    
    def __init__(self):
        self.name = "IWF"
        
        assert self.bundle is Bundle # "You dun goofed the Bundle!"
        
        """
        Set initial values for power
        """
        if self.bundle.type == "ADSL_DOWNSTREAM":
            p_initial = 0.111
        else:
            p_initial = 0.02818
        
        self.p = numpy.tile(p_initial, self.bundle.K) 		#size bundle.K
        self.bundle.calculate_snr()
        
        self.rate_targets = numpy.tile()
        
        #TODO start timer system here
        self.preamble()

        iterations = rate_iterations = 0
        """
        Iteratively Execute am_load_X and check rate targets
        """
        while iterations < 30 and rate_iterations < 10:
            while self.bundle.check_all_margins(2):
                for i,line in enumerate(self.bundle.lines):
                    if line.rate_target == None:    #TODO Where does line.rate_target[] get set?
                        self.am_load_ra(line)       #am_load_X functions in parent class
                    else:
                        self.am_load_fm(line)
                    self.bundle.calculate_snr()
                iterations+=1
            
            #assume the best
            targets_met=True
            
            """
            Target Checking
            """
            for i,line in enumerate(self.bundle.lines): #TODO List Comprehension/Reduce?
                if line.rate_target == None:
                    continue
                if line.rate < line.rate_target:
                    targets_met=False
            
            if targets_met: #we're done here
                utility.log.info("Rate Targets satisfied after %d rate iterations"%rate_iterations)
                break
            
            rate_iterations+=1 #start a new attempt
            
        else: #will only run after 30 failed iterations or 10 failed 
            utility.log.warning("Convergence Failure after "+iterations+" iterations, "+rate_iterations+" rate iterations")
        
        self.postscript()
"""
Class that describes a DSL Line object
"""
#Global Imports
import numpy

#Local Imports
import utility

class Line(object):
    def __init__(self,nt,lt,id,bundle):
        self.MINPSD=-60     #from am_load.c
        
        # nt and lt are network termination (CO or RT) and line termination (CPE)
        self.bundle = bundle
        self.nt = int(nt)
        self.lt = int(lt)
        self.length = self.nt - self.lt
        self.gain = numpy.zeros(bundle.K) 
        self.p = numpy.tile(self.MINPSD,bundle.K) #PSD of this line
        self.id = id #could this be removed as an array of lines?
        self.noise = -140 #Assuming standard background noise
        self.type = 2   #I'm assuming this declares the material of the line
                        #so in theory it can be removed from transfer_fn

        self.p_total = 0    #configured in bundle.calculate_snr
        self.b_total = 0    #ditto
        
        self.p_max = 1 #P=Total AFE power / delta_Freq #TODO
        
        self.snr = numpy.zeros(bundle.K) 
        self.cnr = numpy.zeros(bundle.K) 
        self.gamma_m = numpy.zeros(bundle.K) 
        self.b = numpy.tile(2,bundle.K) #bit loading on each subchannel
        
        
        
    def __str__(self):
        """
        Super cool graphics
        """
        s = ""
        s += "-" * (self.nt / 250)
        s += "|"
        s += "-" * ((self.lt-self.nt) / 250)
        s += "|"
        return s
    
    """
    Return the RLCG Parameterised Transfer Function
    """   
    def transfer_fn(self, freq):
        #C Version uses bundle as a linked list of lines so most of this function is redundant.
        return utility.do_transfer_function(abs(self.length), freq ,type=self.type)
    
    """
    Line Sanity Check
    :from snr.c check_line_sanity
    Called from bundle
    """
    def sanity(self):
        for k in range(self.bundle.K):
            assert self.b[k] >= 0, utility.log.error("Tone %d less than zero:%e"%(k,self.b[k]))
            assert self.gain[k] >= 0, utility.log.error("Gain %d less than zero:%e"%(k,self.gain[k]))
        utility.log.debug("Line %d is sane"%self.id)
        return
    
    """
    Calculate Far End XT noise on a channel in watts
    :from snr.c
    """
    def calc_fext_noise(self,channel):
        noise=0
        for xtalker in self.bundle.lines:
            if xtalker == self : continue   #stop hitting yourself!
            else:
                noise += utility.dbmhz_to_watts(self.p[channel])*self.bundle._h2(xtalker,self,channel)
        return noise
    
    """
    Return AlienXT noise
    :from snr.c
    """
    def alien_xtalk(self,channel): #TODO
        return 0
    
    """
    Return current rate on line
    """
    def rate(self):
        return sum(self.b)"""
The main module for the multiuserdsm framework
"""

from optparse import OptionParser
from bundle import Bundle
from utility import *
from osb import OSB

log.info("Starting Up...")

parser = OptionParser()

parser.add_option("-n","--network", dest="network", help="read network configuration from FILE",metavar="NETWORK",default="test.net")
parser.add_option("-K","--tones", dest="K", help="specify number of DMTCHANNELS",metavar="K",default=512)

(options,args) = parser.parse_args()

def run():
    bundle = Bundle(options.network,options.K)
    """
    Perform algorithm selection and option passing
    """
    algo = OSB(bundle)
    algo.run()"""
Algorithm Modules
"""

#Global Imports
import numpy
import math

#Local Imports
from bundle import Bundle
from algorithm import Algorithm
import utility

class OSB(Algorithm):
    """
    Optimum Spectrum Balancing
    """
    def run(self):
        self.name = "OSB"
        #Aimed-for rates per-line
        self.rate_targets = numpy.tile(100,self.bundle.N)
        #Tolerance of rate-targets (how close is close enough)
        self.target_tolerance=1
        
        self.preamble
        self.w=[]
        #How the fuck does this iterate over the bundle?
        for linea in self.bundle.lines:
            for lineb in self.bundle.lines:
                #Ignore comparing lines to themselves
                if linea.id != lineb.id:
                    self._optimise_w(linea,lineb)
                    utility.log.info("Optimised W for %d,%d"%(linea.id,lineb.id))
        self.postscript
        return  
    
    """
    Optimise W -Currently not used
    """
    def optimise_weights(self,linea,lineb):
        
        #FIXME There is no way this makes sense for multi-user system >2
        if (self.rate_targets[lineb.id]==0):
            self._optimise_w(linea,lineb)
        elif (self.rate_targets[lineb.id]==0):
            #Swap and carry on regardless
            self._optimise_w(lineb,linea)
        else:
            utility.log.error("What the fuck am I supposed to do now?!")
        
    
    """
    Optimise W- working function
    
    """
    def _optimise_w(self,linea,lineb):
        w_max=1
        w_min=0
        while (abs(linea.rate()-self.rate_targets[linea.id]) > self.target_tolerance):
            w=(w_max+w_min)/2
            self.optimise_l1(w,linea,lineb)
            if linea.rate() > self.rate_targets[linea.id]:
                w_max=w
            else:
                w_min=w
    
    """
    Optimise Lambda 1
    :from OSB_original.pdf paper
    """
    def optimise_l1(self,w,linea,lineb):
        l1_max=1
        l1_min=0
        power=-1.0 #initialised so that converged will work
        utility.log.info("optimise_l1(w:%f)"%(w))
        #First need to find the max permissable value of lambda_n (l1_max)
        while sum(linea.p) > linea.p_max:
            l1_max = 2 * l1_max  #This could be replaced by a bitshift?
            self.optimise_l2(w,l1_max,linea,lineb)
        while not self._converged(linea,power): #TODO
            l1 = (l1_max + l1_min)/2
            self.optimise_l2(w,l1,linea,lineb)
            power=sum(linea.p)
            if power > linea.p_max:
                l1_min=l1
            else:
                l1_max=l1
            
                
    """
    Optimise Lambda 2
    :from OSB_original.pdf paper
    """
    def optimise_l2(self,w,l1,linea,lineb):
        l2_max=1
        l2_min=0
        power=-1.0 #initialised so that converged will work
        utility.log.info("optimise_l2(w:%f,l1:%f)"%(w,l1))

        #First need to find the max permissable value of lambda_n (l1_max)
        while sum(lineb.p) > lineb.p_max:
            l2_max = 2 * l2_max  #This could be replaced by a bitshift?
            self.optimise_s(w,l1,l2_max,linea,lineb)
            
        while not self._converged(lineb,power): #TODO
            l2 = (l2_max + l2_min)/2
            self.optimise_p(w,l1,l2,linea,lineb)
            power=sum(lineb.p)
            if power > lineb.p_max:
                l2_min=l2
            else:
                l2_max=l2
                
    """
    Optimise Power (aka optimise_s)
    :from OSB_original.pdf paper
    """   
    def optimise_p(self,w,l1,l2,linea,lineb):
        L_max=-1.0 #max L_k seen
        b1_max = 0 #overall max for b1
        b2_max = 0 #overall max for b2
        #for each subchannel
        for k in range(self.bundle.K):
            #and each bit possibility
            for b1 in xrange(self.MAXBITSPERTONE):
                #for both lines
                for b2 in xrange(self.MAXBITSPERTONE):
                    L_k=self._L(w, l1, l2, linea, lineb, k, b1, b2)
                    if L_k > L_max:
                        L_max = L_k
                        #Store these new maxes
                        linea.b[k] = b1_max = b1
                        lineb.b[k] = b2_max = b2
                        
                    #Powers need to be updated at some point. see PSD_vector.h
            #end of search double loop on b's
            #Update 'best' bit loads
            linea.b[k]=b1_max
            lineb.b[k]=b2_max
        #Update powers
        self.bundle.calculate_snr()
        utility.log.info("optimised_p(w:%f,l1:%f,l2:%f)=(%d,%d)"%(w,l1,l2,b1_max,b2_max))
        
        #The reason why this seems to make sense HERE is that
        #in theory, the above double loop would be 'instant'.
            

    """
    Convergence check for both lambda optimisations
    """                 
    def _converged(self,line,lastpower):
        #How much convergence is satisfactory?
        e=0.01
        return (abs(sum(line.p)-lastpower)<e)
    
    """
    Calculate the Lagrangian
    :from osb_original.pdf and osb.c
    """
    def _L(self,w,l1,l2,linea,lineb,k,b1,b2):
        result=0
        #Weighted Sections
        result+= (w*b1+(1-w)*b2)
        #Lagrangian Sections
        result-= (l1*linea.p[k])*(self.bundle.calc_psd(b1,b2,k)[0])#this isn't supposed to be xtalk, its p_matrix. No idea where the fuck it comes from tho.
        result-= (l2*lineb.p[k])*(self.bundle.calc_psd(b1,b2,k)[1])
        return result
# This is all for experiments at the minute and not being called by anything... yet

__author__="bolster"
__date__ ="$02-Dec-2010 18:48:38$"

import math, cmath, numpy, scipy.special as sps

import logging

# Log everything, and send it to stderr.
log = logging.getLogger('multiuserdsm')
log.setLevel(logging.DEBUG)

h = logging.StreamHandler()
f = logging.Formatter('%(levelname)-7s %(module)s %(lineno)d %(message)s')
h.setFormatter(f)
log.addHandler(h)

material={#taken from type 2 in transfer_fn.c
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
    #log.debug("Returned Transfer Function: %g",cmath.polar(H)[0])
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
    upper=(material["l_0"]+material["l_inf"]*math.pow(freq*1e-3/material["f_m"],material["b"]))
    lower=(1+math.pow(freq*1e-3/material["f_m"],material["b"]))
    return (upper/lower)

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
    try:
        return math.pow(10,input)
    except ValueError:
        log.error("Caught Exception on UndB(%f)"%input)
        raise ValueError

def watts_to_dbmhz(psd):
    return TodB(psd*1e3)

def TodB(input):
    try:
        if input != 0.0: return 10*math.log10(input)
        else: return 0
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

if __name__ == "__main__":

    print do_transfer_function(2,8.5e4)
