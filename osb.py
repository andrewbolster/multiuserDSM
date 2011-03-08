"""
Algorithm Modules
"""

#Global Imports
import numpy
import math

#Local Imports
from bundle import Bundle
from algorithm import Algorithm

class OSB(Algorithm):
    """
    Optimum Spectrum Balancing
    """
    
    def __init__(self):
        self.name = "OSB"
        self.preamble
        self.w=[]
        #How the fuck does this iterate over the bundle?
        for linea in self.bundle.lines:
            for lineb in self.bundle.lines:
                #Ignore comparing lines to themselves
                if linea.id != lineb.id:
                    self.optimise_l1(w,linea,lineb)
        self.postscript
        return  
    
    #Optimise Lambda1 (from osb_original)
    def optimise_l1(self,w,linea,lineb):
        l1_max=1
        l1_min=0
        power=-1.0 #initialised so that converged will work
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
            
                
    #Optimise Lambda2 (from osb_original)
    def optimise_l2(self,w,l1,linea,lineb):
        l2_max=1
        l2_min=0
        power=-1.0 #initialised so that converged will work
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
        
    def optimise_p(self,w,l1,l2,linea,lineb):
        L_max=-1.0 #max L_k seen
        b1_max = 0 #overall max for b1
        b2_max = 0 #overall max for b2
        
        #for each subchannel
        for k in range(self.bundle.K):
            #and each bit possibility
            for b1 in range(self.MAXBITSPERTONE):
                #for both lines
                for b2 in range(self.MAXBITSPERTONE):
                    L_k=self._L(w, l1, l2, linea, lineb, k, b1, b2)
                    if L_k > L_max:
                        L_max = L_k
                        linea.b[k] = b1_max = b1
                        lineb.b[k] = b2_max = b2
                    #Powers need to be updated at some point. see PSD_vector.h
            #end of search double loop on b's
            linea.b[k]=b1_max
            lineb.b[k]=b2_max
            
                     
    #has a line power converged?
    def _converged(self,line,lastpower):
        #How much convergence is satisfactory?
        e=0.01
        return (abs(sum(line.p),lastpower)<e)
    
    #Lagrangian
    def _L(self,w,l1,l2,linea,lineb,k,b1,b2):
        result=0
        #Weighted Sections
        result+= (w*b1+(1-w)*b2)
        #Lagrangian Sections
        result-= (l1*linea.p[k])*(self.bundle.xtalk_gain[b1][b2][k])##It wasnt
        result-= (l2*lineb.p[k])*(self.bundle.xtalk_gain[b1][b2][k])
        return result
