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
        self.p=numpy.zeros((self.bundle.N,self.MAXBITSPERTONE+1,self.MAXBITSPERTONE+1))
        self.b=numpy.asmatrix(numpy.zeros((self.bundle.N,self.bundle.K)))
        
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
            w=(w_max+w_min)/2.0
            
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
        while sum(linea.p) < linea.p_max :
            l1_max = 2 * l1_max  #This could be replaced by a bitshift?
            self.optimise_l2(w,l1_max,linea,lineb)

        utility.log.info("optimised_l1:max %d"%l1_max)
            
        while not self._converged(linea,power): #TODO
            utility.log.debug("Not Converged:l1(%e,%e)"%(l1_min,l1_max))
            l1 = (l1_max + l1_min)/2.0
            self.optimise_l2(w,l1,linea,lineb)
            power=sum(linea.p)
            if power > linea.p_max:
                l1_min=l1
            else:
                l1_max=l1
                
        utility.log.debug("Converged:l1(%e,%e)"%(l1_min,l1_max))

            
                
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
        while sum(lineb.p) < lineb.p_max:
            l2_max = 2 * l2_max  #This could be replaced by a bitshift?
            self.optimise_p(w,l1,l2_max,linea,lineb)

        utility.log.info("optimised_l2:max %d"%l2_max)
        
        while not self._converged(lineb,power): #TODO
            utility.log.debug("Not Converged:l2(%e,%e)"%(l2_min,l2_max))
            l2 = (l2_max + l2_min)/2.0
            self.optimise_p(w,l1,l2,linea,lineb)
            power=sum(lineb.p)
            if power > lineb.p_max:
                l2_min=l2
            else:
                l2_max=l2
                
        utility.log.debug("Converged:l2(%e,%e)"%(l2_min,l2_max))
                
    """
    Optimise Power (aka optimise_s)
    :from OSB_original.pdf paper
    """   
    def optimise_p(self,w,l1,l2,linea,lineb):
        utility.log.info("optimise_p(w:%f,l1:%f,l2:%f)"%(w,l1,l2))
        #for each subchannel
        for k in range(self.bundle.K):
            utility.log.info("optimise_p(w:%f,l1:%f,l2:%f,k:%d)"%(w,l1,l2,k))
            b1_max = 0 #overall max for b1
            b2_max = 0 #overall max for b2
            """
            #and each bit possibility
            for b1 in xrange(self.MAXBITSPERTONE):
                #for both lines
                for b2 in xrange(self.MAXBITSPERTONE):
                    L_k=self._L(w, l1, l2, linea, lineb, k, b1, b2)
                    if L_k > L_max:
                        L_max = L_k
                        #Store these new maxes
                        b1_max = b1
                        b2_max = b2
                        
                    #Powers need to be updated at some point. see PSD_vector.h
            #end of search double loop on b's
            """
            #Generate full matrix of bitrate potentials
            L=self._L_mat(w, l1, l2, linea, lineb, k)
            #find the argmax
            max=numpy.matrix.argmax(L)
            #turn it into a useable tuple index and assign to b_max's
            (b1_max,b2_max)=numpy.unravel_index(max, L.shape)
            
            utility.log.debug("Shape:%s"%str(L.shape))
            
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
        #return (abs(sum(line.p)-lastpower)<e)
        return True
    
    """
    Calculate the Lagrangian
    :from osb_original.pdf and osb.c
    """
    def _L(self,w,l1,l2,linea,lineb,k,b1,b2):
        result=0
        #Weighted Sections
        result+= (w*b1+(1-w)*b2)
        #Lagrangian Sections
        result-= (l1*linea.p[k])*(self.bundle.calc_psd(k)[0])#this isn't supposed to be xtalk, its p_matrix. No idea where the fuck it comes from tho.
        result-= (l2*lineb.p[k])*(self.bundle.calc_psd(k)[1])
        return result
    
    """
    Calculate the full lagrangian matrix for all bitrates
        L=W-L0-L1
        W_ij=wi+(1-w)j
        L0_ij=l1*linea.p[k]*P[0]
        L1_ij=l2*lineb.p[k]*P[1]
    """
    def _L_mat(self,w,l1,l2,linea,lineb,k):
        bitmax=self.MAXBITSPERTONE
        #Initialise Destination matrix
        W=numpy.asmatrix([[ w*i+(1-w)*j for i in xrange(bitmax)] for j in xrange(bitmax)])
        P=self.bundle.calc_psd(k)
        L0=l1*linea.p[k]*P[0]
        L1=l2*lineb.p[k]*P[1]
        L=W-L0-L1
        
        return L
        
        
        
        
