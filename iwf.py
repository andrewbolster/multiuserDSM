"""
Algorithm Modules
"""

class IWF(Algorithm):
    """
    Iterative Water Filling
    """
    
    def __init__(self, bundle):
        super(Algorithm,self).preamble
        
        assert bundle is Bundle, "You dun goofed the Bundle!"
        if bundle.type == "ADSL_DOWNSTREAM":
        	p_initial=0.111
        else:
        	p_initial=0.02818
        
        # @type bundle Bundle
        self.bundle=bundle
        self.p = tile(p_initial, bundle.K) 		#size bundle.K
        bundle.calculate_snr()
        
        #TODO start timer system here
        iterations = rate_iterations = 0
        while iterations < 30 and rate_iterations < 10:
            while bundle.check_all_margins(2):
                for i,line in enumerate(bundle.lines):
                    if line.rate_target == None:
                        self.am_load_ra(line)
                    else:
                        self.am_load_fm(line)
                    bundle.calculate_snr()
                iterations+=1
            
            #assume the best
            targets_met=True
            
            for i,line in enumerate(bundle.lines):
                if line.rate_target == None:
                    continue
                if line.rate < line.rate_target:
                    targets_met=False
            
            if targets_met:
                break
            
        else: #will only run after 30 failed iterations or 10 failed 
            print "Convergence Failure after "+iterations+" iterations, "+rate_iterations+" rate iterations"
            return -1
        
        super(Algorithm,self).postscript
        return
    """
    def am_load_ra(self,line):
        line.b = zeroes(bundle.K)
        tone_full = zeroes(bundle.K)
        gamma_hat = pow(10,(gamma+margin-c_g))
        
        for tone in bundle.K:
            if tone_full[tone] != 0:
                delta_p[tone] = (pow(2,(line.b[tone]+1))-1) *  
    """