"""
Algorithm Parent module 
"""

import sys

class Algorithm(object):
    def __init__(self):
        self.name="Default Algorithm Name; You should never see this!"
        pass
    
    def am_load_ra(self,line):
        line.b = zeroes(bundle.K)
        tone_full = tile(False,bundle.K)
        gamma_hat = pow(10,(gamma+margin-c_g)/10)
        p = zeroes(bundle.K)
        p_total=0
        b_total=0
        
        self._calc_delta_p()
        
        while (0 <= (tone = min(delta_p))): #TODO I know this doesn't work, Don't understand the breaking state see am_load.c find_min()
            line.b[tone]+=1
            p_total += delta_p[tone]
            p[tone] += delta_p[tone]
            
            if (line.b[tone] >= self.MAXBITSPERTONE):
                tone_full[tone] == True
            
            if (p_total > p_):
                line.b[tone]-=1
                p_tot -=delta_p[tone]
                break
            elif (p[tone] > dbmhz_to_watts(-20)) #if this tone is greater than the masking power
                line.b[tone]-=1
                p[tone]-=delta_p[tone]
                p_total -= delta_p[tone]
                tone_full[tone] = True
            
            self._calc_delta_p()
        else:
            print("All Tones are full!") #Breaking statement where min(delta_p) < 0
        
        for tone in xrange(bundle.K):
            line.psd[tone] = watts_to_dbmhz((math.pow(2,line.b[tone])-1)*(gamma_hat/line.cnr[tone]))
            if (line.psd[tone] < min_psd) and (line.b[tone] != 0):
                line.psd[tone] = min_psd
        
        b_total = sum(line.b)
        
        line.service = zeroes(bundle.K) #Still have no idea what this does, but reset it to zero anyway
        
        return b_total
    
    def am_load_fm(self,line):
        """
        Need to discuss with AM merging am_load_* and passing in a 'satisfaction' function;
            am_load_fm == am_load(sumofline.b[]=b_target
            am_load_ra == am_load(sumofline.p[]=p_target
        """
        return self.am_load_ra(line)
        
    def _calc_delta_p(self): #TODO:Used in several am_load related functions; need to make variable safe
        for tone in bundle.K:
            if not tone_full[tone]:
                self.delta_p[tone] = (pow(2,(line.b[tone]-1)) * 3 - 2 )* gamma_hat/line.cnr[tone]
            else:
                self.delta_p[tone] = 100