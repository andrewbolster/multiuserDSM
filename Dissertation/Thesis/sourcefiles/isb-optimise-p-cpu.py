def optimise_p_k(self,lambdas,weights):
    '''
    Optimise Power per tone, CPU bound version
    '''
    for k in range(self.bundle.K):
        #Convergence value check
        b_this=np.tile(0,self.bundle.N)
        b_last=np.tile(-1,self.bundle.N)
        log.debug("Launched channel %d search"%k)
        
        #Until convergence of this channels bitload
        while not (b_last==b_this).all():
            b_last=b_this.copy()       
            for line in xrange(self.bundle.N):
                lk_max=-self.defaults['maxval']
                b_max=[]
                #for each bit modification
                b_this[line]=0
                while b_this[line] <= self.MAXBITSPERTONE:
                    #The lagrangian value for this bit combination
                    lk=self._l_k(b_this,lambdas,weights,k)
                    if lk >= lk_max:
                        lk_max=lk
                        b_max=b_this[line]
                    b_this[line]+=1
                        
                #By now we have b_max for this user on this channel
                b_this[line]=b_max
            #at this point we are hopefully maximised
        self.b[k]=b_this
        self.p[k]=self.bundle.calc_psd(b_this,k)
    #end while
    
def optimise_p_k_alt(self,lambdas, weights):
    '''
    Optimise Power: Alternative version used to verify GPU algorithmic loop folding
    '''
    #Convergence value check
    b_this=np.tile(0,(self.bundle.K,self.bundle.N))
    b_last=np.tile(-1,(self.bundle.K,self.bundle.N))
    #Until convergence of these channels bitload
    while (b_last!=b_this).any():
        b_last=b_this.copy()       
        for k in range(self.bundle.K):
            log.debug("Launched channel %d search"%k)
            for line in xrange(self.bundle.N):
                lk_max=-self.defaults['maxval']
                b_max=[]
                #for each bit modification
                b_this[k,line]=0
                while b_this[k,line] <= self.MAXBITSPERTONE:
                    #The lagrangian value for this bit combination
                    lk=self._l_k(b_this[k],lambdas,weights,k)
                    if lk >= lk_max:
                        lk_max=lk
                        b_max=b_this[k,line]
                    b_this[k,line]+=1
                        
                #By now we have b_max for this user on this channel
                b_this[k,line]=b_max
            #at this point we are hopefully maximised
    self.b=b_this
    for k in range(self.bundle.K):
        self.p[k]=self.bundle.calc_psd(b_this[k],k)
    #end while
    
