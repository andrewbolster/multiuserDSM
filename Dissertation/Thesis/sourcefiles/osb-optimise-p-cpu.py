def optimise_p_k(self,lambdas,weights,K,Kmax):
    '''
    Optimise Power per tone, CPU bound version
    '''
    for k in range(K,Kmax):
        log.debug("Launched channel %d search"%k)
        lk_max=-self.defaults['maxval']
        b_max=[]
        #for each bit combination
        b_combinator=combinations(range(self.MAXBITSPERTONE), self.bundle.N)
        
        for b_combo in b_combinator:
            b_combo=np.asarray(b_combo)
            #The lagrangian value for this bit combination
            lk=self._l_k(b_combo,lambdas,weights,k)
            if lk >= lk_max:
                lk_max=lk
                b_max=b_combo
        #By now we have b_max[k]

        assert len(b_max)>0, "No Successful Lk's found,%s"%b_max
        self.p[k]=self.bundle.calc_psd(b_max,k)
        self.b[k]=b_max
        #print "CPU LKmax %d:%s:%s:%s"%(k,str(lk_max),str(b_max),str(self.p[k]))

    #end for
   
