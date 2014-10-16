def calc_fext_xtalk_gain(self,victim,xtalker,freq,dir):
    '''
    Calculate Far End XT Gain
    :adapted from channel_matrix.c 
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
            victim nt > xtalker nt (3)
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
    
    '''
    (vid,xid)=(victim.id,xtalker.id)
    #Check first if there is any shared sector (swapping screws with abs() operation)
    #Check if either v.lt or v.nt is between x.lt/x.nt
    if  (xtalker.lt<victim.lt and victim.lt<xtalker.nt) or (xtalker.lt<victim.nt and victim.nt<xtalker.nt):
        log.error("(V/X):(%d/%d):No Shared Sector"%(vid,xid))
        return 0

    
    if dir == "DOWNSTREAM":                 #If DOWNSTREAM, swap  and negate nt/lt's
        (D1,D2)=(-xtalker.lt,-xtalker.nt)
        (A,B)=(-victim.lt,-victim.nt)
    else: #Upstream
        (D1,D2)=(xtalker.nt,xtalker.lt)
        (A,B)=(victim.nt,victim.lt)

    
    #Shared length is the H2 Span
    shared_length=abs(max(A,D1)-min(B,D2))/1000.0
    #Head Length is the H1/H4 Span
    head_length=(A-D1)/1000.0
    #Tail Length is the H3 Span
    tail_length=(B-D2)/1000.0
    
    h1 = self.insertion_loss(head_length, freq)
    h2 = self.insertion_loss(shared_length, freq)
    h3 = self.insertion_loss(tail_length, freq)
    '''
    This _H takes away the biggest pain of the fext gain cases; 
    H1 > 0 in cases 1,3,4,6,7
    H2 active in all cases, but using max(v.nt,x.lt)
    H3 > 0 in cases 1,2,3
    H4 Not relevent for FEXT
    
    NB, insertion_loss(<0,f)=1
    
    '''
    H = h1*h2*h3
    
    try:
        gain = H * self.fext(freq, shared_length)
    except ValueError:
        log.error("Failed on (V/X):(%d/%d):(A,B,D1,D2)=(%d,%d,%d,%d):Shared:%f"%(vid,xid,A,B,D1,D2,shared_length))
        raise
    
    return gain  
