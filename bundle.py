"""
Class that describes the DSL bundle object
"""

import sys
import cmath
from line import Line

class Bundle(object):
	def __init__(self,network_file,K):
		self.K = K				# the number of DMT channels
		self.lines = []			# the DSL line objects
		self.xtalk_gain = None	# XT Matrix (NB Must be FULLY instantiated before operating)

        #Assuming that each bundle is going to be one medium
        """
        Material Parameters
        """
        self.material={
            "r_0c"	:0,
            "a_c"	:0,
            "r_0s"	:0,
            "a_s"	:0,

            "l_0"	:0,
            "l_inf"	:0,

            "b"		:0,
            "f_m"	:0,

            "c_inf"	:0,
            "c_0"	:0,
            "c_e"	:0,

            "g_0"	:0,
            "g_e"	:0
            }

		"""
		Try to open and parse the network configuration file
		"""
		try:
			with open(network_file,"r") as nf:
				for line in nf:
					# nt and lt are network termination (CO or RT) and line termination (CPE)
					nt,lt = line.split(",")
					self.lines.append(Line(nt,lt))

			self.N = len(self.lines)

		except IOError:
			print "Cannot open the network file",network_file
			sys.exit(1)


		"""
		Calculate the channel matrix
		"""
		self.xtalk_gain = zeroes(self.N,self.N,self.K)
		self.channel_matrix = self.calc_channel_matrix()


	def calc_channel_matrix(self): #TODO
		for K in range(self.K):
			for i,li in enumerate(self.lines):
				for j,lj in enumerate(self.lines):
					if i == j:
						"""
						TODO:
						Could the li.gain assignment be moved inside the object?
						"""
						li.gain = self.xtalk_gain[i,j,k] = li.transfer_fn(freq_on_tone(K))
					else:
						self.xtalk_gain[i,j,k] = self.calc_fext_xtalk_gain(i,j,K,DOWN)



	def freq_on_tone(self,K): #TODO
		"""
		Assume ADSL downstream for now
		"""
		return K * 4312.5 + 140156.25;

	def calc_fext_xtalk_gain(self,victim,xtalker,freq,dir): #TODO
        """
        Calculate Far End XT Gain

        TODO: Look into "Impact of Crosstalk Channel Estimation on the DSM
        Performance for DSL Networks" for an alternative

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


        4 bundle segment types for 2 lines (xtalker and victim);
        H1:xtalker before victim (1,4,7)
        H2:xtalker and victim (all)
        H3:victim after xtalker(1,2,3)
        H4:victim before xtalker(3,6)
		"""




		pass

    def fext(self,freq,length):
        """
        Calculate FEXT for given Freq and Length
        Uses:
            UndB
        """

        """
        Constants taken from
            “Transmission and Multiplexing (TM); Access transmission systems on metallic
            access cables; Very high speed Digital Subscriber Line (VDSL); Part 1: Functional
            requirements,”
        """
        K_xn=0.0032 #NOTUSED
        K_xf=0.0056 #NOTUSED

        """
        Model from BT's simulation parameters document cp38-2
        """
        x=-55
        n=(lines-1)/(max_bundle_size-1) #n disturbers REFACTOR lines max_bundle_size
        x+= 20*log10(freq/90e3)         #f in Hz
        x+= 6*log10(n)                  #shared length in km
        x+= 10*log10(length)

        try:
            return UndB(x)
        except ZeroDivisionError:
            return 0


        """
        Transfer function makes more sense in the bundle as most of the time,
        it is looking at a common length between two lines
        """
    def insertion_loss(self,length,freq): #TODO: Reintroduce factor removed from C version?
        return do_transfer_fn(length,freq )


        """
        Let the transfer function default to type 2;
            Allows for easy default change later
            Allows for easy 'case-based' changes
        """

    def do_transfer_function(self,length,freq,type=2):#TODO Where on earth to Zs/Zl come from? They do nothing!
        """
        Z=impedance/l, Y=admittance/l
        Z=R+jwL, Y=G+jwC
        Z0=Characteristic Impedence, gamma=propagation constant
        Z0=sqrt(Z/Y), gamma=sqrt(Z*Y)
        Should Use PyGSL, but start off with cmath
        """

        Z = complex(self._R(freq),self._L(freq))
        Y = complex(self._G(freq),self._C(freq))

        Z0 = cmath.sqrt(Z/Y)
        gamma = cmath.sqrt(Z*Y)

        Zs = complex(100,0)
        Zl = complex(100,0)

        upper = Z0 * ( 1/cmath.cosh(gamma)) #sech=1/cosh
        lower = Zs * ( (Z0/Zl)+ cmath.tanh(gamma) ) +
                Z0 * ( 1+ (Z0/Zl)*cmath.tanh(gamma) )

        H = upper/lower
        return cmath.polar(H)[0]        #polar returns (abs(h),real(h))

    def _R(self,freq):
		"""
		Return R Parameter for transfer function
		"""
		c_partial = self._R_partial_root(	freq,
											self.material["r_0c"],
											self.material["a_c"]
										)

		if self.material["r_0s"] > 0:
			s_partial = self._R_partial_root(freq,self.material["r_0s"],self.material["a_s"] )
			return (c_partial*s_partial)/(c_partial+s_partial)
		else:
			return c_partial

		return


	def _R_partial_root(self,freq,resistance,constant):
		"""
		Partial root for R
		"""
		return math.pow(constant*math.pow(freq,2)+math.pow(resistance,4),(0.25))

	def _L(self,freq):
		"""
		Return L Parameter for transfer function
		"""
		return (self.material["l_0"]+self.material["l_inf"]*math.pow(self._L_partial(freq),self.material["b"]))/(1+self._L_partial(freq))

	def _L_partial(self,freq):
		return freq*1e-3/self.material["f_m"]

	def _C(self,freq):
		return self.material["c_inf"]+self.material["c_0"]*math.pow(freq,-self.material["c_e"])

	def _G(self,freq):
		return self.material["g_0"]*math.pow(freq,self.material["g_e"])
