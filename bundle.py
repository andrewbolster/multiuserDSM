"""
Class that describes the DSL bundle object
"""

import sys
from line import Line

class Bundle(object):
	def __init__(self,network_file,K):
		self.K = K				# the number of DMT channels
		self.lines = []			# the DSL line objects
		self.xtalk_gain = None	# XT Matrix (NB Must be FULLY instantiated before operating)


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


