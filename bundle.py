"""
Class that describes the DSL bundle object
"""

import sys
from line import Line

class Bundle(object):
	def __init__(self,network_file,K):
		self.K = K			# the number of DMT channels
		self.lines = [] 	# the DSL line objects
		

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
		self.channel_matrix = self.calc_channel_matrix()
		
		
	def calc_channel_matrix(self):
		for K in range(self.K):
			for i,li in enumerate(self.lines):
				for j,lj in enumerate(self.lines):
					if i == j:
						li.transfer_fn(freq_on_tone(K))
						
		
	def freq_on_tone(self,K):
		"""
		Assume ADSL downstream for now
		"""
		return K * 4312.5 + 140156.25;
	
	def calc_fext_xtalk_gain(self,v,x,freq,dir):
		pass
	