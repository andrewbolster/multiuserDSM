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
					
		except IOError:
			print "Cannot open the network file",network_file
			sys.exit(1)


		"""
		Calculate the channel matrix
		"""
		self.channel_matrix = self.calc_channel_matrix()
		
		
	def calc_channel_matrix(self):
		pass
		
	def transfer_fn(self,line_id,freq):
		pass
	
	def calc_fext_xtalk_gain(self,v,x,freq,dir):
		pass
	