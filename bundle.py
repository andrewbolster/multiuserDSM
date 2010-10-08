"""
Class that describes the DSL bundle object
"""

import sys
from line import Line

class Bundle(object):
	def __init__(self,network_file):

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


		for line in self.lines:
			print line
