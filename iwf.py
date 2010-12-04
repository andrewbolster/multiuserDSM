"""
Iterative Water Filling Algorithm Module
"""
from wx.py.buffer import Buffer

class IWF(algorithm):
	def __init__(self, bundle):
		assert bundle is Bundle, "You dun goofed the Bundle!"
		if bundle.type == "ADSL_DOWNSTREAM":
			p_initial=0.111
		else:
			p_initial=0.02818

		# @type bundle Bundle
		self.p = tile(p_initial, bundle.K) 		#size bundle.K

		utility








