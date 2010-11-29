


class Line(object):
	def __init__(self,nt,lt):
        # nt and lt are network termination (CO or RT) and line termination (CPE)
		self.nt = int(nt)
		self.lt = int(lt)
		self.gain = None

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

	def __str__(self):
		"""
		Super cool graphics
		"""
		s = ""
		s += "-" * (self.nt / 250)
		s += "|" 
		s += "-" * ((self.lt-self.nt) / 250)
		s += "|"
		return s
	
	def transfer_fn(self,freq): #TODO
		"""
		Return the RLCG Parameterised Transfer Function
		"""
		pass

	def _R(self,freq):
		"""
		Return R Parameter for transfer function
		"""
		c_partial = self._R_partial_root(	freq,
											self.material["r_0c"],
											self.material["a_c"]
										)

		if self.material["r_s"] > 0:
			s_partial = self._R_partial_root(freq,self.material["r_0s"],self.material["a_s"] )
			return (c_partial*s_partial)/(c_partial+s_partial)
		else:
			return c_partial

		return


	def _R_partial_root(self,freq,resistance,constant):
		"""
		Partial root for R
		TODO: MEMOISE
		"""
		return math.pow(constant*math.pow(freq,2)+math.pow(resistance),(0.25))

	def _L(self,freq):
		"""
		Return L Parameter for transfer function
		"""
		return (self.material["l_0"]+self.material["l_inf"]*math.pow(self._L_partial(freq),self.material["b"]))

	def _L_partial(self,freq):
		return freq/self.material["f_m"]

	def _C(self,freq):
		return self.material["c_inf"]+self.material["c_0"]*math.pow(freq,-self.material["c_e"])

	def _G(self,freq):
		return self.material["g_0"]*math.pow(freq,self.material["g_e"])
		
