


class Line(object):
	def __init__(self,nt,lt):
		self.nt = int(nt)
		self.lt = int(lt)

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
