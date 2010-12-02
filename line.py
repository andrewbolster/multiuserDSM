
import utility

class Line(Bundle):
	def __init__(self,nt,lt):
        # nt and lt are network termination (CO or RT) and line termination (CPE)
		self.nt = int(nt)
		self.lt = int(lt)
        self.length = self.nt - self.lt
        self.gain = None
        self.id = none #could this be removed as an array of lines?
        self.type = 2   #I'm assuming this declares the material of the line
                        #so in theory it can be removed from transfer_fn

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

	def transfer_fn(self, freq): #TODO
		"""
		Return the RLCG Parameterised Transfer Function
		"""
        #C Version uses bundle as a linked list of lines so most of this function is redundant.
        return super(Bundle, self).do_transfer_fn(type = self.type, length=self.length/1000, freq )
