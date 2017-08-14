import numpy as np
import scipy as scp

class SupernovaModel(object, ):
	
	def __init__(self, ):
		pass
		
	def SignalTimeFunction(self, t, ):
		value = (1.+t/.1)**-1.
		return value	
	
	def SignalTimePDF(self, ):
		norm = scp.integrate.quad(self.SignalTimeFunction, 0., 1.)[0]
		f = lambda t: self.SignalTimeFunction(t)/norm
		return f	
	
	def BackgroundTimePDF(self, ):
		f = lambda t: np.ones_like(t) / 1.
		return f