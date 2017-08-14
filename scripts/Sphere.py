import numpy as np
from icecube import astro	


class Sphere():
	
	def __init__(self, ):
		pass    
	
	def great_circle_distance(self, ra_1, dec_1, ra_2, dec_2, ):
		value = astro.angular_distance(ra_1, dec_1, ra_2, dec_2)
		return value