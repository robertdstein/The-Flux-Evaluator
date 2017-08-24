import numpy as np
from icecube import astro	


class Sphere:
    """Class for calculating great circle distances.

    """

    def __init__(self, ):
        pass

    def great_circle_distance(self, ra_1, dec_1, ra_2, dec_2, ):
        """Calculates Great Circle Distance between two points on the
        celestial sphere.

        :param ra_1: First Right Ascension
        :param dec_1: First Declination
        :param ra_2: Second Right Ascension
        :param dec_2: Second Declination
        :return: Great Circle distance
        """
        value = astro.angular_distance(ra_1, dec_1, ra_2, dec_2)
        return value
