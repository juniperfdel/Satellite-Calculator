from typing import Union

import numpy
from numpy import array


def two_object_distance(
		ra1: Union[float, array], dec1: Union[float, array],
		ra2: Union[float, array], dec2: Union[float, array]
) -> Union[float, array]:
	"""
	Calculate the angular distance between two coordinates.

	Parameters
	----------
	ra1 : float, array
		Right ascension of the first object in degrees.
	dec1 : float, array
		Declination of the first object in degrees.
	ra2 : float, array
		Right ascension of the second object in degrees.
	dec2 : float, array
		Declination of the second object in degrees.
	Returns
	-------
	Angle : float, array
		The angular distance in DEGREES between the first
		and second coordinate in the sky.

	"""
	delta_lon = numpy.deg2rad(ra1 - ra2)
	delta_lat = numpy.deg2rad(dec1 - dec2)
	# Haversine formula
	dist = \
		2.0 * numpy.arcsin(
			numpy.sqrt(
				numpy.sin(delta_lat / 2.0) ** 2 +
				numpy.cos(dec1 * numpy.pi / 180.) *
				numpy.cos(dec2 * numpy.pi / 180.) *
				numpy.sin(delta_lon / 2.0) ** 2
			)
		)
	return numpy.rad2deg(dist)
