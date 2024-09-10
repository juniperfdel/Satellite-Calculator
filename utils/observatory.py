from typing import Tuple

from skyfield.toposlib import GeographicPosition, iers2010
from skyfield.vectorlib import VectorSum

from utils.skyfield_utils import SkyfieldConstants


class Observatories:
    def __init__(self, lat: float, long: float, alt: float, in_name: str):
        self.name = in_name
        self.lat = lat
        self.long = long
        self.alt = alt
        self.sf: GeographicPosition = iers2010.latlon(self.lat, self.long, elevation_m=float(self.alt))
        self.sf_vector_sum: VectorSum = SkyfieldConstants.earth + self.sf
        self.lat_long = (self.lat, self.long)
        self.values = [self.lat, self.long, self.alt]

    def __iter__(self):
        yield from self.values
