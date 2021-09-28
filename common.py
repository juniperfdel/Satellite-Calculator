from __future__ import annotations

# Standard python imports
import sys
from datetime import datetime, timedelta
from typing import Dict, List, Tuple, Sequence, Optional, Union, Callable, Any

import numpy
# Non-standard python imports
import yaml
from geopy import distance
from numpy import array as nparray
from pytz import timezone
# Skyfield imports
from skyfield import almanac
from skyfield.api import Loader, iers2010
from skyfield.jpllib import SpiceKernel
from skyfield.nutationlib import iau2000b
from skyfield.positionlib import Barycentric
from skyfield.sgp4lib import EarthSatellite
from skyfield.timelib import Time
from skyfield.toposlib import GeographicPosition
# For typing, not enforced by python; but useful for humans
from skyfield.units import Angle
from skyfield.vectorlib import VectorSum

# load into local directory
load = Loader('skyfield_data')

# Script Globals
timescale = load.timescale()
ephemeris = load('de440s.bsp')
earth = ephemeris['earth']
moon = ephemeris['moon']


# Stolen from https://stackoverflow.com/questions/27531718/datetime-timezone-conversion-using-pytz
def timezone_converter(
		input_dt: datetime,
		current_tz: Optional[str] = 'UTC',
		target_tz: Optional[str] = 'UTC'
) -> datetime:
	current_tz = timezone(current_tz)
	target_tz = timezone(target_tz)
	target_dt = None
	if (getattr(input_dt, "tzinfo", None) is not None) and (input_dt.tzinfo.utcoffset(input_dt) is not None):
		target_dt = input_dt.astimezone(target_tz)
	else:
		target_dt = current_tz.localize(input_dt).astimezone(target_tz)
	return target_tz.normalize(target_dt)

# converting between skyfield.Time, datetime, and numpy.datetime to do stuff is annoying and bad design
# object which handles all the time shenanigans until it is needed for a function input
class TimeObj:
	def __init__(self, input_time: Any) -> None:
		if isinstance(input_time, datetime):
			self.dt = timezone_converter(input_time)
			self.sf = timescale.from_datetime(self.dt)
			self.np = numpy.datetime64(self.dt)
		elif isinstance(input_time, Time):
			self.sf = input_time
			self.dt = input_time.utc_datetime()
			self.np = numpy.datetime64(self.dt)
		else:
			raise TypeError(
				"Only datetime or Skyfield Time objects are allowed!, "
				f"inputTime was of type {type(input_time)}"
			)
	
	def __repr__(self) -> str:
		return "timeObj(datetime=%r, skyfield.Time=%r)" % (self.dt, self.sf)
	
	def __str__(self) -> str:
		return str(self.dt)
	
	def get_start_day(self) -> TimeObj:
		return TimeObj(datetime.combine(self.dt.date(), datetime.min.time()))
	
	def get_off(self, **kwargs: int) -> TimeObj:
		return TimeObj(self.dt + timedelta(**kwargs))
	
	def get_datetime(self) -> datetime:
		return self.dt
	
	def get_skyfield(self) -> Time:
		return self.sf
	
	def get_numpy(self) -> numpy.datetime64:
		return self.np
	
	def iso_format(self) -> str:
		return self.dt.isoformat()
	
	def get_file_format(self) -> str:
		return self.dt.strftime("%Y%m%d")
	
	def get_off_list(self, t_step: numpy.timedelta64, t_end: TimeObj, final_type: int = 0):
		num_arr = numpy.arange(self.np, t_end.np, t_step)
		if final_type == 0:
			return num_arr
		elif final_type == 1:
			return [timezone_converter(t) for t in num_arr.astype(datetime)]
		elif final_type == 2:
			dt_list = [timezone_converter(t) for t in num_arr.astype(datetime)]
			return timescale.from_datetimes(dt_list)
		else:
			raise TypeError("Please select a type 0 (numpy), 1 (datetime), or 2 (Time)! ")

# Time Globals
utc = timezone("UTC")

# ------- Sun rise/set functions -------

# DAYLENGTH parameters specifing the distance below
# the horizon the sun needs to be for the event to be considered
# occured; the values are from
# https://www.weather.gov/fsd/twilight
# https://ase.tufts.edu/cosmos/view_picture.asp?id=1331
DAYLENGTH_CENTER_HORIZON = 0.0
DAYLENGTH_TOP_HORIZON = 0.26667
DAYLENGTH_TOP_HORIZON_APPARENTLY = 0.8333
DAYLENGTH_CIVIL_TWILIGHT = 6.0
DAYLENGTH_NAUTICAL_TWILIGHT = 12.0
DAYLENGTH_ASTRONOMICAL_TWILIGHT = 18.0


# a helper function for finding the sun set/rise times with skyfield
# from https://stackoverflow.com/questions/54485777/finding-twilight-times-with-skyfield
def daylength(ephemeris_: SpiceKernel, topos_: GeographicPosition, degrees_: float) -> Callable[[Time], bool]:
	"""Build a function of time that returns the daylength.

	The function that this returns will expect a single argument that is a
	:class:`~skyfield.timelib.Time` and will return ``True`` if the sun is up
	or twilight has started, else ``False``.
	"""
	sun = ephemeris_['sun']
	topos_at = (ephemeris_['earth'] + topos_).at
	
	def is_sun_up_at(t: Time):
		"""Return `True` if the sun has risen by time `t`."""
		t._nutation_angles = iau2000b(t.tt)
		return topos_at(t).observe(sun).apparent().altaz()[0].degrees > -degrees_
	
	is_sun_up_at.rough_period = 0.5  # twice a day
	return is_sun_up_at


# From https://github.com/sczesla/PyAstronomy/blob/master/src/pyasl/asl/angularDistance.py
def two_object_distance_explict(
		ra1: Union[float, nparray], dec1: Union[float, nparray],
		ra2: Union[float, nparray], dec2: Union[float, nparray]
) -> Union[float, nparray]:
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
	delt_lon = (ra1 - ra2) * numpy.pi / 180.
	delt_lat = (dec1 - dec2) * numpy.pi / 180.
	# Haversine formula
	dist = \
		2.0 * numpy.arcsin(
			numpy.sqrt(
				numpy.sin(delt_lat / 2.0) ** 2 +
				numpy.cos(dec1 * numpy.pi / 180.) *
				numpy.cos(dec2 * numpy.pi / 180.) *
				numpy.sin(delt_lon / 2.0) ** 2
			)
		)
	return dist / numpy.pi * 180.


# ------- Moon functions -------

def moon_fraction(in_time: TimeObj) -> float:
	tt = in_time.get_skyfield()
	return almanac.fraction_illuminated(ephemeris, "moon", tt)


def convert_to_deg(in_object: Any) -> Any:
	return in_object.degrees if type(in_object) is Angle else in_object


def moon_position(in_time: TimeObj, earth_observer: VectorSum) -> Tuple[float, float]:
	tt = in_time.get_skyfield()
	astrometric = earth_observer.at(tt).observe(moon)
	altitude_tmp, azimuth_tmp, dist_au = map(convert_to_deg, astrometric.apparent().altaz())
	return altitude_tmp, azimuth_tmp


def solve_for_moon_position(in_time: TimeObj, earth_observer: VectorSum) -> Tuple[float, float, float]:
	moon_frac = moon_fraction(in_time)
	moon_alt, moon_az = moon_position(in_time, earth_observer)
	return moon_frac, moon_alt, moon_az


def dis_from_moon(in_time: TimeObj, current_alt: float, current_az: float, earth_observer: VectorSum) -> Tuple[
	float, float]:
	moon_frac, moon_alt, moon_az = solve_for_moon_position(in_time, earth_observer)
	moon_dis = two_object_distance_explict(current_az, current_alt, moon_az, moon_alt)
	return moon_dis, moon_frac


# Counts the number of objects which are not "None" in a list
def num_not_none(lst) -> int:
	return sum(x is not None for x in lst)


# Helper class which stores the transit times and associated information
class SatTransit:
	def __init__(
			self,
			in_alt_: float,
			in_sat_name_: str,
			in_obs_loc_: Tuple[float, float],
			rise_time_: TimeObj,
			culm_time_: TimeObj,
			set_time_: TimeObj
	) -> None:
		self.alt = in_alt_
		self.satName = in_sat_name_
		self.obsLoc = in_obs_loc_
		
		self.rise = rise_time_
		self.culm = culm_time_
		self.set = set_time_
	
	def get_duration(self) -> float:
		return float((self.set.get_skyfield() - self.rise.get_skyfield()) * 24 * 60)
	
	def get_culmination(self) -> TimeObj:
		return self.culm
	
	def get_rise(self) -> TimeObj:
		return self.rise
	
	def get_set(self) -> TimeObj:
		return self.set
	
	def get_alt(self) -> float:
		return self.alt
	
	def add_to_dict(self, in_dict: Dict[str, Any]) -> None:
		in_dict["culm"] = self.get_culmination().get_skyfield()
		in_dict["rise"] = self.rise.get_skyfield()
		in_dict["set"] = self.set.get_skyfield()
		in_dict["dur"] = self.get_duration()
		in_dict["rise_set_alt"] = self.alt
		in_dict["sat_name"] = self.satName
		in_dict["obs_lat"] = self.obsLoc[0]
		in_dict["obs_long"] = self.obsLoc[1]


# The calculator class which does all the heavy lifting
class SatCalculator:
	def __init__(
			self,
			satellite_: EarthSatellite,
			observatory_name_: str,
			observatory_: GeographicPosition,
			earth_observatory_: VectorSum
	) -> None:
		
		self.sat = satellite_
		self.obs = observatory_
		self.obsName = observatory_name_
		self.earthobs = earth_observatory_
		
		self.diff = self.sat - self.obs
	
	def get_sat_name(self) -> str:
		return self.sat.name.replace("/", "|")
	
	def get_obs_name(self) -> str:
		return self.obsName
	
	def get_epoch(self) -> Time:
		return self.sat.epoch
	
	def get_obs_coord_at(self, in_time: TimeObj) -> Tuple[float, float, float, float]:
		toc = self.diff.at(in_time.get_skyfield())
		alt, az, _ = toc.altaz()
		ra, dec, _ = toc.radec()
		return alt._degrees, az._degrees, ra._degrees, dec._degrees
	
	def get_obs_coord_between(
			self,
			start_time: TimeObj,
			duration: float,
			step: float,
			include_altaz: bool = False
	) -> nparray:
		
		rv = None
		if include_altaz:
			rv = numpy.zeros((int(duration // step), 5), dtype=numpy.float64)
		else:
			rv = numpy.zeros((int(duration // step), 3), dtype=numpy.float64)
		
		us_step = int(step * 1000000)
		us_duration = int(duration * 1000000)
		
		np_step = numpy.timedelta64(us_step, 'us')
		np_start_time = numpy.datetime64(start_time.get_datetime())
		
		end_time = start_time.get_off(microseconds=us_duration)
		
		np_time_list = start_time.get_off_list(np_step, end_time)
		
		sf_time_list = timescale.from_datetimes([timezone_converter(t) for t in np_time_list.astype(datetime)])
		toc_list = self.diff.at(sf_time_list)
		ra, dec, _ = toc_list.radec()
		rv[:, 0] = (np_time_list - np_start_time).astype(numpy.float64) / 1000000.0
		rv[:, 1] = ra._degrees
		rv[:, 2] = dec._degrees
		if include_altaz:
			alt, az, _ = toc_list.altaz()
			rv[:, 3] = alt._degrees
			rv[:, 4] = az._degrees
		return rv
	
	def get_day_transits(self, in_date: TimeObj, input_alt: float) -> Sequence[SatTransit]:
		obs_loc = (self.obs.latitude.degrees, self.obs.longitude.degrees)
		in_day = in_date.get_start_day()
		t0 = in_day.get_skyfield()
		t1 = in_day.get_off(days=1).get_skyfield()
		event_times, event_list = self.sat.find_events(self.obs, t0, t1, altitude_degrees=input_alt)
		
		transits = []
		
		if len(event_list) >= 3:
			event_track = [None, None, None]
			
			for eventTime, eventIndex in zip(event_times, event_list):
				event_track[eventIndex] = eventTime
				if num_not_none(event_track) == 3:
					event_track = [TimeObj(x) for x in event_track]
					transits.append(SatTransit(input_alt, self.get_sat_name(), obs_loc, *event_track))
					event_track = [None, None, None]
			
			if num_not_none(event_track) > 0:
				event_track = [TimeObj(x) if type(x) is Time else None for x in event_track]
				transits.append(SatTransit(input_alt, self.get_sat_name(), obs_loc, *event_track))
		
		return transits
	
	def _in_day_at_obs(self, in_time: TimeObj) -> bool:
		in_day = in_time.get_start_day()
		t0 = in_day.get_skyfield()
		t1 = in_day.get_off(days=1).get_skyfield()
		astro_twil_time, astro_twil_up = \
			almanac.find_discrete(t0, t1, daylength(ephemeris, self.obs, DAYLENGTH_ASTRONOMICAL_TWILIGHT))
		
		sun_rise = None
		sun_set = None
		t_arr = astro_twil_time.utc_datetime()
		for time, isRise in zip(t_arr, astro_twil_up):
			if isRise:
				sun_rise = time
			else:
				sun_set = time
		
		return sun_set >= in_time.get_datetime() or in_time.get_datetime() >= sun_rise
	
	def _get_sat_earth_data(self, geocentric_: Barycentric) -> Tuple[Tuple[float, float], Tuple[float, float], float]:
		subpoint = iers2010.subpoint(geocentric_)
		
		sat_loc = (subpoint.latitude.degrees, subpoint.longitude.degrees)
		obs_loc = (self.obs.latitude.degrees, self.obs.longitude.degrees)
		great_c_dist = distance.distance(sat_loc, obs_loc)
		return sat_loc, obs_loc, great_c_dist.km
	
	def get_satillite_data_at(self, in_time: TimeObj) -> Dict[str, Any]:
		geocentric = self.sat.at(in_time.get_skyfield())
		
		ret_data = {
			"obs_name": self.get_obs_name(),
			"name": self.get_sat_name(),
			"sun": geocentric.is_sunlit(ephemeris),
			"is_day": self._in_day_at_obs(in_time)
		}
		
		(ret_data['lat'], ret_data['long']), (ret_data["obs_lat"], ret_data["obs_long"]), ret_data[
			'approx_dis'] = self._get_sat_earth_data(geocentric)
		
		ret_data['alt'], ret_data['az'], ret_data['ra'], ret_data['dec'] = self.get_obs_coord_at(
			in_time)  # ICRF ("J2000")
		
		ret_data['moon_dis'], ret_data['moon_frac'] = \
			dis_from_moon(in_time, ret_data['alt'], ret_data['az'], self.earthobs)
		
		return ret_data


def build_sat_list(reload_sat: bool, use_all_: bool = False) -> List[EarthSatellite]:
	today_utc = TimeObj(datetime.now())
	sat_url = "https://celestrak.com/NORAD/elements/active.txt"
	satellites = None
	if reload_sat:
		satellites = load.tle_file(sat_url, reload=True)
	else:
		satellites = load.tle_file(sat_url)
	
	with open("config/sat_list.txt", "r") as fp:
		sat_names = [x.strip() for x in fp.readlines()]
	
	satellite_array = []
	for sat in satellites:
		days = today_utc.get_skyfield() - sat.epoch
		if days < 14 and (sum(bool(satname in sat.name) for satname in sat_names) > 0 or use_all_):
			satellite_array.append(sat)
	
	if not satellite_array:
		print(
			"No inputted satellite from config/sat_list.txt was not found, "
			"please retry with --reload or check to ensure each are still active"
		)
		sys.exit(1)
	
	return satellite_array


def build_sat_calcs(reload_sat_: bool, use_all_: bool = False) -> List[SatCalculator]:
	satellite_array = build_sat_list(reload_sat_, use_all_)
	
	sat_calculators = []
	
	with open("config/obs_data.yaml", "r") as fp:
		obs_data = yaml.load(fp, Loader=yaml.SafeLoader)
	
	for obsName, obsArr in obs_data.items():
		observatory = iers2010.latlon(*obsArr)
		earth_observatory = earth + observatory
		for satObject in satellite_array:
			sat_calculators.append(SatCalculator(satObject, obsName, observatory, earth_observatory))
	
	return sat_calculators
