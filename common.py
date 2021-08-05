from __future__ import annotations

# Standard python imports
import sys
import os
from datetime import datetime, timedelta, date
import argparse
import configparser
import pathlib

# Non-standard python imports
import yaml

from tqdm import tqdm

from pytz import timezone
from tzlocal import get_localzone

import numpy
from numpy import around

from geopy import distance

#Skyfield imports
from skyfield import almanac
from skyfield.nutationlib import iau2000b
from skyfield.api import Star, EarthSatellite, Loader, iers2010, N, S, E, W


# For typing, not enforced by python; but useful for humans
from skyfield.units import Angle
from skyfield.jpllib import SpiceKernel
from skyfield.toposlib import GeographicPosition
from skyfield.vectorlib import VectorSum
from skyfield.timelib import Time
from skyfield.sgp4lib import EarthSatellite
from skyfield.positionlib import Barycentric

from numpy import array as nparray

from typing import Dict, List, Tuple, Sequence, Optional, Union, Callable, Any

#load into local directory
load = Loader('skyfield_data')


# Script Globals
timescale = load.timescale()
ephemeris = load('de440s.bsp')
earth = ephemeris['earth']
moon = ephemeris['moon']

#Stolen from https://stackoverflow.com/questions/27531718/datetime-timezone-conversion-using-pytz
def timezone_converter(input_dt: datetime, current_tz: Optional[str]='UTC', target_tz: Optional[str]='UTC') -> datetime:
	current_tz = timezone(current_tz)
	target_tz = timezone(target_tz)
	target_dt = None
	if (getattr(input_dt,"tzinfo",None) is not None) and (input_dt.tzinfo.utcoffset(input_dt) is not None):
		target_dt = input_dt.astimezone(target_tz)
	else:
		target_dt = current_tz.localize(input_dt).astimezone(target_tz)
	return target_tz.normalize(target_dt)

# converting between skyfield.Time and datetime to do stuff is annoying and bad design
# object which handles all the time shinnanigans until it is needed for a function input
class timeObj:
	def __init__(self, inputTime: Any) -> None:
		if type(inputTime) == datetime:
			self.dt = timezone_converter(inputTime)
			self.sf = timescale.from_datetime(self.dt)
		elif type(inputTime) == Time:
			self.sf = inputTime
			self.dt = inputTime.utc_datetime()
		else:
			 raise TypeError(f"Only datetime or Skyfield Time objects are allowed!, inputTime was of type {type(inputTime)}")

	def __repr__(self) -> str:
		return "timeObj(datetime=%r, skyfield.Time=%r)" % (self.dt, self.sf)

	def __str__(self) -> str:
		return str(self.dt)

	def getStartDay(self) -> timeObj:
		return timeObj(datetime.combine(self.dt.date(), datetime.min.time()))

	def getOff(self, **kwargs:int) -> timeObj:
		offdt = self.dt + timedelta(**kwargs)
		return timeObj(offdt)

	def getDateTime(self) -> datetime:
		return self.dt

	def getSkyfield(self) -> Time:
		return self.sf

	def isoformat(self) -> str:
		return self.dt.isoformat()

	def getfileformat(self) -> str:
		return self.dt.strftime("%Y%m%d")

# Time Globals
utc = timezone("UTC")
todayUTC = timeObj(datetime.now())

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
def daylength(ephemeris_ : SpiceKernel, topos_ : GeographicPosition, degrees_ : float) -> Callable[[Time], bool]:
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

# ------- Moon position functions -------

# From https://github.com/sczesla/PyAstronomy/blob/master/src/pyasl/asl/angularDistance.py
def twoObjectDistanceExplict(ra1: Union[float, nparray], dec1: Union[float, nparray], ra2: Union[float, nparray], dec2: Union[float, nparray]) -> Union[float, nparray]:
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
	delt_lon = (ra1 - ra2)*numpy.pi/180.
	delt_lat = (dec1 - dec2)*numpy.pi/180.
	# Haversine formula
	dist = 2.0*numpy.arcsin( numpy.sqrt( numpy.sin(delt_lat/2.0)**2 + \
		 numpy.cos(dec1*numpy.pi/180.)*numpy.cos(dec2*numpy.pi/180.)*numpy.sin(delt_lon/2.0)**2 ) )
	return dist/numpy.pi*180.


def moonFraction(inTime: timeObj) -> float:
	 Tt = inTime.getSkyfield()
	 return almanac.fraction_illuminated(ephemeris,"moon",Tt)


def moonPosition(inTime: timeObj, earthObserver: VectorSum) -> Tuple[float, float]:
	Tt = inTime.getSkyfield()
	degreeFunc =  lambda x: x.degrees if type(x) is Angle else x
	astrometric = earthObserver.at(Tt).observe(moon)
	altitudeTmp, azimuthTmp, dist_AU = map(degreeFunc,astrometric.apparent().altaz())
	return (altitudeTmp,azimuthTmp)


def solveForMoonPosition(inTime: timeObj, earthObserver: VectorSum)-> Tuple[float, float, float]:
	moonFrac = moonFraction(inTime)
	moonAlt, moonAz = moonPosition(inTime, earthObserver)
	return (moonFrac, moonAlt, moonAz)


def disFromMoon(inTime: timeObj, currentAlt: float, currentAz: float, earthObserver: VectorSum) -> Tuple[float, float]:
	moonFrac, moonAlt, moonAz = solveForMoonPosition(inTime, earthObserver)
	moonDis = twoObjectDistanceExplict(currentAz, currentAlt, moonAz, moonAlt)
	return (moonDis, moonFrac)



# Counts the number of objects which are not "None" in a list
def nNotNones(lst) -> int:
	return sum(x is not None for x in lst)

# Helper class which stores the transit times and associated information
class satTransit:
	def __init__(self, inAlt_: float, inSatName_: str, inObsLoc_: Tuple[float, float], riseTime_: timeObj, culmTime_: timeObj, setTime_: timeObj) -> None:
		self.alt = inAlt_
		self.satName = inSatName_
		self.obsLoc = inObsLoc_

		self.rise = riseTime_
		self.culm = culmTime_
		self.set = setTime_

	def getDuration(self) -> float:
		return float((self.set.getSkyfield() - self.rise.getSkyfield()) * 24 * 60)

	def getCulmination(self) -> timeObj:
		return self.culm

	def getRise(self) -> timeObj:
		return self.rise

	def getSet(self) -> timeObj:
		return self.set

	def getAlt(self) -> float:
		return self.alt

	def addToDict(self, inDict: Dict[str, Any]) -> None:
		inDict["culm"] = self.getCulmination().getSkyfield()
		inDict["rise"] = self.rise.getSkyfield()
		inDict["set"] = self.set.getSkyfield()
		inDict["dur"] = self.getDuration()
		inDict["rise_set_alt"] = self.alt
		inDict["sat_name"] = self.satName
		inDict["obs_lat"] = self.obsLoc[0]
		inDict["obs_long"] = self.obsLoc[1]


#The calculator class which does all the heavy lifting
class satCalculator:
	def __init__(self, satellite_: EarthSatellite, observatoryName_: str, observatory_: GeographicPosition, earthObservatory_: VectorSum) -> None:
		self.sat = satellite_
		self.obs = observatory_
		self.obsName = observatoryName_
		self.earthobs = earthObservatory_

		self.diff = self.sat - self.obs

	def getSatName(self) -> str:
		return self.sat.name

	def getObsName(self) -> str:
		return self.obsName

	def getEpoch(self) -> Time:
		return self.sat.epoch

	def getObsCoordAt(self, inTime: timeObj) -> Tuple[float, float, float, float]:
		toc = self.diff.at( inTime.getSkyfield() )
		alt, az, _ = toc.altaz()
		ra, dec, _ = toc.radec()
		return (alt._degrees, az._degrees, ra._degrees, dec._degrees)

	def getObsCoordBetween(self, startTime: timeObj, duration: float, step: float, include_altaz: bool=False) -> nparray:
		rv = None
		if include_altaz:
			rv = numpy.zeros((int(duration//step),5), dtype=numpy.float64)
		else:
			rv = numpy.zeros((int(duration//step),3), dtype=numpy.float64)

		#Cacluation step
		iStep = 0
		for tStep in numpy.arange(0, duration, step):
			tStepMicro = int(numpy.floor(tStep * 1000000))
			nextStepT = startTime.getOff(microseconds=tStepMicro)

			alt, az, ra, dec = self.getObsCoordAt(nextStepT)
			rv[iStep][0] = tStep
			rv[iStep][1] = ra
			rv[iStep][2] = dec
			if include_altaz:
				rv[iStep][3] = alt
				rv[iStep][4] = az
			iStep += 1

		return rv

	def getDayTransits(self, inDate: timeObj, input_alt: float) -> Sequence[satTransit]:
		obsLoc = (self.obs.latitude.degrees, self.obs.longitude.degrees)
		inDay = inDate.getStartDay()
		t0 = inDay.getSkyfield()
		t1 = inDay.getOff(days=1).getSkyfield()
		eventTimes, eventList = self.sat.find_events(self.obs, t0, t1, altitude_degrees=input_alt)

		transits = []

		if (len(eventList) >= 3):
				eventTrack = [None, None, None]

				for eventTime, eventIndex in zip(eventTimes, eventList):
					eventTrack[eventIndex] = eventTime
					if nNotNones(eventTrack) == 3:
						eventTrack = [timeObj(x) for x in eventTrack]
						transits.append( satTransit(input_alt, self.getSatName(), obsLoc, *eventTrack) )
						eventTrack = [None, None, None]

				if nNotNones(eventTrack) > 0:
					eventTrack = [timeObj(x) if type(x) is Time else None for x in eventTrack]
					transits.append( satTransit(input_alt, self.getSatName(), obsLoc, *eventTrack) )

		return transits

	def _inDayAtObs(self, inTime: timeObj) -> bool:
		inDay = inTime.getStartDay()
		t0 = inDay.getSkyfield()
		t1 = inDay.getOff(days=1).getSkyfield()
		astro_twil_time, astro_twil_up = almanac.find_discrete(t0, t1, daylength(ephemeris, self.obs, DAYLENGTH_ASTRONOMICAL_TWILIGHT))

		sunRise = None
		sunSet = None
		tArr = astro_twil_time.utc_datetime()
		for time, isRise in zip(tArr, astro_twil_up):
			if isRise:
				sunRise = time
			else:
				sunSet = time

		if ((sunSet < inTime.getDateTime()) and (inTime.getDateTime() < sunRise)):
			return False
		return True


	def _getSatEarthData(self, geocentric_: Barycentric) -> Tuple[Tuple[float, float], Tuple[float, float], float]:
		subpoint = iers2010.subpoint(geocentric_)

		satLoc = (subpoint.latitude.degrees, subpoint.longitude.degrees)
		verLoc = (self.obs.latitude.degrees, self.obs.longitude.degrees)
		greatCDist = distance.distance(satLoc,verLoc)
		return (satLoc, verLoc, greatCDist.km)

	def getSatilliteDataAt(self, inTime: timeObj) -> Dict[str, Any]:
		geocentric = self.sat.at(inTime.getSkyfield())

		retData = {}
		retData["obs_name"] = self.getObsName()


		retData['name'] = self.getSatName()

		retData['sun'] = geocentric.is_sunlit(ephemeris)
		retData['is_day'] = self._inDayAtObs(inTime)


		(retData['lat'], retData['long']), (retData["obs_lat"], retData["obs_long"]), retData['approx_dis'] = self._getSatEarthData(geocentric)

		retData['alt'], retData['az'], retData['ra'], retData['dec'] = self.getObsCoordAt(inTime) # ICRF ("J2000")

		retData['moon_dis'], retData['moon_frac']  = disFromMoon(inTime, retData['alt'], retData['az'], self.earthobs)
		return retData




def buildSatList(reloadSat: bool) -> List[EarthSatellite]:
	global todayUTC
	satURL = "https://celestrak.com/NORAD/elements/active.txt"
	satellites = None
	if reloadSat:
		satellites = load.tle_file(satURL, reload=True)
	else:
		satellites = load.tle_file(satURL)

	with open("config/sat_list.txt", "r") as fp:
		satNames = [x.strip() for x in fp.readlines()]

	satelliteArray = []
	for sat in satellites:
		days = todayUTC.getSkyfield() - sat.epoch
		if (days < 14) and (0 < sum([bool(satname in sat.name) for satname in satNames])):
			satelliteArray.append(sat)

	if len(satelliteArray) == 0:
		print("No inputted satellite from config/sat_list.txt was not found, please retry with --reload or check to ensure each are still active")
		sys.exit(1)

	return satelliteArray


def buildSatCalcs(reloadSat_: bool) -> List[satCalculator]:
	satelliteArray = buildSatList(reloadSat_)

	satCalculators = []

	with open("config/obs_data.yaml","r") as fp:
		obsData = yaml.load(fp, Loader=yaml.SafeLoader)


	for obsName, obsArr in obsData.items():
		observatory = iers2010.latlon(*obsArr)
		earthObservatory = earth + observatory
		for satObject in satelliteArray:
			satCalculators.append(satCalculator(satObject, obsName, observatory, earthObservatory))

	return satCalculators
