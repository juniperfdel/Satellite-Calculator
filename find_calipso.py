from __future__ import annotations

# Standard python imports
import sys
import dateutil
from datetime import datetime, timedelta, date
import argparse

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
from skyfield.api import Star, EarthSatellite, load, iers2010, N, S, E, W

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

# Input Globals
reloadSat = False
obsData = None
nDays = None
localTZ = None


# Script Globals
timescale = load.timescale()
ephemeris = load('de440s.bsp')
earth = ephemeris['earth']
moon = ephemeris['moon']

calipso = None

calipsoNs = []

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
def twoObjectDistanceExplict(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
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

# Helper class which stores the transit times
class calipsoTransit:
	def __init__(self, riseTime_: timeObj, culmTime_: timeObj, setTime_: timeObj) -> None:
		self.rise = riseTime_
		self.culm = culmTime_
		self.set = setTime_

	def getDuration(self) -> float:
		return float((self.set.getSkyfield() - self.rise.getSkyfield()) * 24 * 60)

	def getCulmination(self) -> timeObj:
		return self.culm

	def addToDict(self, inDict: Dict[str, Any]) -> None:
		inDict["culminate"] = self.getCulmination().getSkyfield()
		inDict["rise"] = self.rise.getSkyfield()
		inDict["fall"] = self.set.getSkyfield()
		inDict["dur"] = self.getDuration()


#The calculator class which does all the heavy lifting
class calipsoCalculator:
	def __init__(self, satillite_: EarthSatellite, observatoryName_: str, observatory_: GeographicPosition, earthObservatory_: VectorSum) -> None:
		self.sat = satillite_
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

	def getDayTransits(self, inDate: timeObj, input_alt: float) -> Sequence[calipsoTransit]:
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
						transits.append( calipsoTransit(*eventTrack) )
						eventTrack = [None, None, None]

				if nNotNones(eventTrack) > 0:
					eventTrack = [timeObj(x) if type(x) is Time else None for x in eventTrack]
					transits.append(calipsoTransit(*eventTrack))

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

	def getSatilliteDataAt(self, inTime: timeObj) -> dict[str, Any]:
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




def turnToCSV(inFileName: str, inData:  Dict[str, Any]) -> None:
	global localTZ
	localAbbr = datetime.now().astimezone(localTZ).strftime('%Z')
	csvRowDefs = [
		("Observatory Latitude", lambda x: x["obs_lat"]),
		("Observatory Longitude", lambda x: x["obs_long"]),
		("Satellite name", lambda x: x["name"]),
		("Culmination Time (UTC)", lambda x: x["culminate"].astimezone(utc).strftime('%Y %b %d %H:%M:%S')),
		(f"Culmination Time ({localAbbr})", lambda x: x["culminate"].astimezone(localTZ).strftime('%Y %b %d %H:%M:%S')),
		("Daytime at Culmination?", lambda x: 'Y' if x['is_day'] else 'N'),
		("Satellite Shadow Latitude (At Culmination)", lambda x: x["lat"]),
		("Satellite Shadow Longitude (At Culmination)", lambda x: x["long"]),
		("Approx. Distance Between Shadow and Observatory (km)", lambda x: x["approx_dis"]),
		("Satellite RA (At Culmination)", lambda x: x["ra"]),
		("Satellite Dec (At Culmination)", lambda x: x["dec"]),
		("Satellite Altitude (At Culmination)", lambda x: x["alt"]),
		("Satellite Azimuth (At Culmination)", lambda x: x["az"]),
		("Duration Above 30 deg alt [minutes]", lambda x: around(x['dur'],2)),
		("Is the Satellite Sunlit?", lambda x: 'Y' if x['sun'] else 'N'),
		("Moon Distance (At Culmination)", lambda x: around(x['moon_dis'],2)),
		("Moon Fraction (At Culmination)", lambda x: around(x['moon_frac'],3))
	]

	csvHead = ",".join([x[0] for x in csvRowDefs]) + "\n"
	csvTransforms = [x[1] for x in csvRowDefs]
	with open(inFileName,"w") as fp:
		fp.write(csvHead)
		for data in inData:
			rowS = ",".join([str(csvTrans(data)) for csvTrans in csvTransforms]) + "\n"
			fp.write(rowS)


def main(calipsoNs: Sequence[calipsoCalculator], nDays: int, solveForSatPos: Union[Sequence[float], None]) -> None:
	global todayUTC

	finalDay = todayUTC.getOff(days=nDays)

	print("Calculating CALIPSO culminations")

	for calipsoCalc in calipsoNs:
		print(f"Looking for CALIPSO at {calipsoCalc.getObsName()} from {todayUTC} to {finalDay}")

		if solveForSatPos is not None:
			satStore = h5py.File(calipsoCalc.getObsName() + '.hdf5', 'w')

		culmData = []
		for dayNum in tqdm(range(nDays)):

			nextDay = todayUTC.getOff(days=dayNum)
			dayTransits = calipsoCalc.getDayTransits(nextDay, 30.)

			for dayTransit in dayTransits:

					culmination = dayTransit.getCulmination()
					satData = calipsoCalc.getSatilliteDataAt(culmination)
					dayTransit.addToDict(satData)
					culmData.append(satData)

					if solveForSatPos is not None:
							microsecondsBefore = solveForSatPos[0] * -1000000
							secondsDur = solveForSatPos[1]
							secondsStep = solveForSatPos[2]

							startT = culmination.getOff(microseconds=microsecondsBefore)
							storePos = calipsoCalc.getObsCoordBetween(startT, secondsDur, secondsStep, bool(int(solveForSatPos[3])))

							#Storage step
							datasetName = startT.isoformat()
							satStore.create_dataset(datasetName, data=storePos)

		if len(culmData) > 0:
			print(f"Found CALIPSO culminating {len(culmData)} times at {calipsoCalc.getObsName()}, for more info open {calipsoCalc.getObsName()}.csv")
			print("")
			turnToCSV(calipsoCalc.getObsName() + ".csv", culmData)
		else:
			print(f"Did not find CALIPSO culminating at {obsName} from {todayUTC} to {finalDay}")

		if solveForSatPos:
			satStore.close()

def buildSatGlobals(reloadSat: bool) -> None:
	global calipso, todayUTC
	satURL = "https://celestrak.com/NORAD/elements/active.txt"
	satellites = None
	if reloadSat:
		satellites = load.tle_file(satURL, reload=True)
	else:
		satellites = load.tle_file(satURL)

	for sat in satellites:
		days = todayUTC.getSkyfield() - sat.epoch
		if (days < 14) and ("CALIPSO" in sat.name):
			calipso = sat
			break

	if calipso is None:
		print("CALIPSO was not found, please retry with --reload or check NASA to ensure it is still active")
		sys.exit(1)



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Calculate CALIPSO culminations for observatorys found in obs_data.yaml')
	parser.add_argument('days', type=int, nargs='?', default=120, help='The number of days to calculate for, default of 120 days')
	parser.add_argument('-tz','--set-timezone', type=str, help="The timezone to cacluate with respect to, default is your local timezone")
	parser.add_argument('-r','--reload', action='store_true', help="Re-download the TLE file from the internet, default is false")
	parser.add_argument('-s','--solve', nargs=4, type=float, help="Solve for the satillite positions around culmination, Format: <start time in seconds before culmination> <duration in seconds> <timestep in seconds> <alt/az flag>. If <alt/az flag> is set to 1, then alt/az data will be added. Setting any other number will prevent this.")

	args = parser.parse_args()
	print(args)

	nDays = args.days
	reloadSat = args.reload

	inputTZ = getattr(args,"set_timezone", None)
	if inputTZ is None:
		localTZ = get_localzone()
	else:
		localTZ = timezone(inputTZ)

	solveForSatPos = getattr(args,"solve", None)

	if solveForSatPos is not None:
		try:
			import h5py
		except ImportError as err:
			print("In order to store satillite positons please install h5py https://pypi.org/project/h5py/\n\n")
			raise err


	buildSatGlobals(reloadSat)

	with open("obs_data.yaml","r") as fp:
		obsData = yaml.load(fp, Loader=yaml.SafeLoader)

	for obsName, obsArr in obsData.items():
		observatory = iers2010.latlon(*obsArr)
		earthObservatory = earth + observatory
		calipsoNs.append(calipsoCalculator(calipso, obsName, observatory, earthObservatory))

	main(calipsoNs, nDays, solveForSatPos)
