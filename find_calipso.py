import sys
import dateutil
import numpy
from tqdm import tqdm

from datetime import datetime, timedelta, date
from pytz import timezone
from tzlocal import get_localzone
from numpy import around

from skyfield import almanac
from skyfield.nutationlib import iau2000b
from skyfield.api import Star, EarthSatellite, load, iers2010, N, S, E, W
from skyfield.units import Angle as skyAngle

from geopy import distance


# Input Globals
reloadSat = False
obsData = None
nDays = None
localTZ = None


# Script Globals
startTime = None
timescale = None
ephemeris = None
utc = None

calipso = None

earth = None
moon = None

todayUTC = None

satStore = None

calipsoNs = []

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
def daylength(ephemeris_, topos_, degrees_):
	"""Build a function of time that returns the daylength.

	The function that this returns will expect a single argument that is a
	:class:`~skyfield.timelib.Time` and will return ``True`` if the sun is up
	or twilight has started, else ``False``.
	"""
	sun = ephemeris_['sun']
	topos_at = (ephemeris_['earth'] + topos_).at

	def is_sun_up_at(t):
		"""Return `True` if the sun has risen by time `t`."""
		t._nutation_angles = iau2000b(t.tt)
		return topos_at(t).observe(sun).apparent().altaz()[0].degrees > -degrees_

	is_sun_up_at.rough_period = 0.5  # twice a day
	return is_sun_up_at

# ------- Moon position functions -------

# From https://github.com/sczesla/PyAstronomy/blob/master/src/pyasl/asl/angularDistance.py
def twoObjectDistanceExplict(ra1,dec1,ra2,dec2):
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


def moonFraction(dTt):
	 Tt = timescale.utc(dTt)
	 return almanac.fraction_illuminated(ephemeris,"moon",Tt)


def moonPosition(dTt, earthObserver):
	Tt = timescale.utc(dTt)
	degreeFunc =  lambda x: x.degrees if type(x) is skyAngle else x
	astrometric = earthObserver.at(Tt).observe(moon)
	altitudeTmp, azimuthTmp, dist_AU = map(degreeFunc,astrometric.apparent().altaz())
	return (altitudeTmp,azimuthTmp)


def solveForMoonPosition(currentDT, earthObserver):
	moonFrac = moonFraction(currentDT)
	moonAlt, moonAz = moonPosition(currentDT, earthObserver)
	return (moonFrac, moonAlt, moonAz)


def disFromMoon(currentDT, currentAlt, currentAz, earthObserver):
	degreeFunc =  lambda x: x.degrees if type(x) is skyAngle else x
	currentAlt = degreeFunc(currentAlt)
	currentAz = degreeFunc(currentAz)
	moonFrac, moonAlt, moonAz = solveForMoonPosition(currentDT, earthObserver)
	moonDis = twoObjectDistanceExplict(currentAz, currentAlt, moonAz, moonAlt)
	return moonDis, moonFrac



# Counts the number of objects which are not "None" in a list
def nNotNones(lst):
	return sum(x is not None for x in lst)

#Stolen from https://stackoverflow.com/questions/27531718/datetime-timezone-conversion-using-pytz
def timezone_converter(input_dt, current_tz='UTC', target_tz='UTC'):
	current_tz = timezone(current_tz)
	target_tz = timezone(target_tz)
	target_dt = None
	if (getattr(input_dt,"tzinfo",None) is not None) and (input_dt.tzinfo.utcoffset(input_dt) is not None):
		target_dt = input_dt.astimezone(target_tz)
	else:
		target_dt = current_tz.localize(input_dt).astimezone(target_tz)
	return target_tz.normalize(target_dt)


class timeObj:
	def __init__(self, inputTime):
		if type(inputTime) == datetime:
			self.dt = timezone_converter(inputTime)
			self.sf = timescale.from_datetime(self.dt)
		else:
			self.sf = inputTime
			self.dt = inputTime.utc_datetime()

	def getStartDay(self):
		return timeObj(datetime.combine(self.dt.date(), startTime))

	def getOff(self, **kwargs):
		offdt = self.dt + timedelta(**kwargs)
		return timeObj(offdt)

	def getDateTime(self):
		return self.dt

	def getSkyfield(self):
		return self.sf

# Helper class which stores the transit times
class calipsoTransit:
	def __init__(self, riseTime_, culmTime_, setTime_):
		self.rise = riseTime_
		self.culm = culmTime_
		self.set = setTime_

	def getDuration(self):
		return (self.set - self.rise) * 24 * 60

	def getCulmination(self):
		return self.culm

	def addToDict(self, inDict):
		inDict["culminate"] = self.getCulmination()
		inDict["rise"] = self.rise
		inDict["fall"] = self.set


#The calculator class which does all the heavy lifting
class calipsoCalculator:
	def __init__(self, satillite_, observatoryName_, observatory_, earthObservatory_):
		self.sat = satillite_
		self.obs = observatory_
		self.obsName = observatoryName_
		self.earthobs = earthObservatory_

		self.diff = self.sat - self.obs

	def getSatName(self):
		return self.sat.name

	def getObsName(self):
		return self.obsName

	def getEpoch(self):
		return self.sat.epoch

	def getObsCoordAt(self, inTime):
		inTimesf = ts.from_datetime(inTime)
		toc = self.diff.at(inTimesf)
		alt, az, _ = toc.altaz()
		ra, dec, _ = toc.radec()
		return (alt, az, ra, dec)

	def getObsCoordBetween(self, startTime, duration, step, include_altaz=False):
		rv = None
		if include_altaz:
			rv = numpy.zeros((duration//step,5), dtype=numpy.float64)
		else:
			rv = numpy.zeros((duration//step,3), dtype=numpy.float64)

		#Cacluation step
		iStep = 0
		for tStep in numpy.arange(0, duration, step):
			tStepMicro = int(numpy.floor(tStep * 1000000))
			nextStepT = startTime + timedelta(microseconds=tStepMicro)

			alt, az, ra, dec = self.getObsCoordAt(nextStepT)
			rv[iStep][0] = tStep
			rv[iStep][1] = ra._degrees
			rv[iStep][2] = dec._degrees
			if include_altaz:
				rv[iStep][3] = alt._degrees
				rv[iStep][4] = az._degrees
			iStep += 1

		return rv

	def getDayTransits(self, inDate, input_alt):
		inDay = datetime.combine(inDate.date(), startTime)
		inDay = timezone_converter(inDay)
		t0 = timescale.utc(inDay)
		t1 = timescale.utc(inDay + timedelta(days=1))
		eventTimes, eventList = self.sat.find_events(self.obs, t0, t1, altitude_degrees=input_alt)

		transits = []

		if (len(eventList) >= 3):
				eventTrack = [None, None, None]

				for eventTime, eventIndex in zip(eventTimes, eventList):
					eventTrack[eventIndex] = eventTime
					if nNotNones(eventTrack) == 3:
						transits.append( calipsoTransit(*eventTrack) )
						eventTrack = [None, None, None]

				if nNotNones(eventTrack) > 0:
					transits.append(calipsoTransit(*eventTrack))

		return transits

	def _inDayAtObs(self, inTime):
		inDate = datetime.combine(inTime.date(), startTime)
		inDate = timezone_converter(inDate)
		t0 = timescale.utc(inDate)
		t1 = timescale.utc(inDate + timedelta(days=1))
		astro_twil_time, astro_twil_up = almanac.find_discrete(t0, t1, daylength(ephemeris, self.obs, DAYLENGTH_ASTRONOMICAL_TWILIGHT))

		sunRise = None
		sunSet = None
		tArr = astro_twil_time.utc_datetime()
		for time, isRise in zip(tArr, astro_twil_up):
			if isRise:
				sunRise = time
			else:
				sunSet = time

		if ((sunSet < inTime) and (inTime < sunRise)):
			return False
		return True


	def _getSatEarthData(self, geocentric_):
		subpoint = iers2010.subpoint(geocentric_)

		satLoc = (subpoint.latitude.degrees, subpoint.longitude.degrees)
		verLoc = (self.obs.latitude.degrees, self.obs.longitude.degrees)
		greatCDist = distance.distance(satLoc,verLoc)
		return satLoc, verLoc, greatCDist.km

	def getSatilliteDataAt(self, inTime):
		inTime = timezone_converter(inTime)
		inTimesf = ts.from_datetime(inTime)
		geocentric = self.sat.at(inTimesf)

		retData = {}
		retData["obs_name"] = self.getObsName()


		retData['name'] = self.getSatName()

		retData['sun'] = geocentric.is_sunlit(ephemeris)
		retData['is_day'] = self._inDayAtObs(inTime)

		retData['moon_dis'], retData['moon_frac']  = disFromMoon(inTime, retData['alt'], retData['az'], self.earthobs)

		(retData['lat'], retData['long']), (retData["obs_lat"], retData["obs_long"]), satData['approx_dis'] = self._getSatEarthData(geocentric)

		retData['alt'], retData['az'], retData['ra'], retData['dec'] = self.getObsCoordAt(inTime) # ICRF ("J2000")

		return retData




def turnToCSV(inFileName, inData):
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
		("Satellite RA (At Culmination)", lambda x: x["ra"]._degrees),
		("Satellite Dec (At Culmination)", lambda x: x["dec"]._degrees),
		("Satellite Altitude (At Culmination)", lambda x: x["alt"]._degrees),
		("Satellite Azimuth (At Culmination)", lambda x: x["az"]._degrees),
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
			rowS = ",".join([csvTrans(data) for csvTrans in csvTransforms]) + "\n"
			fp.write(rowS)


def main():
	global calipsoNs, todayUTC, nDays, solveForSatPos

	finalDay = todayUTC + timedelta(days=nDays)

	print("Calculating CALIPSO culminations")

	for calipsoCalc in calipsoNs:
		print(f"Looking for CALIPSO at {calipsoCalc.getObsName()} from {todayUTC} to {finalDay}")

		if solveForSatPos:
			satStore = h5py.File(calipsoCalc.getObsName() + '.hdf5', 'w')

		culmData = []
		for dayNum in tqdm(range(nDays)):

			nextDay = todayUTC + timedelta(days=dayNum)
			dayTransits = calipsoCalc.getDayTransits(nextDay, 30)

			for dayTransit in dayTransits:

					culmination = dayTransit.getCulmination()
					satData = calipsoCalc.getSatilliteDataAt(culmination)
					dayTransit.addToDict(satData)
					culmData.append(satData)

					if solveForSatPos:
							startT = culmination.utc_datetime() - timedelta(seconds=10)
							storePos = calipsoCalc.getObsCoordBetween(startT, 20, 1)
							#Storage step
							datasetName = startT.isoformat()
							satStore.create_dataset(datasetName, data=storePos)

		turnToCSV(calipsoCalc.getObsName() + ".csv", culmData)

		if solveForSatPos:
			satStore.close()


	if solveForSatPos:
		satStore.close()

	if len(culmData) > 0:
		print(f"Found CALIPSO culminating {len(culmData)} times at {obsName}, for more info open {obsName}_CALIPSO_Info.csv")
		print("")

	else:
		print(f"Did not find CALIPSO culminating at {obsName} from {todayUTC} to {finalDay}")

def loadSkyfieldAPI():
	global timescale, ephemeris, earth, moon
	timescale = load.timescale()
	ephemeris = load('de421.bsp')
	earth = ephemeris['earth']
	moon = ephemeris['moon']

def buildTimeGlobals():
	global startTime, todayUTC, utc
	startTime = datetime.min.time()
	todayUTC = timezone_converter(datetime.combine(datetime.today().date(), startTime))
	utc = timezone("UTC")

def buildSatGlobals():
	global reloadSat, calipso, todayUTC
	satURL = "https://celestrak.com/NORAD/elements/active.txt"
	satellites = None
	if reloadSat:
		satellites = load.tle_file(satURL, reload=True)
	else:
		satellites = load.tle_file(satURL)

	todaySFT = timescale.utc(todayUTC)
	for sat in satellites:
		days = todaySFT - sat.epoch
		if (days < 14) and ("CALIPSO" in sat.name):
			calipso = sat
			break

	if calipso is None:
		print("CALIPSO was not found, please retry with reloadSat=True or check NASA to ensure it is still active")
		sys.exit(1)


def buildScriptGlobals():
	loadSkyfieldAPI()
	buildTimeGlobals()
	buildSatGlobals()


def buildScriptInputs():
	global obsData, nDays, localTZ, solveForSatPos, reloadSat, solveForSatPos
	#Script Inputs
	# -----------------------
	# Reload Satillite Data?
	reloadSat = False

	# Location Inputs
	# "<Obsevatory Name>" : ['<input latitude>', '<input longitude>', <altitude>]
	obsData = {
		"Lick": [37.341 * N, 121.642 * W, 1282],
		"VERITAS": [31.675 * N, 110.9519 * W, 1280]
	}

	# Calculate for the next n days from today
	nDays = 120
	# timezone output
	localTZ = get_localzone()
	# localTZ = timezone('America/Phoenix')

	#If for some reason you want to specifically
	#calculate for the positions of CALIPSO around
	#culmination
	solveForSatPos = False
	# -----------------------

	if solveForSatPos:
		try:
			import h5py
		except ImportError:
			print("In order to store satillite positons please install h5py https://pypi.org/project/h5py/")
			sys.exit(1)

if __name__ == "__main__":
	buildScriptInputs()
	buildScriptGlobals()

	for obsName, obsArr in obsData.items():
		observatory = iers2010.latlon(*obsArr)
		earthObservatory = earth + observatory
		calipsoNs.append(calipsoCalculator(calipso, obsName, observatory, earthObservatory))

	main()
