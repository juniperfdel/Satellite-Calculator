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
from skyfield.api import Star, EarthSatellite, Topos, load, wgs84
from skyfield.units import Angle as skyAngle

from geopy import distance


#Stolen from https://stackoverflow.com/questions/27531718/datetime-timezone-conversion-using-pytz
def timezone_converter(input_dt, current_tz='UTC', target_tz='America/Phoenix'):
	current_tz = timezone(current_tz)
	target_tz = timezone(target_tz)
	target_dt = None
	if input_dt.tzinfo is not None and input_dt.tzinfo.utcoffset(input_dt) is not None:
		target_dt = input_dt.astimezone(target_tz)
	else:
		target_dt = current_tz.localize(input_dt).astimezone(target_tz)
	return target_tz.normalize(target_dt)

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
	 Tt = ts.utc(dTt)
	 return almanac.fraction_illuminated(planets,"moon",Tt)


def moonPosition(dTt, earthObserver):
	Tt = ts.utc(dTt)
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

# a helper function for finding the sun set/rise
# from https://stackoverflow.com/questions/54485777/finding-twilight-times-with-skyfield
def daylength(ephemeris, topos, degrees):
	"""Build a function of time that returns the daylength.

	The function that this returns will expect a single argument that is a
	:class:`~skyfield.timelib.Time` and will return ``True`` if the sun is up
	or twilight has started, else ``False``.
	"""
	sun = ephemeris['sun']
	topos_at = (ephemeris['earth'] + topos).at

	def is_sun_up_at(t):
		"""Return `True` if the sun has risen by time `t`."""
		t._nutation_angles = iau2000b(t.tt)
		return topos_at(t).observe(sun).apparent().altaz()[0].degrees > -degrees

	is_sun_up_at.rough_period = 0.5  # twice a day
	return is_sun_up_at

def sunRiseSetFinder(startDT, numDays, obsLocation):
	global startTime
	sunSetDict = {}
	for nDay in tqdm(range(numDays)):
		currentDT = startDT + timedelta(days=nDay)
		currentDayKey = datetime.strftime(currentDT, "%Y-%m-%d")

		currentDay = datetime.combine(currentDT.date(), startTime)
		currentDay = timezone_converter(currentDay,target_tz="UTC")
		t0 = ts.utc(currentDay)
		t1 = ts.utc(currentDay + timedelta(days=1))

		astro_twil_time, astro_twil_up = almanac.find_discrete(t0, t1, daylength(planets, obsLocation, DAYLENGTH_ASTRONOMICAL_TWILIGHT))

		sunSetDict[currentDayKey] = {}
		tArr = astro_twil_time.utc_datetime()
		for i,v in enumerate(astro_twil_up):
			key = "set"
			if v:
				key = "rise"
			sunSetDict[currentDayKey][key] = tArr[i]
	return sunSetDict

# Uses the data above to figure out if a given datetime is in the day or not
def isInDay(currentDT, sunTimeDict):
	curDayK = datetime.strftime(currentDT, "%Y-%m-%d")
	curSunD = sunTimeDict[curDayK]
	riseDT = curSunD["rise"]
	setDT = curSunD["set"]
	if ((setDT < currentDT) and (currentDT < riseDT)):
		return False
	return True


# Helper class used in all scripts calculating CALIPSO information
# to reduce the amount of boiler-plate code
class calipsoNode:
	def __init__(self, satt, observatory):
		self.sat = satt
		self.obs = observatory
		self.diff = self.sat - self.obs

	def getEpoch(self):
		return self.sat.epoch

	def getObsCoord(self, inTime):
		toc = self.diff.at(inTime)
		alt, az, _ = toc.altaz()
		ra, dec, _ = toc.radec()
		return (alt, az, ra, dec)

	def getEvents(self, inDate, input_alt):
		startDay = datetime.combine(datetime(inDate[0],inDate[1],inDate[2]), startTime)
		startDay = timezone_converter(startDay, current_tz='UTC', target_tz='UTC')

		t0 = ts.utc(startDay)
		t1 = ts.utc(startDay + timedelta(days=1))
		return self.sat.find_events(self.obs, t0, t1, altitude_degrees=input_alt)

def nNones(lst):
	return sum(x is not None for x in lst)

#Script Inputs
# -----------------------

# Location Inputs
# "<Obsevatory Name>" : ['<input latitude>', '<input longitude>', <altitude>]
obsData = {
	"Lick": ['37.341 N', '121.642 W', 1282],
	"VERITAS": ['31.675 N', '110.9519 W', 1280]
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

#Useful Globals
startTime = datetime.min.time()
ts = load.timescale()
planets = load('de421.bsp')
utc = timezone("UTC")

satURL = "https://celestrak.com/NORAD/elements/active.txt"
satellites = load.tle_file(satURL)#, reload=True)

earth = planets['earth']
moon = planets['moon']

todayUTC = timezone_converter(datetime.combine(datetime.today().date(), startTime), target_tz="UTC")

for obsName, obsInfo in obsData.items():
	finalDay = todayUTC + timedelta(days=nDays)
	print(f"Looking for CALIPSO at {obsName} from {todayUTC} to {finalDay}")
	lat = obsInfo[0]
	lon = obsInfo[1]
	elv = obsInfo[2]

	obsLocation = Topos(lat, lon, elevation_m=elv)
	earthObserver = earth + obsLocation

	print(f"Solving for Sun rises and sets at {obsName}")
	sunRiseSetDict = sunRiseSetFinder(todayUTC, nDays, obsLocation)

	now = localTZ.localize(datetime.now())
	t = ts.from_datetime(now)
	activeSat = []
	for sat in satellites:
		days = t - sat.epoch
		if days < 14:
			activeSat.append(sat)

	calSat = None
	for sat in activeSat:
		if ("CALIPSO" in sat.name):
			print("Found CALIPSO in active database!")
			calSat = sat
	if calSat is None:
		print("CALIPSO was not found, please retry with reload=True or check NASA to ensure it is still active")
		sys.exit(1)

	culmData = []
	cal = calipsoNode(calSat, obsLocation)

	print("Calculating CALIPSO culminations")
	for dayNum in tqdm(range(nDays)):
		nextDay = todayUTC + timedelta(days=dayNum)
		ti, eev = cal.getEvents([nextDay.year, nextDay.month, nextDay.day], 30)
		if (len(eev) >= 3):
				rF = {}
				dd = [None, None, None]
				for ti, event in zip(ti, eev):
					name = ('rise', 'culminate', 'fall')[event]
					rF[name] = ti
					dd[event] = event
					numNotNone = nNones(dd)
					if numNotNone == 3:
						satData = {}
						satData = rF
						satData['name'] = calSat.name
						minA = (rF['fall'] - rF['rise']) * 24 * 60
						satData['dur'] = minA

						geocentric = sat.at(rF['culminate'])
						subpoint = wgs84.subpoint(geocentric)

						satLoc = (subpoint.latitude.degrees, subpoint.longitude.degrees)
						verLoc = (obsLocation.latitude.degrees, obsLocation.longitude.degrees)
						greatCDist = distance.distance(satLoc,verLoc)
						satData['lat'] = subpoint.latitude
						satData['long'] = subpoint.longitude

						satData['approx_dis'] = greatCDist.km

						sunlit = geocentric.is_sunlit(planets)
						difference = calSat - obsLocation
						topocentric = difference.at(rF['culminate'])

						ra, dec, _ = topocentric.radec()  # ICRF ("J2000")
						alt, az, _ = topocentric.altaz()

						satData['ra'] = ra
						satData['dec'] = dec

						satData['alt'] = alt
						satData['az'] = az

						satData['sun'] = sunlit
						satData['isDay'] = isInDay(rF['culminate'].utc_datetime(),sunRiseSetDict)
						moonDis, moonFrac = disFromMoon(rF['culminate'].utc_datetime(), alt, az, earthObserver)
						satData['moonDis'] = moonDis
						satData['moonFrac'] = moonFrac

						if solveForSatPos:
							startT = rF['culminate'].utc_datetime() - timedelta(seconds=10)
							secA = 20
							secI = numpy.ceil(secA).astype('int')
							posArr = []
							for tStep in range(0,secI,1):
								nearFut = startT + timedelta(seconds=tStep)
								tNearFut = ts.from_datetime(nearFut)
								curSatLoc = difference.at(tNearFut)
								curSatRA, curSatDec, _ = curSatLoc.radec()
								posArr.append([nearFut, curSatRA, curSatDec])
							satData['satTrack'] = list(posArr)

						culmData.append(satData)
						rF = {}
						dd = [None, None, None]

	if len(culmData) > 0:
		print(f"Found CALIPSO culminating {len(culmData)} times at {obsName}, for more info open {obsName}_CALIPSO_Info.csv")
		print("")

		tzNameAb = culmData[0]['culminate'].astimezone(localTZ).strftime('%Z')
		head = ('Observatory Latitude', 'Observatory Longitude', 'Satellite name', 'Culmination Time (UTC)' ,f'Culmination Time ({tzNameAb})', 'Daytime at Culmination?', 'Satellite Shadow Latitude (At Culmination)' , 'Satellite Shadow Longitude (At Culmination)','Approx. Distance Between Shadow and Observatory (km)' , 'Satellite RA (At Culmination)', 'Satellite Dec (At Culmination)', 'Satellite Altitude (At Culmination)', 'Satellite Azimuth (At Culmination)', 'Duration Above 30 deg alt [minutes]', 'Is the Satellite Sunlit?', 'Moon Distance (At Culmination)', 'Moon Fraction (At Culmination)')
		with open(f"{obsName}_CALIPSO_Info.csv","w") as satFile:
			satFile.write(f"{','.join(head)}\n")
			for satD in culmData:
				satFile.write(f"{obsLocation.latitude}, {obsLocation.longitude}, {satD['name']}, {satD['culminate'].astimezone(utc).strftime('%Y %b %d %H:%M:%S')} , {satD['culminate'].astimezone(localTZ).strftime('%Y %b %d %H:%M:%S')}, {'Y' if satD['isDay'] else 'N'}, {satD['lat']}, {satD['long']}, {around(satD['approx_dis'])}, {satD['ra']._degrees}, {satD['dec']._degrees}, {satD['alt']}, {satD['az']}, {around(satD['dur'],2)}, {'Y' if satD['sun'] else 'N'}, {around(satD['moonDis'],2)},  {around(satD['moonFrac'],3)} \n")
	else:
		print(f"Did not find CALIPSO culminating at {obsName} from {todayUTC} to {finalDay}")
