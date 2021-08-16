from common import *

# Input Globals
reloadSat = False
obsData = None
nDays = None
localTZ = None

satCalculators = []

# from https://codereview.stackexchange.com/questions/156144/get-value-from-dictionary-given-a-list-of-nested-keys
def nested_get(input_dict: Dict[str, Any], nested_key: Sequence[str]):
	internal_dict_value = input_dict
	for k in nested_key:
		if k is None:
			break
		internal_dict_value = internal_dict_value.get(k, None)
		if internal_dict_value is None:
			return None
	return internal_dict_value

def isOkay(test_: str, okayRows_: Sequence[Tuple[str, bool]]) -> bool:
	for testS, testB in okayRows_:
		if test_.lower() in testS:
			return testB

def turnToCSV(inFileName: str, inConfigSection_: str, inData:  List[Dict[str, Any]]) -> None:
	global localTZ

	inFileName = os.path.join("output",inFileName)
	noop = lambda x: x
	# Build the .csv column definitions; put the key
	# definitions outside of the lambda to prevent
	# unexpected behavior
	csvColDefs = [
		("Observatory Latitude", ("meta", "obs_lat"), noop),
		("Observatory Longitude", ("meta", "obs_long"), noop),
		("Satellite name", ("meta", "sat_name"), noop),
		("Rise/Set Altitude", ("meta", "rise_set_alt"), noop),
		("Duration between Rise and Set [minutes]", ("meta", "dur"), noop),
		("Local Timezone", (None, None), lambda x: datetime.now().astimezone(localTZ).strftime('%Z'))
	]

	evTN2N = {"culm": "Culmination", "rise": "Rise", "set":"Set"}
	for evTimeName in ["culm", "rise", "set"]:
		aName = evTN2N[evTimeName]
		csvColDefs.extend([
			(f"{aName} Time (UTC)", ("meta", evTimeName) , lambda x: x.astimezone(utc).strftime('%Y %b %d %H:%M:%S')),
			(f"{aName} Time (Local)", ("meta", evTimeName), lambda x: x.astimezone(localTZ).strftime('%Y %b %d %H:%M:%S')),
			(f"Daytime at {aName}?", (evTimeName, 'is_day'), lambda x: 'Y' if x else 'N'),
			(f"Satellite Shadow Latitude (At {aName})", (evTimeName, 'lat'), noop),
			(f"Satellite Shadow Longitude (At {aName})", (evTimeName, 'long'), noop),
			(f"Approx. Distance Between Shadow and Observatory (km) (At {aName})", (evTimeName, 'approx_dis'), noop),
			(f"Satellite RA (At {aName})", (evTimeName, 'ra'), noop),
			(f"Satellite Dec (At {aName})", (evTimeName, 'dec'), noop),
			(f"Satellite Altitude (At {aName})", (evTimeName, 'alt'), noop),
			(f"Satellite Azimuth (At {aName})", (evTimeName, 'az'), noop),
			(f"Is the Satellite Sunlit? (At {aName})", (evTimeName, 'sun'), lambda x: 'Y' if x else 'N'),
			(f"Moon Distance (At {aName})", (evTimeName, 'moon_dis'), lambda x: around(x,2)),
			(f"Moon Fraction (At {aName})", (evTimeName, 'moon_frac'), lambda x: around(x,3))
		])
	#Filter the columns using the .ini file
	csvConfigParser = configparser.ConfigParser()
	csvConfigParser.read_file(open('config/csv_config.ini'))
	configSec = csvConfigParser[inConfigSection_]
	okayCols = [(k, configSec.getboolean(k)) for k in configSec.keys()]
	csvColDefs = [x for x in csvColDefs if isOkay(x[0], okayCols)]

	# Sort based on config file
	cSortKey = list(csvConfigParser[inConfigSection_].keys())
	csvColDefs = sorted(csvColDefs, key=lambda x: cSortKey.index(x[0].lower()))

	#Create the .csv file
	csvHead = ",".join([x[0] for x in csvColDefs]) + "\n"
	csvTransforms = [(x[1], x[2]) for x in csvColDefs]
	with open(inFileName,"w") as fp:
		fp.write(csvHead)
		for data in inData:
			#Meta-programming to the rescue!
			rowS = ",".join([f"{csvT[1](nested_get(data, csvT[0]))}" for csvT in csvTransforms]) + "\n"
			fp.write(rowS)


def main(reloadSat_: bool, useAll_:bool, solveForSatPos: Union[Sequence[float], None], inAltAmt: float,  nDays: int, inConfigSection: str) -> None:
	global todayUTC

	satCalculators = buildSatCalcs(reloadSat_, useAll_)
	finalDay = todayUTC.getOff(days=nDays)

	print(f"Calculating culminations of {len(satCalculators)} satellites")

	for satCalc in satCalculators:
		print(f"Looking for {satCalc.getSatName()} at {satCalc.getObsName()} from {todayUTC} to {finalDay}")

		if solveForSatPos is not None:
			pathlib.Path(os.path.join("output","sat_pos")).mkdir(parents=True, exist_ok=True)
			satFileName = os.path.join("output","sat_pos",f"{satCalc.getSatName()}_{satCalc.getObsName()}_{todayUTC.getfileformat()}-{finalDay.getfileformat()}.hdf5")
			satStore = h5py.File(satFileName, 'w')

		culmData = []
		for dayNum in tqdm(range(nDays)):

			nextDay = todayUTC.getOff(days=dayNum)
			dayTransits = satCalc.getDayTransits(nextDay, inAltAmt)

			for dayTransit in dayTransits:
					satData = {}
					satData["meta"] = {}
					dayTransit.addToDict(satData["meta"])
					for evTimeName in ["culm", "rise", "set"]:
						evTime = getattr(dayTransit, evTimeName)
						evData = satCalc.getSatilliteDataAt(evTime)
						satData[evTimeName] = evData
					culmData.append(satData)
					if solveForSatPos is not None:
							microsecondsBefore = int(solveForSatPos[0] * -1000000)
							secondsDur = solveForSatPos[1]
							secondsStep = solveForSatPos[2]

							startT =  dayTransit.getCulmination().getOff(microseconds=microsecondsBefore)
							storePos = satCalc.getObsCoordBetween(startT, secondsDur, secondsStep, bool(int(solveForSatPos[3])))

							#Storage step
							datasetName = startT.isoformat()
							satStore.create_dataset(datasetName, data=storePos)

		if len(culmData) > 0:
			print(f"Found {satCalc.getSatName()} culminating {len(culmData)} times at {satCalc.getObsName()}, for more info open 'output/{satCalc.getSatName()}_{satCalc.getObsName()}_{todayUTC.getfileformat()}-{finalDay.getfileformat()}.csv'")
			print("")
			turnToCSV(f"output/{satCalc.getSatName()}_{satCalc.getObsName()}_{todayUTC.getfileformat()}-{finalDay.getfileformat()}.csv", inConfigSection, culmData)
		else:
			print(f"Did not find {satCalc.getSatName()} culminating at {obsName} from {todayUTC} to {finalDay}")

		if solveForSatPos:
			satStore.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Calculate Satellite rises, sets, and culminations for observatorys found in config/obs_data.yaml')
	parser.add_argument('days', type=int, nargs='?', default=120, help='The number of days to calculate for, default of 120 days')
	parser.add_argument('-tz','--set-timezone', type=str, help="The timezone to cacluate with respect to, default is your local timezone")
	parser.add_argument('-r','--reload', action='store_true', help="Re-download the TLE file from the internet, default is false")
	parser.add_argument('-s','--solve', nargs=4, type=float, help="Solve for the satellite positions around culmination, Format: <start time in seconds before culmination> <duration in seconds> <timestep in seconds> <alt/az flag>. If <alt/az flag> is set to 1, then alt/az data will be added. Setting any other number will prevent this.")
	parser.add_argument('-a', '--alt', nargs=1, type=float, default=30., help="Change the altitude at which the calculations for rising and falling occur")
	parser.add_argument('-c', '--csv', type=str, default="DEFAULT", help="Specifies the section in csv_config.ini to use when building the csv files")
	parser.add_argument('-a','--all', action='store_true', help="Use all satellites from database, default is false")
	
	args = parser.parse_args()


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
			print("In order to store satellite positons please install h5py https://pypi.org/project/h5py/\n\n")
			raise err

	pathlib.Path(os.path.join("output")).mkdir(parents=True, exist_ok=True)
	main(reloadSat, getattr(args,"all",False),  solveForSatPos, getattr(args,"alt", 30.), nDays, getattr(args,"csv", "DEFAULT"))
