from common import *

try:
	import h5py
except ImportError as err:
	print("In order to store satellite positons please install h5py https://pypi.org/project/h5py/\n\n")
	raise err

localTZ = None

def main(reloadSat_: bool, useAll_: bool, nDays: int, inCoord_: Tuple[float, float], inRad_: float, faltaz_: int) -> None:
	global todayUTC

	startDay = todayUTC.getStartDay()
	finalDay = startDay.getOff(days=nDays)

	satCalculators = buildSatCalcs(reloadSat_, useAll_)

	inputC1 = numpy.full((1, 86400), inCoord_[0])
	inputC2 = numpy.full((1, 86400), inCoord_[1])
	for satCalc in  satCalculators:
		outFileName = os.path.join("output","sat_pos",f"{satCalc.getSatName()}_{satCalc.getObsName()}_{int(inRad_*100)}_{startDay.getfileformat()}-{finalDay.getfileformat()}.hdf5")
		satStore = h5py.File(outFileName, 'w')
		print(f"Looking for {satCalc.getSatName()} at {satCalc.getObsName()} from {startDay} to {finalDay} around {inCoord_} at a radius of {inRad_}")
		for dayNum in tqdm(range(nDays)):
			startT = startDay.getOff(days=dayNum)
			dayPos = satCalc.getObsCoordBetween(startT, 86400., 1., True)
			
			if bool(faltaz_):
				satInd = (3,4)
			else:
				satInd = (1,2)
			satC1 = dayPos[:,satInd[0]]
			satC2 = dayPos[:, satInd[1]]
			disArr = twoObjectDistanceExplict(satC1, satC2, inputC1[0], inputC2[0])
			indArr = numpy.nonzero(disArr < inRad_)[0]
			goodArr = dayPos[indArr, :]

			#Storage step
			datasetName = startT.isoformat()
			satStore.create_dataset(datasetName, data=goodArr)
		satStore.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='search the inputted satellites track for how close each will get to your input.')
	parser.add_argument('days', type=int, help='The number of days to calculate for')
	parser.add_argument('coordinate_1', type=float, help="RA or Alt depending on <flag>")
	parser.add_argument('coordinate_2', type=float, help="Dec or Az depending on <flag>")
	parser.add_argument('radius', type=float, help="The radius around the coordinate")
	parser.add_argument('flag', type=int, help="0 = RA/Dec; 1 = Alt/Az")
	parser.add_argument('-tz','--set-timezone', type=str, help="The timezone to cacluate with respect to, default is your local timezone")
	parser.add_argument('-r','--reload', action='store_true', help="Re-download the TLE file from the internet, default is false")
	parser.add_argument('-a','--all', action='store_true', help="Use all satellites from database, default is false")
	args = parser.parse_args()

	nDays = args.days
	reloadSat = args.reload

	inputTZ = getattr(args,"set_timezone", None)
	if inputTZ is None:
		localTZ = get_localzone()
	else:
		localTZ = timezone(inputTZ)

	pathlib.Path(os.path.join("output","sat_pos")).mkdir(parents=True, exist_ok=True)
	main(reloadSat, args.all, nDays, (args.coordinate_1, args.coordinate_2), args.radius, args.flag)
