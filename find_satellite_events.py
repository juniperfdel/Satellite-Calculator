import argparse
import configparser
import os
import pathlib

from numpy import around
from tqdm import tqdm
from tzlocal import get_localzone

from common import *

# Input Globals
reloadSat = False
obsData = None
nDays = None
localTZ = None

# Time Globals
utc = timezone("UTC")


# from https://codereview.stackexchange.com/questions/156144/get-value-from-dictionary-given-a-list-of-nested-keys
def nested_get(input_dict: Dict[str, Any], nested_key: Sequence[str]):
	internal_dict_value = input_dict
	for k in nested_key:
		if k is None:
			break
		internal_dict_value = internal_dict_value.get(k)
		if internal_dict_value is None:
			return None
	return internal_dict_value


class ColumnDefinition:
	def __init__(self, column_key, data_location, column_parser, ini_loc=-1):
		self.column_key = column_key
		self.data_location = data_location
		self.column_parser = column_parser
		self.location = ini_loc
	
	def set_location(self, loc):
		self.location = loc
	
	def get_ini_key(self):
		return self.column_key.lower()
	
	def get_parsed_value_from_row(self, in_data_row):
		return self.column_parser(nested_get(in_data_row, self.data_location))
	
	def copy(self):
		return ColumnDefinition(self.column_key, self.data_location, self.column_parser, self.location)


def is_okay(test_: str, okay_rows: Sequence[Tuple[str, bool]]) -> bool:
	"""
	tests if the row is allowed by the csv_config file
	
	Returns
	-------
	bool
	"""
	for testS, testB in okay_rows:
		if test_.lower() in testS:
			return testB
	return False


def noop(x: Any):
	return x


def turn_to_csv(in_file_name: str, in_config_section: str, in_data: List[Dict[str, Any]], ignore_day: bool) -> None:
	global localTZ
	
	in_file_name = os.path.join("output", in_file_name)
	# Build the .csv column definitions; put the key
	# definitions outside of the lambda to prevent
	# unexpected behavior
	csv_col_defs = [
		("Observatory Latitude", ("meta", "obs_lat"), noop),
		("Observatory Longitude", ("meta", "obs_long"), noop),
		("Satellite name", ("meta", "sat_name"), noop),
		("Rise/Set Altitude", ("meta", "rise_set_alt"), noop),
		("Duration between Rise and Set [minutes]", ("meta", "dur"), noop),
		("Local Timezone", (None, None), lambda x: datetime.now().astimezone(localTZ).strftime('%Z'))
	]
	
	evTN2N = {"culm": "Culmination", "rise": "Rise", "set": "Set"}
	for evTimeName in ["culm", "rise", "set"]:
		aName = evTN2N[evTimeName]
		csv_col_defs.extend([
			(f"{aName} Time (UTC)", ("meta", evTimeName), lambda x: x.astimezone(utc).strftime('%Y %b %d %H:%M:%S')),
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
			(f"Moon Distance (At {aName})", (evTimeName, 'moon_dis'), lambda x: around(x, 2)),
			(f"Moon Fraction (At {aName})", (evTimeName, 'moon_frac'), lambda x: around(x, 3))
		])
	
	csv_col_defs = [ColumnDefinition(k, dl, fn) for k, dl, fn in csv_col_defs]
	# Filter the columns using the .ini file
	csv_config_parser = configparser.ConfigParser(allow_no_value=True)
	csv_config_parser.read_file(open(os.path.join('config', 'csv_config.ini')))
	
	good_col_defs = []
	for option_ind, option in enumerate(csv_config_parser.options(in_config_section)):
		for csv_col_def in csv_col_defs:
			if csv_col_def.get_ini_key() == option:
				csv_col_def.set_location(option_ind)
				good_col_defs.append(csv_col_def.copy())
	
	# Sort based on config file
	good_col_defs = sorted(good_col_defs, key=lambda x: x.location)
	
	# Create the .csv file
	csv_head = ",".join(x.column_key for x in good_col_defs) + "\n"
	with open(in_file_name, "w") as fp:
		fp.write(csv_head)
		for data in in_data:
			possible_row = [
				(csv_datum.column_key, str(csv_datum.get_parsed_value_from_row(data)))
				for csv_datum in good_col_defs
			]
			if ignore_day:
				
				ignore_cell = any("sunlit" in possible_cell[0].lower() and possible_cell[1].lower() == 'y' for possible_cell in possible_row)
				if ignore_cell:
					continue
			
			row_s = (",".join(cell[1] for cell in possible_row) + "\n")
			fp.write(row_s)


def main(
		reloadSat_: bool, useAll_: bool, inAltAmt: float,
		n_days: int, inConfigSection: str, ignore_day: bool
) -> None:
	today_utc = TimeObj(datetime.now())
	
	sat_calculators = build_sat_calcs(reloadSat_, useAll_)
	finalDay = today_utc.get_off(days=n_days)
	
	print(f"Calculating culminations of {len(sat_calculators)} satellites")
	
	for sat_calc in sat_calculators:
		print(f"Looking for {sat_calc.get_sat_name()} at {sat_calc.get_obs_name()} from {today_utc} to {finalDay}")
		
		culm_data = []
		for day_num in tqdm(range(n_days)):
			
			next_day = today_utc.get_off(days=day_num)
			day_transits = sat_calc.get_day_transits(next_day, inAltAmt)
			
			for day_transit in day_transits:
				sat_data = {"meta": {}}
				day_transit.add_to_dict(sat_data["meta"])
				for ev_time_name in ["culm", "rise", "set"]:
					ev_time = getattr(day_transit, ev_time_name)
					ev_data = sat_calc.get_satillite_data_at(ev_time)
					sat_data[ev_time_name] = ev_data
				culm_data.append(sat_data)
		
		if culm_data:
			print(
				f"Found {sat_calc.get_sat_name()} culminating {len(culm_data)} "
				f"times at {sat_calc.get_obs_name()}, for more info open "
				f"'output/{sat_calc.get_sat_name()}_{sat_calc.get_obs_name()}_"
				f"{today_utc.get_file_format()}-{finalDay.get_file_format()}.csv'"
			)
			print("")
			turn_to_csv(
				f"{sat_calc.get_sat_name()}_{sat_calc.get_obs_name()}_{today_utc.get_file_format()}-{finalDay.get_file_format()}.csv",
				inConfigSection, culm_data, ignore_day)
		else:
			print(f"Did not find {sat_calc.get_sat_name()} culminating at {sat_calc.get_obs_name()} from {today_utc} to {finalDay}")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description='Calculate Satellite rises, sets, and culminations for observatorys found in config/obs_data.yaml')
	parser.add_argument('days', type=int, nargs='?', default=120,
		help='The number of days to calculate for, default of 120 days')
	parser.add_argument('-tz', '--set-timezone', type=str,
		help="The timezone to calculate with respect to, default is your local timezone")
	parser.add_argument('-r', '--reload', action='store_true', help="Re-download the TLE file from the internet")
	parser.add_argument('-a', '--alt', nargs=1, type=float, default=30.,
		help="Change the altitude at which the calculations for rising and falling occur")
	parser.add_argument('-c', '--csv', type=str, default="ALL",
		help="Specifies the section in csv_config.ini to use when building the csv files")
	parser.add_argument('-i', '--ignore-daytime', action='store_true',
		help='ignore the data points which occur during the day')
	parser.add_argument('-l', '--all', action='store_true', help="Use all satellites from database")
	
	args = parser.parse_args()
	
	nDays = args.days
	reloadSat = args.reload
	
	inputTZ = getattr(args, "set_timezone", None)
	if inputTZ is None:
		localTZ = get_localzone()
	else:
		localTZ = timezone(inputTZ)
	
	solveForSatPos = getattr(args, "solve", None)
	
	if solveForSatPos is not None:
		try:
			import h5py
		except ImportError as err:
			print("In order to store satellite positions please install h5py https://pypi.org/project/h5py/\n\n")
			raise err
	
	pathlib.Path(os.path.join("output")).mkdir(parents=True, exist_ok=True)
	main(reloadSat, getattr(args, "all", False), getattr(args, "alt", 30.), nDays,
		getattr(args, "csv", "ALL"), getattr(args, "ignore_daytime", False))
