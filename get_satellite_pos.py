import argparse
import os
import pathlib

from tqdm import tqdm

from common import *

try:
	import h5py
except ImportError as err:
	print("In order to store satellite positions please install h5py https://pypi.org/project/h5py/\n\n")
	raise err

localTZ = None


class HDF5FileNotSetup(IOError):
	pass


class HDF5FileHandler:
	def __init__(self):
		self.file_cache = {}
		self.cur_file = ""
		self.cur_group = None
		self.is_working = False
	
	def open_file(self, file_name):
		self.is_working = True
		self.cur_file = file_name
		self.cur_group = None
		if file_name in self.file_cache and self.file_cache[file_name]:
			return
		else:
			self.file_cache[file_name] = h5py.File(file_name, mode='a')
	
	def add_group(self, group_name):
		if not self.is_working:
			raise HDF5FileNotSetup
		if self.cur_group is None:
			self.cur_group = self.file_cache[self.cur_file].require_group(group_name)
		else:
			self.cur_group = self.cur_group.require_group(group_name)
	
	def add_dataset(self, dataset_name, data_input):
		if not self.is_working:
			raise HDF5FileNotSetup
		if self.cur_group is None:
			self.file_cache[self.cur_file].create_dataset(dataset_name, data=data_input)
		else:
			self.cur_group.create_dataset(dataset_name, data=data_input)
	
	def finish_all(self):
		for file_name in self.file_cache.keys():
			if self.file_cache[file_name]:
				self.file_cache[file_name].close()
	
	def close_working_file(self):
		if self.is_working:
			self.file_cache[self.cur_file].close()
		self.cur_file = ""
		self.cur_group = None
		self.is_working = False


def main(reloadSat_: bool, useAll_: bool, merge_: bool, nDays: int, inCoord_: Tuple[float, float], inRad_: float, faltaz_: int, ignore_empty: bool) -> None:
	today_utc = TimeObj(datetime.now())
	
	startDay = today_utc.get_start_day()
	finalDay = startDay.get_off(days=nDays)
	
	sat_calculators = build_sat_calcs(reloadSat_, useAll_)
	
	inputC1 = numpy.full((1, 86400), inCoord_[0])
	inputC2 = numpy.full((1, 86400), inCoord_[1])
	
	file_handler = HDF5FileHandler()
	
	for sat_calc in sat_calculators:
		if merge_:
			out_file = os.path.join("output", "calculated_satellite_data.hdf5")
			file_handler.open_file(out_file)
			file_handler.add_group(sat_calc.get_obs_name())
			file_handler.add_group(sat_calc.get_sat_name())
			file_handler.add_group(str(numpy.around(inRad_, 3)))
		else:
			out_file = os.path.join("output", "sat_pos", f"{sat_calc.get_sat_name()}_{sat_calc.get_obs_name()}_{int(inRad_ * 100)}_{today_utc.get_file_format()}-{finalDay.get_file_format()}.hdf5")
			file_handler.open_file(out_file)
		
		print(
			f"Looking for {sat_calc.get_sat_name()} at "
			f"{sat_calc.get_obs_name()} from {today_utc} to {finalDay} "
			f"around {inCoord_} at a radius of {inRad_}"
		)
		
		for dayNum in tqdm(range(nDays)):
			start_t = startDay.get_off(days=dayNum)
			day_pos = sat_calc.get_obs_coord_between(start_t, 86400., 1., True)
			
			sat_ind = (3, 4) if bool(faltaz_) else (1, 2)
			sat_c1 = day_pos[:, sat_ind[0]]
			sat_c2 = day_pos[:, sat_ind[1]]
			dis_arr = two_object_distance_explict(sat_c1, sat_c2, inputC1[0], inputC2[0])
			ind_arr = numpy.nonzero(dis_arr < inRad_)[0]
			good_arr = day_pos[ind_arr, :]
			
			# Storage step
			if ignore_empty and (0 in good_arr.shape):
				continue
			dataset_name = start_t.iso_format()
			file_handler.add_dataset(dataset_name, good_arr)
	
	file_handler.finish_all()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description='search the inputted satellites track for how close each will get to your input.')
	parser.add_argument('days', type=int, help='The number of days to calculate for')
	parser.add_argument('coordinate_1', type=float, help="RA or Alt depending on <flag>")
	parser.add_argument('coordinate_2', type=float, help="Dec or Az depending on <flag>")
	parser.add_argument('radius', type=float, help="The radius around the coordinate")
	parser.add_argument('flag', type=int, help="0 = RA/Dec; 1 = Alt/Az")
	parser.add_argument('-r', '--reload', action='store_true', help="Re-download the TLE file from the internet")
	parser.add_argument('-a', '--all', action='store_true', help="Use all satellites from database")
	parser.add_argument('-m', '--merge', action='store_true', help="write all data into a single file")
	parser.add_argument('-i','--ignore-empty', action='store_true', help="If the dataset is empty, it will not be written to the file")
	
	args = parser.parse_args()
	
	nDays = args.days
	reloadSat = args.reload
	
	pathlib.Path(os.path.join("output", "sat_pos")).mkdir(parents=True, exist_ok=True)
	main(reloadSat, args.all, args.merge, nDays, (args.coordinate_1, args.coordinate_2), args.radius, args.flag, args.ignore_empty)
