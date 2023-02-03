import argparse
import os
from pathlib import Path
from typing import Tuple

import numpy
from tqdm import tqdm

from utils import (
    build_obs_sat,
    get_obs_coord_between,
    two_object_distance,
    TimeDeltaObj,
    today,
)

try:
    import h5py
except ImportError as err:
    print(
        "In order to store satellite positions please install h5py https://pypi.org/project/h5py/\n\n"
    )
    raise err


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
            self.file_cache[file_name] = h5py.File(file_name, mode="a")

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


def make_coord_check(in_radius, coord_1, coord_2, use_alt_az):
	df_ind = ('alt', 'az') if bool(use_alt_az) else ('ra', 'dec')
	
	def coord_check(in_df):
		return two_object_distance(in_df[df_ind[0]], in_df[df_ind[1]], coord_1, coord_2) < in_radius
	
	return coord_check

def main(pargs) -> None:
	today_utc = today()
	
	start_day = today_utc.get_start_day()
	final_day = start_day.get_off(days=pargs.days)

	step_size = pargs.step if (pargs.step is not None) and (pargs.step > 1) else 1
	step_t = TimeDeltaObj(seconds=step_size)

	if args.search:
		coord_checker = make_coord_check(
			args.search[2], 
			args.search[1], 
			args.search[0], 
			bool(args.search[3])
		)
	else:
		coord_checker = None
	
	file_handler = HDF5FileHandler()
	tles = [None] if pargs.tles is None else pargs.tles

	for tle in tles:
		tle_name = "default" if tle is None else tle.stem
		observatory_sats = build_obs_sat(pargs.reload, pargs.sat_cache, pargs.all, tle)
		for obs_sat_obj in observatory_sats:
			if pargs.merge:
				out_file = os.path.join("output", "calculated_satellite_data.hdf5")
				file_handler.open_file(out_file)
				if tle is not None:
					file_handler.add_group(tle_name)
				file_handler.add_group(obs_sat_obj.obs_name)
				file_handler.add_group(obs_sat_obj.sat_name)
				if coord_checker is not None:
					file_handler.add_group(str(numpy.around(args.search[2], 3)))
			else:
				rad_str = "" if coord_checker is None else str(int(args.search[2] * 100))
				out_file = os.path.join(
					"output", "sat_pos",
					f"{tle_name}_{obs_sat_obj.sat_name}_{obs_sat_obj.obs_name}_"
					f"{rad_str}_{today_utc.get_file_format()}-"
					f"{final_day.get_file_format()}.hdf5"
				)
				file_handler.open_file(out_file)
			
			user_out = "" if coord_checker is None else f"around {(args.search[0], args.search[1])} at a radius of {args.search[2]} "
			print(
				f"Looking for {obs_sat_obj.sat_name} at "
				f"{obs_sat_obj.obs_name} from {today_utc} to {final_day} " + user_out +
				f"using the {tle_name} TLE file"
			)
			
			for day_num in tqdm(range(pargs.days)):
				start_t = start_day.get_off(days=day_num)
				end_t = start_t.get_off(days=1)
				pd_sat_data = get_obs_coord_between(obs_sat_obj, start_t, end_t, step_t)

				if coord_checker is not None:
					pd_sat_data = pd_sat_data[coord_checker(pd_sat_data)]
				
				# Storage step
				if pargs.ignore_empty and (0 in pd_sat_data.shape):
					continue

				dataset_name = start_t.iso_format()
				file_handler.add_dataset(dataset_name, pd_sat_data.to_numpy())
	
	file_handler.finish_all()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description='search the inputted satellites track for how close each will get to your input.')
	parser.add_argument('days', type=int, help='The number of days to calculate for')
	parser.add_argument("--search", nargs=4, type=float, help="<RA/Alt> <Dec/Az> <Radius> <Flag>")
	parser.add_argument('-r', '--reload', action='store_true', help="Re-download the TLE file from the internet")
	parser.add_argument('-sc', '--sat-cache', action='store_true', help="Automatically Cache TLE files")
	parser.add_argument('-a', '--all', action='store_true', help="Use all satellites from database")
	parser.add_argument('-m', '--merge', action='store_true', help="write all data into a single file")
	parser.add_argument('-i', '--ignore-empty', action='store_true', help="If the dataset is empty, it will not be written to the file")
	parser.add_argument('--tles', type=Path, nargs="+", help="Specify tle files to use instead of downloading them")
	parser.add_argument('--step', type=int, nargs=1, help="Step size in seconds", deault=1)
	
	args = parser.parse_args()
	
	Path(os.path.join("output", "sat_pos")).mkdir(parents=True, exist_ok=True)
	main(args)
