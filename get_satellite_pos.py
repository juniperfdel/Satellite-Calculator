import argparse
from itertools import pairwise
from pathlib import Path
from typing import Callable, Optional

import numpy
import pandas
from numpy.typing import ArrayLike
from tqdm import tqdm

from utils.arg_utils import add_common_params
from utils.common import (ObservatorySatelliteFactory, get_obs_coord_between,
                          make_bounded_time_list)
from utils.math_utils import two_object_distance
from utils.time_utils import TimeDeltaObj

try:
    import h5py
except ImportError as err:
    print(
        "In order to store satellite positions please install h5py https://pypi.org/project/h5py/\n\n"
    )
    raise err

output_path = Path("output") / "sat_pos"


class HDF5FileNotSetup(IOError):
    pass


class HDF5FileHandler:
    def __init__(self):
        self.file_cache: dict[str, h5py.File] = {}
        self.cur_file: str = ""
        self.cur_group: Optional[h5py.Group] = None
        self.is_working: bool = False

    def open_file(self, file_name: str):
        self.is_working = True
        self.cur_file = file_name
        self.cur_group = None
        if file_name in self.file_cache and self.file_cache[file_name]:
            return
        else:
            fpath = Path(file_name)
            if fpath.exists():
                fpath.unlink()
            self.file_cache[file_name] = h5py.File(file_name, mode="a")

    def add_group(self, group_name: str):
        if not self.is_working:
            raise HDF5FileNotSetup
        if self.cur_group is None:
            self.cur_group = self.file_cache[self.cur_file].require_group(group_name)
        else:
            self.cur_group = self.cur_group.require_group(group_name)

    def add_dataset(self, dataset_name: str, data_input: ArrayLike, **iattrs):
        if not self.is_working:
            raise HDF5FileNotSetup
        if self.cur_group is None:
            new_ds: h5py.Dataset = self.file_cache[self.cur_file].create_dataset(
                dataset_name, data=data_input
            )
        else:
            new_ds: h5py.Dataset = self.cur_group.create_dataset(
                dataset_name, data=data_input
            )
        for k, v in iattrs.items():
            new_ds.attrs[k] = v

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


def make_coord_check(
    in_radius, coord_1, coord_2, use_alt_az
) -> Callable[[pandas.DataFrame], pandas.DataFrame]:
    df_ind = ("alt", "az") if bool(use_alt_az) else ("ra", "dec")

    def coord_check(in_df: pandas.DataFrame) -> pandas.DataFrame:
        return (
            two_object_distance(in_df[df_ind[0]], in_df[df_ind[1]], coord_1, coord_2)
            < in_radius
        )

    return coord_check


def main(pargs: argparse.Namespace) -> None:
    obs_sat_fact = ObservatorySatelliteFactory(
        local_tz=pargs.timezone,
        reload_sat=pargs.reload,
        cache_sat=pargs.cache,
        use_all=pargs.all,
        start_date=pargs.start_date,
        end_days=pargs.days,
        ignore_limit=pargs.ignore_limit,
        tles=pargs.tles,
        chain_tles=pargs.chain,
    )

    step_size = int(pargs.step) if (pargs.step is not None) and (pargs.step > 1) else 1
    step_t = TimeDeltaObj(seconds=step_size)
    day_step = TimeDeltaObj(days=1)

    if pargs.flag < 2:
        coord_checker = make_coord_check(
            pargs.radius, pargs.coord1, pargs.coord2, bool(pargs.flag)
        )
    else:
        coord_checker = None

    file_handler = HDF5FileHandler()

    obs_sat_fact_bar = tqdm(obs_sat_fact, position=0)

    for start_utc, days_to_calc, obs_sat_obj in obs_sat_fact_bar:
        final_day = start_utc + TimeDeltaObj(days=days_to_calc)
        start_end_pairs = list(
            pairwise(make_bounded_time_list(start_utc, final_day, day_step))
        )

        out_file = str(output_path / f"{pargs.output_file}.hdf5")
        file_handler.open_file(out_file)
        file_handler.add_group(obs_sat_obj.sat_epoch_obj.iso_format())
        file_handler.add_group(obs_sat_obj.obs_name)
        file_handler.add_group(obs_sat_obj.sat_name)
        if coord_checker is not None:
            file_handler.add_group(str(numpy.around(pargs.search[2], 3)))

        user_out = (
            ""
            if coord_checker is None
            else f"{(pargs.search[0], pargs.search[1])};{args.search[2]};"
        )
        obs_sat_fact_bar.set_description(
            f"{obs_sat_obj.sat_name};{obs_sat_obj.obs_name};"
            f"{start_utc.get_compact_fmt()}--{final_day.get_compact_fmt()};"
            + user_out
            + f"Ep@{obs_sat_obj.compact_sat_epoch_str}"
        )

        for start_t, end_t in tqdm(start_end_pairs, position=1, leave=False):
            days_since_epoch = (start_t - obs_sat_obj.sat_epoch_obj).total_days()

            if (pargs.force_tle_limit is not None) and (
                (days_since_epoch < 0) or (pargs.force_tle_limit < days_since_epoch)
            ):
                continue

            pd_sat_data = get_obs_coord_between(obs_sat_obj, start_t, end_t, step_t)

            if coord_checker is not None:
                pd_sat_data = pd_sat_data[coord_checker(pd_sat_data)]

            # Storage step
            if pargs.ignore_empty and (0 in pd_sat_data.shape):
                continue

            dataset_name = start_t.iso_format()
            file_handler.add_dataset(
                dataset_name,
                pd_sat_data.to_numpy(),
                column_names=numpy.array(pd_sat_data.columns, dtype="S"),
                step_size_sec=step_size,
            )

    file_handler.finish_all()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="search the inputted satellites track for how close each will get to your input."
    )

    parser.add_argument("coord1", type=float, help="RA or Alt depending on <flag>")
    parser.add_argument("coord2", type=float, help="Dec or Az depending on <flag>")
    parser.add_argument("radius", type=float, help="The radius around the coordinate")
    parser.add_argument("flag", type=int, help="0 = RA/Dec; 1 = Alt/Az; 2 = ignore")

    add_common_params(parser, "calculated_satellite_data")

    parser.add_argument(
        "-i",
        "--ignore-empty",
        action="store_true",
        help="If the dataset is empty, it will not be written to the file",
    )

    parser.add_argument("--step", type=int, default=1, help="Step size in seconds")

    parser.add_argument(
        "--force-tle-limit",
        type=float,
        help="Don't calculate past this amount of days from the model epoch",
    )

    parser.add_argument(
        "--chain",
        action="store_true",
        help="In the specific case of one satellite having many sequential models - chain them together",
    )

    args = parser.parse_args()

    output_path.mkdir(parents=True, exist_ok=True)
    main(args)
