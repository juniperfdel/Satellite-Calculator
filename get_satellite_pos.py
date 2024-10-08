import argparse
from itertools import pairwise
from pathlib import Path
from typing import Callable, Optional

import numpy
import pandas
from numpy.typing import ArrayLike
from tqdm import tqdm

from utils.arg_utils import add_common_params
from utils.common import (
    ObservatorySatelliteFactory,
    get_obs_coord_between,
    make_bounded_time_list,
)
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


def make_coord_check(
    in_radius: float, coord_1: float, coord_2: float, use_alt_az: bool
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
        end_date_or_days=pargs.end,
        ignore_limit=pargs.ignore_limit,
        tles=pargs.tles,
        chain_tles=pargs.chain,
    )

    step_size = int(pargs.step) if (pargs.step is not None) and (pargs.step > 1) else 1
    step_t = TimeDeltaObj(seconds=step_size)
    day_step = TimeDeltaObj(days=1)

    coord_group_name = ""
    user_out = ""
    if pargs.filter_radius:
        coord_checker = make_coord_check(
            pargs.filter_radius[2],
            pargs.filter_radius[0],
            pargs.filter_radius[1],
            bool(int(pargs.filter_radius[3])),
        )
        ns = lambda x: numpy.around(pargs.filter_radius[x], 3)
        coord_group_name = f"/c1_{ns(1)}_c2_{ns(1)}_r_{ns(2)}".replace(".", "")
        user_out = f"{(pargs.filter_radius[0], pargs.filter_radius[1])};{args.filter_radius[2]};"
    else:
        coord_checker = None

    obs_sat_fact_bar = tqdm(obs_sat_fact, position=0)

    file_name = str(output_path / f"{pargs.output_file}.hdf5")
    hdf5_file = h5py.File(file_name, mode="a")

    for start_utc, days_to_calc, obs_sat_obj in obs_sat_fact_bar:
        final_day = start_utc + TimeDeltaObj(days=days_to_calc)
        start_end_pairs = list(
            pairwise(make_bounded_time_list(start_utc, final_day, day_step))
        )

        if pargs.chain:
            cur_group: h5py.Group = hdf5_file.require_group(
                f"{coord_group_name}/{obs_sat_obj.obs_name}/{obs_sat_obj.sat_name}"
            )
        else:
            cur_group: h5py.Group = hdf5_file.require_group(
                f"{coord_group_name}/{obs_sat_obj.obs_name}/{obs_sat_obj.sat_name}/{obs_sat_obj.sat_epoch_obj.iso_format()}"
            )

        obs_sat_fact_bar.set_description(
            f"{obs_sat_obj.sat_name};{obs_sat_obj.obs_name};"
            f"{start_utc.get_compact_fmt()}--{final_day.get_compact_fmt()};"
            f"{user_out}"
            f"Ep@{obs_sat_obj.compact_sat_epoch_str}"
        )

        for start_t, end_t in tqdm(start_end_pairs, position=1, leave=False):
            days_since_epoch = (start_t - obs_sat_obj.sat_epoch_obj).total_days()

            if pargs.force_tle_limit and (
                (days_since_epoch < 0) or (pargs.force_tle_limit < days_since_epoch)
            ):
                continue

            pd_sat_data = get_obs_coord_between(obs_sat_obj, start_t, end_t, step_t)

            if coord_checker is not None:
                pd_sat_data = pd_sat_data[coord_checker(pd_sat_data)]

            # Storage step
            if pargs.ignore_empty and (0 in pd_sat_data.shape):
                continue

            new_ds: h5py.Dataset = cur_group.create_dataset(
                start_t.iso_format(), data=pd_sat_data.to_numpy()
            )

            new_ds.attrs["column_names"] = (
                numpy.array(pd_sat_data.columns, dtype="S"),
            )
            new_ds.attrs["step_size_sec"] = step_size
            if pargs.chain:
                new_ds.attrs["epoch"] = obs_sat_obj.sat_epoch_obj.iso_format()

    hdf5_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate the inputted satellites position data"
    )

    add_common_params(parser, "calculated_satellite_data")

    parser.add_argument(
        "--filter-radius",
        type=float,
        nargs=4,
        help="Filter output to only include positions inside a radius around a sky position; Format: <RA/Alt> <Dec/Az> <Radius in Degrees> <0/1> where 0 = RA/Dec; 1 = Alt/Az",
    )

    parser.add_argument(
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

    args = parser.parse_args()

    output_path.mkdir(parents=True, exist_ok=True)
    main(args)
