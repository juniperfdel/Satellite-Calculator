import argparse
import configparser
from itertools import pairwise
import json
import os


from operator import attrgetter
from pathlib import Path
from typing import List

from tqdm import tqdm

from utils.arg_utils import add_common_params
from utils.common import ObservatorySatelliteFactory, get_day_transits, make_bounded_time_list
from utils.time_utils import TimeDeltaObj
from utils.transit import DayTransits, SingleCulmination


class ColumnDefinition:
    """Creates the column name and value getter from a definition"""

    def __init__(self, *column_definition):
        self.column_name = column_definition[0]
        self.attr_getter = attrgetter(".".join(column_definition[1:]))
        self.index = -1

    def get_ini_key(self) -> str:
        return self.column_name.lower()

    def get_value(self, in_transit_data: SingleCulmination):
        return self.attr_getter(in_transit_data)


def write_to_csv_file(
    in_file_name: str,
    in_config_section: str,
    in_days: List[DayTransits],
    ignore_day: bool,
    ignore_day_rising: bool,
    ignore_day_set: bool,
) -> None:
    in_file_name = os.path.join("output", in_file_name)

    with open("col_defs.json", "r") as fp:
        csv_col_defs = json.load(fp)

    csv_col_defs = [ColumnDefinition(*arg_defs) for arg_defs in csv_col_defs]
    # Filter the columns using the .ini file
    csv_config_parser = configparser.ConfigParser(allow_no_value=True)
    csv_config_parser.read_file(open(os.path.join("config", "csv_config.ini")))

    good_col_defs: list[ColumnDefinition] = []
    for option_ind, option in enumerate(csv_config_parser.options(in_config_section)):
        for csv_col_def in csv_col_defs:
            if csv_col_def.get_ini_key() == option:
                csv_col_def.index = option_ind
                good_col_defs.append(csv_col_def)

    # Sort based on config file
    good_col_defs.sort(key=lambda x: x.index)

    # Create the .csv file
    csv_head = ",".join(x.column_name for x in good_col_defs) + "\n"
    with open(in_file_name, "w") as fp:
        fp.write(csv_head)
        for day_transits in in_days:
            for single_culmination in day_transits:
                if single_culmination.culmination.is_daytime and ignore_day:
                    continue

                if single_culmination.start.is_daytime and ignore_day_rising:
                    continue

                if single_culmination.end.is_daytime and ignore_day_set:
                    continue

                row = [
                    str(csv_column.get_value(single_culmination))
                    for csv_column in good_col_defs
                ]
                row_s = ",".join(row) + "\n"
                fp.write(row_s)


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
    )
    print(f"Calculating culminations of {len(obs_sat_fact.active_sats)} satellites")

    day_step = TimeDeltaObj(days=1)
    
    culmination_data = []
    for start_utc, days_to_calc, sat_obs in obs_sat_fact:
        final_day = start_utc + TimeDeltaObj(days=days_to_calc)
        start_end_pairs = list(
            pairwise(make_bounded_time_list(start_utc, final_day, day_step))
        )

        print(
            f"Looking for {sat_obs.sat_name} at "
            f"{sat_obs.obs_name} from {start_utc} to {final_day} " 
            f" with a TLE model whose epoch is {sat_obs.sat_epoch_str}"
        )

        for day_start, day_end in tqdm(start_end_pairs):
            for alt_lim in pargs.alt:
                day_transits = get_day_transits(sat_obs, day_start, day_end, alt_lim)
                culmination_data.extend(day_transits)

    write_to_csv_file(
        f"{pargs.output_file}.csv",
        pargs.csv,
        culmination_data,
        pargs.ignore_daytime,
        pargs.ignore_rise_day,
        pargs.ignore_set_day,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate Satellite rises, sets, and culminations for observatories found in config/obs_data.yaml"
    )

    add_common_params(parser, "culmination_output")

    parser.add_argument(
        "-a",
        "--alt",
        nargs="+",
        type=float,
        default=[30.0],
        help="Change the altitude at which the calculations for rising and falling occur",
    )

    parser.add_argument(
        "-c",
        "--csv",
        type=str,
        default="ALL",
        help="Specifies the section in csv_config.ini to use when building the csv files",
    )

    parser.add_argument(
        "-i",
        "--ignore-daytime",
        action="store_true",
        help="ignore the data points which occur during the day",
    )

    parser.add_argument(
        "-ir",
        "--ignore-rise-day",
        action="store_true",
        help="ignore the data points for which the rising occurs during the day",
    )

    parser.add_argument(
        "-is",
        "--ignore-set-day",
        action="store_true",
        help="ignore the data points for which the setting occurs during the day",
    )

    args = parser.parse_args()

    Path("output").mkdir(parents=True, exist_ok=True)
    main(args)
