import argparse
import configparser
import json
import os
import pathlib
from operator import attrgetter
from typing import List

import pytz
from tqdm import tqdm
from tzlocal import get_localzone

from utils import DayTransits, SingleCulmination, build_obs_sat, get_day_transits, today


class ColumnDefinition:
    """Creates the column name and value getter from a definition"""

    def __init__(self, *column_definition):
        self.column_name = column_definition[0]
        self.attr_getter = attrgetter(".".join(column_definition[1:]))
        self.index = -1

    def get_ini_key(self):
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

    good_col_defs = []
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


def main(in_args) -> None:
    today_utc = today()
    today_utc.set_local_timezone(
        (
            get_localzone()
            if in_args.timezone == "local"
            else pytz.timezone(in_args.timezone)
        ).zone
    )

    sat_obs_objects = build_obs_sat(in_args.reload, in_args.cache, in_args.all)
    final_day = today_utc.get_off(days=in_args.days)

    print(f"Calculating culminations of {len(sat_obs_objects)} satellites")

    culmination_data = []
    for sat_obs in sat_obs_objects:
        print(
            f"Looking for {sat_obs.sat_name} at {sat_obs.obs_name} from {today_utc} to {final_day}"
        )

        for day_num in tqdm(range(in_args.days)):
            for alt_lim in in_args.alt:
                day_start = today_utc.get_off(days=day_num).get_start_day()
                day_end = day_start.get_off(days=1)

                day_transits = get_day_transits(sat_obs, day_start, day_end, alt_lim)
                culmination_data.extend(day_transits)

    write_to_csv_file(
        f"{in_args.output_file[0]}.csv",
        in_args.csv,
        culmination_data,
        in_args.ignore_daytime,
        in_args.ignore_rise_day,
        in_args.ignore_set_day,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate Satellite rises, sets, and culminations for observatories found in config/obs_data.yaml"
    )

    parser.add_argument(
        "days",
        type=int,
        nargs="?",
        default=120,
        help="The number of days to calculate for, default of 120 days",
    )

    parser.add_argument(
        "-tz",
        "--timezone",
        type=str,
        default="local",
        help="The timezone to calculate with respect to, default is your local timezone",
    )

    parser.add_argument(
        "-r",
        "--reload",
        action="store_true",
        help="Re-download the TLE file from the internet",
    )
    parser.add_argument(
        "-ca", "--cache", action="store_true", help="Automatically Cache TLE files"
    )
    parser.add_argument(
        "-a",
        "--alt",
        nargs="+",
        type=float,
        default=30.0,
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

    parser.add_argument(
        "-o", "--output-file", nargs=1, type=str, default="culmination_output"
    )

    parser.add_argument(
        "-l", "--all", action="store_true", help="Use all satellites from database"
    )

    args = parser.parse_args()

    pathlib.Path(os.path.join("output")).mkdir(parents=True, exist_ok=True)
    main(args)
