import sys
import time
import pandas
import shutil
import requests
import yaml

import pytz

from datetime import datetime
from itertools import chain
from pathlib import Path
from typing import Sequence, Union

from tzlocal import get_localzone

from skyfield.sgp4lib import EarthSatellite
from skyfield.toposlib import iers2010

from utils.observatory import Observatories
from utils.satellite import ObservatorySatellite
from utils.skyfield_utils import SkyfieldConstants
from utils.time_utils import TimeDeltaObj, TimeObj, get_off_list, today
from utils.transit import DayTransits

strQ = Union[str, None]
pathQ = Union[Path, None]


def get_day_transits(
    in_obs_sat: ObservatorySatellite,
    start_time: TimeObj,
    end_time: TimeObj,
    input_alt: float,
) -> Sequence[DayTransits]:
    return find_sat_events(in_obs_sat, start_time, end_time, input_alt)


def get_obs_coord_between(
    in_obs_sat: ObservatorySatellite,
    start_time: TimeObj,
    end_time: TimeObj,
    time_step: TimeDeltaObj,
) -> pandas.DataFrame:
    sf_list = get_off_list(start_time, end_time, time_step, 2)
    toc_list = in_obs_sat.diff.at(sf_list)
    ra, dec, _ = toc_list.radec()
    rv_df = {"ra": ra._degrees, "dec": dec._degrees}
    alt, az, _ = toc_list.altaz()
    rv_df["alt"] = alt._degrees
    rv_df["az"] = az._degrees

    subpts = iers2010.subpoint(in_obs_sat.sat.at(sf_list))
    rv_df["sub_lat"] = subpts.latitude._degrees
    rv_df["sub_long"] = subpts.longitude._degrees

    return pandas.DataFrame(data=rv_df)


def find_sat_events(
    in_obs_sat: ObservatorySatellite,
    in_start: TimeObj,
    in_end: TimeObj,
    limit_alt: float,
) -> list[DayTransits]:
    event_times, event_flags = in_obs_sat.sat.find_events(
        in_obs_sat.observatory.sf, in_start.sf, in_end.sf, altitude_degrees=limit_alt
    )

    events_zipped = list(zip(event_times, event_flags))
    rv_transits = []
    for time, flag in events_zipped:
        if flag == 0:
            new_transit = DayTransits(in_obs_sat, limit_alt)
            new_transit.start = TimeObj(time, in_local_tz=in_start.local_tz)
            rv_transits.append(new_transit)

    for time, flag in events_zipped:
        if flag == 2:
            test_time_obj = TimeObj(time, in_local_tz=in_start.local_tz)
            test_i = 0
            while test_i < len(rv_transits):
                if (
                    (rv_transits[test_i].end is None)
                    or (rv_transits[test_i].end.time < test_time_obj)
                ) and (rv_transits[test_i].start.time < test_time_obj):
                    rv_transits[test_i].end = test_time_obj
                    test_i = 0
                else:
                    test_i += 1

    for test_trans in rv_transits:
        for time, flag in events_zipped:
            if flag == 1:
                test_trans.add(time)

    return rv_transits


class ObservatorySatelliteFactory:
    def __init__(
        self,
        local_tz: str = "local",
        reload_sat: bool = True,
        cache_sat: bool = False,
        use_all: bool = False,
        start_date: strQ = "today",
        ignore_limit: bool = False,
        tles: Union[None, list[str]] = None,
    ):
        self.start_date: str = start_date
        self.local_tz: str = local_tz
        self.today_utc: TimeObj = today()
        self.start_utc: Union[TimeObj, None] = None
        self.tles: list[str] = (
            ["https://celestrak.com/NORAD/elements/active.txt"]
            if tles is None
            else tles
        )
        self.reload_sat: bool = reload_sat
        self.cache_sat: bool = cache_sat
        self.use_all: bool = use_all
        self.ignore_limit: bool = ignore_limit
        self.active_sats: list[EarthSatellite] = []
        self.active_observatories: list[Observatories] = []

        self._make_start_utc()
        self._make_active_sats()
        self._make_active_obs()

    def _make_start_utc(self):
        try:
            self.start_utc = TimeObj(
                datetime.strptime(self.start_date, "%Y-%m-%d")
            ).get_start_day()
        except ValueError:
            self.start_utc = self.today_utc.get_start_day()

        self.start_utc.set_local_timezone(
            (
                get_localzone()
                if self.local_tz == "local"
                else pytz.timezone(self.local_tz)
            ).zone
        )

    def _make_active_sats(self):
        self.active_sats = list(
            chain.from_iterable(
                self._make_locally_active_sats(local_tle) for local_tle in self.tles
            )
        )

    def _parse_hist_tle_file(self, in_file: str) -> list[EarthSatellite]:
        itnames_path = Path(SkyfieldConstants.load.path_to("ids_to_names.txt"))
        itnames: dict[str, str] = {}
        if itnames_path.exists():
            itnames: dict[str, str] = dict(map(lambda x: x.split(","), itnames_path.read_text().splitlines(keepends=False)))

        tle_lines = list(filter(Path(SkyfieldConstants.load.path_to(in_file)).read_text().splitlines(keepends=False)))
        tle_curs = 0
        tle_len = len(tle_lines)
        rv: list[EarthSatellite] = []
        while tle_curs < tle_len:
            sat_name = None
            line_one = None
            line_two = None
            c_line = tle_lines[tle_curs]
            if c_line[0] != "1":
                sat_name = c_line
                tle_curs += 1
                c_line = tle_lines[tle_curs]
            if c_line[0] == "1":
                line_one = c_line
                tle_curs += 1
                c_line = tle_lines[tle_curs]
            if c_line[0] == "2":
                line_two = c_line
                tle_curs += 1

            if line_one is None or line_two is None:
                raise ValueError(f"Unable to parse file as TLE file: {in_file} @ line {tle_curs}")
            
            if sat_name is None:
                sat_id = line_two.split(maxsplit=2)[1]
                if sat_id in itnames:
                    sat_name = itnames[sat_id]
                else:
                    print(f"Performing lookup of Catalog number {sat_id}")
                    time.sleep(2.)
                    sat_name = requests.get(f"https://celestrak.org/NORAD/elements/gp.php?CATNR={sat_id}").content.decode("ascii").split()[0]
                    print(f"Found as {sat_name}\n")
                    itnames[sat_id] = sat_name
            
            rv.append(EarthSatellite(line_one, line_two, name=sat_name))

        itnames_path.write_text("\n".join(f"{_sat_id},{_sat_name}" for _sat_id, _sat_name in itnames.items()))
        return rv

    def _make_locally_active_sats(self, sat_url: str = "") -> list[EarthSatellite]:
        is_url = "https://" in sat_url

        try:
            if self.reload_sat:
                loaded_satellites = SkyfieldConstants.load.tle_file(sat_url, reload=True)
            else:
                loaded_satellites = SkyfieldConstants.load.tle_file(sat_url)
        except ValueError as err:
            # Check if Historical TLE
            if not is_url:
                loaded_satellites = self._parse_hist_tle_file(sat_url)
            else:
                raise err

        if self.cache_sat and is_url:
            active_file = Path(SkyfieldConstants.load.path_to("active.txt"))
            cache_loc = Path(
                SkyfieldConstants.load.path_to(
                    f"{self.today_utc.get_file_format()}_sats.txt"
                )
            )
            shutil.copy(active_file, cache_loc)

        satellites = loaded_satellites.copy()
        if not self.use_all:
            with open("config/sat_list.txt", "r") as fp:
                sat_names = [x.strip() for x in fp.readlines()]

            satellites = []
            satellites = [
                sat
                for sat in loaded_satellites
                if (self.ignore_limit or ((self.start_utc.sf - sat.epoch) < 14))
                and any(sat_name in sat.name for sat_name in sat_names)
            ]

        return satellites

    def _make_active_obs(self):
        with open("config/obs_data.yaml", "r") as fp:
            obs_data = yaml.load(fp, Loader=yaml.SafeLoader)

        for obs_name in obs_data:
            obs_lat, obs_long, obs_ele = obs_data[obs_name]
            self.active_observatories.append(
                Observatories(
                    float(obs_lat), float(obs_long), float(obs_ele), str(obs_name)
                )
            )

    def __iter__(self):
        return chain.from_iterable(
            (
                (
                    ObservatorySatellite(current_obs, satellite)
                    for satellite in self.active_sats
                )
                for current_obs in self.active_observatories
            )
        )
