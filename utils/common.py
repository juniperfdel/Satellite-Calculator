from datetime import datetime
from itertools import chain
from pathlib import Path
from typing import Sequence, Union

import sys
import pandas
import shutil
import yaml

import pytz
import numpy as np
from tzlocal import get_localzone

from skyfield.sgp4lib import EarthSatellite

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
    rv_df = get_off_list(start_time, end_time, time_step, 4)

    toc_list = in_obs_sat.diff.at(sf_list)
    ra, dec, _ = toc_list.radec()
    rv_df["ra"] = ra._degrees
    rv_df["dec"] = dec._degrees

    alt, az, _ = toc_list.altaz()
    rv_df["alt"] = alt._degrees
    rv_df["az"] = az._degrees

    return rv_df


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
            if tles is None else tles
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
        self.start_utc = (
            TimeObj(datetime.strptime(self.start_date, "%Y-%m-%d"))
            if bool(self.start_date)
            else self.today_utc
        ).get_start_day()

        self.start_utc.set_local_timezone(
            (
                get_localzone()
                if self.local_tz == "local"
                else pytz.timezone(self.local_tz)
            ).zone
        )

    def _make_active_sats(self):
        self.active_sats = list(chain.from_iterable(
            self._make_locally_active_sats(local_tle)
            for local_tle in self.tles
        ))

    def _make_locally_active_sats(self, sat_url: str = "") -> list[EarthSatellite]:
        is_url = "https://" in sat_url

        if self.reload_sat:
            loaded_satellites = SkyfieldConstants.load.tle_file(sat_url, reload=True)
        else:
            loaded_satellites = SkyfieldConstants.load.tle_file(sat_url)

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
                if (self.ignore_limit or ((self.start_utc.sf - sat.epoch) < 14)) and 
                any(sat_name in sat.name for sat_name in sat_names)
            ]
        
        return satellites

    def _make_active_obs(self):
        with open("config/obs_data.yaml", "r") as fp:
            obs_data = yaml.load(fp, Loader=yaml.SafeLoader)
        
        for obs_name in obs_data:
            obs_lat, obs_long, obs_ele = obs_data[obs_name]
            self.active_observatories.append(Observatories(obs_lat, obs_long, obs_ele, obs_name))

    def __iter__(self):
        return chain.from_iterable(((
                ObservatorySatellite(current_obs, satellite)
                for satellite in self.active_sats
            ) for current_obs in self.active_observatories)
        )