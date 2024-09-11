from __future__ import annotations

from datetime import datetime, timedelta
from enum import auto, Enum
from typing import Any, Union
from zoneinfo import ZoneInfo

import numpy as np
import pandas
from astropy.time import Time as AstropyTime
from astropy.time import TimeDelta as AstropyTimeDelta
from astropy.timeseries import TimeSeries as AstropyTimeSeries
from numpy.typing import ArrayLike
from pandas import Timedelta as PandasTimeDelta
from pandas import Timestamp as PandasTimeStamp
from skyfield.timelib import Time as SkyfieldTime

from utils.skyfield_utils import SkyfieldConstants


def timezone_converter(
    input_dt: datetime,
    current_tz: str = "UTC",
    target_tz: str = "UTC",
) -> datetime:
    return input_dt.replace(tzinfo=ZoneInfo(current_tz)).astimezone(ZoneInfo(target_tz))


class TimeObj:
    def __init__(self, input_time: Any, local_tz: str = "UTC") -> None:
        """
        Create a custom time object which is used to switch between
        inputs to the various packages used, as well as simplify
        common time associated functions

        :param input_time:
        """
        self.dt = None
        self.local_tz: str = local_tz

        if isinstance(input_time, TimeObj):
            self.dt = input_time.dt
            self.local_tz = input_time.local_tz
        elif isinstance(input_time, datetime):
            self.dt = timezone_converter(input_time)
        elif isinstance(input_time, SkyfieldTime):
            self.dt = input_time.utc_datetime()
        elif isinstance(input_time, np.datetime64):
            self.dt = timezone_converter(input_time.astype(datetime))
        elif isinstance(input_time, AstropyTime):
            self.dt = timezone_converter(input_time.to_datetime())
        elif isinstance(input_time, PandasTimeStamp):
            self.dt = timezone_converter(input_time.to_pydatetime())

        if self.dt is None:
            raise TypeError(
                "Only datetime, numpy.datetime64, skyfield.timelib.Time, "
                "pandas.Timestamp or astropy.time.Time objects are allowed!, "
                f"inputTime was of type {type(input_time)}"
            )

        self.sf = SkyfieldConstants.timescale.from_datetime(self.dt)
        self.np = np.datetime64(self.dt)
        self.ap = AstropyTime(self.dt)
        self.pd = PandasTimeStamp(self.dt)

    def __repr__(self) -> str:
        return (
            "TimeObj(datetime=%r, skyfield.timelib.Time=%r, numpy.datetime64=%r, astropy.time.Time=%r, pandas.Timestamp=%r)"
            % (self.dt, self.sf, self.np, self.ap, self.pd)
        )

    def __str__(self) -> str:
        return str(self.dt)

    def __add__(self, other: "TimeDeltaObj") -> TimeObj:
        return TimeObj(self.dt + other.dt, local_tz=self.local_tz)

    def __sub__(self, other: TimeObj) -> "TimeDeltaObj":
        return TimeDeltaObj(self.dt - other.dt)

    def __ge__(self, other: TimeObj) -> bool:
        return self.dt >= other.dt

    def __gt__(self, other: TimeObj) -> bool:
        return self.dt > other.dt

    def __lt__(self, other: TimeObj) -> bool:
        return self.dt < other.dt

    def __le__(self, other: TimeObj) -> bool:
        return self.dt <= other.dt

    def __eq__(self, other: TimeObj) -> bool:
        return self.dt == other.dt

    def get_start_day(self):
        return TimeObj(
            datetime.combine(self.dt.date(), datetime.min.time()),
            local_tz=self.local_tz,
        )

    def get_off(self, **kwargs: int):
        return TimeObj(self.dt + timedelta(**kwargs), local_tz=self.local_tz)

    def get_off_center(self, **kwargs: int):
        time_d = timedelta(**kwargs) / 2
        return TimeObj(self.dt - time_d, local_tz=self.local_tz), TimeObj(
            self.dt + time_d, local_tz=self.local_tz
        )

    def get_datetime(self) -> datetime:
        return self.dt

    def get_skyfield(self) -> SkyfieldTime:
        return self.sf

    def get_numpy(self) -> np.datetime64:
        return self.np

    def get_astropy(self) -> AstropyTime:
        return self.ap

    def get_pandas(self) -> PandasTimeStamp:
        return self.pd

    def iso_format(self) -> str:
        return self.dt.isoformat()

    def get_file_format(self) -> str:
        return self.dt.strftime("%Y%m%d")

    def get_compact_fmt(self) -> str:
        return self.dt.strftime("%Y%m%d %H:%M")

    def get_jd(self) -> float:
        return self.ap.jd

    def get_mjd(self) -> float:
        return self.ap.mjd

    def get_local_time(self) -> datetime:
        return timezone_converter(self.dt, target_tz=self.local_tz)

    def utc_formatted_str(self) -> str:
        return self.dt.strftime("%Y %b %d %H:%M:%S")

    def local_formatted_str(self) -> str:
        return self.get_local_time().strftime("%Y %b %d %H:%M:%S")


class TimeDeltaObj:
    def __init__(self, input_delta=None, **delta_args) -> None:
        """
        Create a custom time delta object which is used to switch between
        inputs to the various packages used, as well as simplify
        common time associated functions

        :param input_delta:
        """
        if delta_args:
            self.dt = timedelta(**delta_args)
        elif isinstance(input_delta, timedelta):
            self.dt = input_delta
        elif isinstance(input_delta, np.timedelta64):
            self.dt = input_delta.astype(timedelta)
        elif isinstance(input_delta, AstropyTimeDelta):
            self.dt = input_delta.to_datetime()
        elif isinstance(input_delta, PandasTimeDelta):
            self.dt = input_delta.to_pytimedelta()
        else:
            raise TypeError(
                "Input delta keyword args, a datetime, numpy.timedelta64, "
                "pandas.Timedelta or astropy.time.TimeDelta! "
                f"inputTime was of type {type(input_delta)} and "
                f"inputted keyword args were {delta_args}"
            )

        self.np = np.timedelta64(self.dt)
        self.pd = pandas.to_timedelta(self.dt)
        self.ap = AstropyTimeDelta(self.dt)

    def total_seconds(self) -> float:
        return self.dt.total_seconds()

    def total_days(self) -> float:
        return self.dt.total_seconds() / 86400

    def get_datetime(self) -> timedelta:
        return self.dt

    def get_numpy(self) -> np.timedelta64:
        return self.np

    def get_astropy(self) -> AstropyTimeDelta:
        return self.ap

    def get_pandas(self) -> PandasTimeDelta:
        return self.pd


class OffListTypes(Enum):
    NumpyTL = auto()
    PythonTL = auto()
    SkyfieldTL = auto()
    AstropyTL = auto()
    PandasTL = auto()
    TimeObjTL = auto()


def get_off_list(
    t_start: TimeObj,
    t_end: TimeObj,
    t_step: TimeDeltaObj,
    final_type: OffListTypes = OffListTypes.NumpyTL,
) -> Union[
    AstropyTimeSeries,
    ArrayLike,
    list[datetime],
    SkyfieldTime,
    pandas.DataFrame,
    TimeObj,
]:
    """
    Create a list of times for use by different astronomical packages

    Parameters
    ----------
    t_start
    t_end
    t_step
    final_type

    Returns
    -------
    """
    if final_type == OffListTypes.AstropyTL or final_type == OffListTypes.PandasTL:
        ap_step = t_step.get_astropy()
        n_steps = (t_end.np - t_start.np) // t_step.np + 1
        ap_ts = AstropyTimeSeries(
            time_start=t_start.ap, time_delta=ap_step, n_samples=n_steps
        )
        return ap_ts if final_type == OffListTypes.AstropyTL else ap_ts.to_pandas()

    num_list = np.arange(t_start.np, t_end.np + t_step.np, t_step.np)
    if final_type == OffListTypes.NumpyTL:
        return num_list

    dt_list = [timezone_converter(t) for t in num_list.astype(datetime)]
    if final_type == OffListTypes.PythonTL:
        return dt_list

    if final_type == OffListTypes.SkyfieldTL:
        return SkyfieldConstants.timescale.from_datetimes(dt_list)

    if final_type == OffListTypes.TimeObjTL:
        assert (
            t_start.local_tz == t_end.local_tz
        ), "Start and End should have the same local timezones!"
        return [TimeObj(x, local_tz=t_start.local_tz) for x in dt_list]

    raise TypeError(
        "Please select a type 0 (numpy.datetime64), 1 (datetime), 2 (skyfield.timelib.Time), 3 (astropy.time.Time), 4 (pandas.DataFrame), or 5 (TimeObj)! "
    )


def today() -> TimeObj:
    dt_today = datetime.combine(datetime.now().date(), datetime.min.time())
    return TimeObj(dt_today)
