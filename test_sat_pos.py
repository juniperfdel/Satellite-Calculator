import pandas as pd
from datetime import datetime

from astropy.time import Time as AstropyTime
from skyfield.sgp4lib import EarthSatellite

from utils.time_utils import TimeObj, TimeDeltaObj
from utils.observatory import Observatories
from utils.satellite import ObservatorySatellite
from utils.common import get_obs_coord_between


def test_satellite_position():
    calipso_sat = EarthSatellite(
        "1 29108U 06016B   21130.88612831  .00000124  00000+0  33840-4 0  9995",
        "2 29108  98.2735  83.5693 0001324  84.0787 276.0565 14.62540865800172",
    )
    veritas_obs = ob = Observatories(31.675, -110.9519, 1280.0, "VERITAS")

    vc_os = ObservatorySatellite(veritas_obs, calipso_sat)
    start_time = TimeObj(datetime(2021, 5, 12, 10, 3, 47))
    end_time = TimeObj(datetime(2021, 5, 12, 10, 3, 55))
    coord_return = get_obs_coord_between(
        vc_os, start_time, end_time, TimeDeltaObj(seconds=1)
    )
    coord_return["time_fmt_str"] = [
        AstropyTime(x, format="mjd").to_datetime().strftime("%Y-%m-%d %H:%M:%S")
        for x in coord_return["mjd"].to_list()
    ]
    coord_return.set_index("time_fmt_str", drop=True, inplace=True)
    test_data = pd.read_csv("test_sat_pos_dat.csv", index_col="time_fmt_str")
    print(coord_return.loc[test_data.index].drop(columns=["alt", "az"]) - test_data)


if __name__ == "__main__":
    test_satellite_position()
