from __future__ import annotations

from datetime import datetime
from typing import Tuple, Union

from geopy import distance
from skyfield.almanac import dark_twilight_day, fraction_illuminated
from skyfield.sgp4lib import EarthSatellite
from skyfield.toposlib import GeographicPosition, iers2010


from utils.math_utils import two_object_distance
from utils.struct_utils import MetaFormatter
from utils.observatory import Observatories
from utils.skyfield_utils import SkyfieldConstants
from utils.time_utils import TimeObj


class ObservatorySatellite(metaclass=MetaFormatter):
    def __init__(self, in_observatory: Observatories, satellite_: EarthSatellite):
        self.observatory = in_observatory
        self.sat = satellite_
        self.diff = self.sat - self.observatory.sf
        self.sat_name = self.sat.name
        self.obs_name = self.observatory.name

    def at(self, in_time: TimeObj):
        return SatellitePosition(self, in_time)

    @property
    def sat_epoch_obj(self) -> TimeObj:
        return TimeObj(self.sat.epoch)

    @property
    def sat_epoch_str(self) -> str:
        return TimeObj(self.sat.epoch).utc_formatted_str()

    @property
    def sat_epoch_mjd(self) -> float:
        return (self.sat.model.jdsatepoch - 2400000.5)

    @property
    def sat_mean_motion(self) -> float:
        return self.sat.model.no_kozai


class SatellitePosition(metaclass=MetaFormatter):
    def __init__(self, in_obs_sat: ObservatorySatellite, in_time: TimeObj):
        self.obs: Observatories = in_obs_sat.observatory
        self.sat: EarthSatellite = in_obs_sat.sat

        self.twilight_checker = dark_twilight_day(
            SkyfieldConstants.ephemeris, self.obs.sf
        )
        self.time = in_time

        self._sat_toc = in_obs_sat.diff.at(self.time.sf)

    @property
    def local_tz(self) -> Union[None, str]:
        return self.time.local_tz

    @property
    def local_datetime(self) -> datetime:
        return self.time.get_local_time()

    @property
    def utc_datetime(self) -> datetime:
        return self.time.dt

    @property
    def utc_str(self) -> str:
        return self.time.utc_formatted_str()

    @property
    def sec_since_midnight(self) -> float:
        return (self.time.get_mjd() - int(self.time.get_mjd())) * 86400

    @property
    def local_time_str(self) -> str:
        return self.time.local_formatted_str()

    @property
    def obs_name(self) -> str:
        return self.obs.name

    @property
    def obs_lat(self) -> float:
        return self.obs.lat

    @property
    def obs_long(self) -> float:
        return self.obs.long

    @property
    def sat_name(self) -> str:
        return self.sat.name

    @property
    def sat_epoch(self) -> TimeObj:
        return TimeObj(self.sat.epoch)

    @property
    def sat_epoch_time_diff_days(self) -> float:
        return self.time.get_mjd() - (self.sat.model.jdsatepoch - 2400000.5)

    @property
    def sat_ra_dec(self) -> Tuple[float, float]:
        sat_ra, sat_dec, _ = self._sat_toc.radec(
            epoch=SkyfieldConstants.timescale.J2000
        )
        return sat_ra._degrees, sat_dec._degrees

    @property
    def sat_alt_az(self) -> Tuple[float, float]:
        sat_alt, sat_az, _ = self._sat_toc.altaz()
        return sat_alt._degrees, sat_az._degrees

    @property
    def sat_ra(self) -> float:
        ra, _ = self.sat_ra_dec
        return ra

    @property
    def sat_dec(self) -> float:
        _, dec = self.sat_ra_dec
        return dec

    @property
    def sat_alt(self) -> float:
        alt, _ = self.sat_alt_az
        return alt

    @property
    def sat_az(self) -> float:
        _, az = self.sat_alt_az
        return az

    @property
    def moon_ra_dec(self) -> Tuple[float, float]:
        moon_ra, moon_dec, _ = (
            self.obs.sf_vector_sum.at(self.time.sf)
            .observe(SkyfieldConstants.moon)
            .apparent()
            .radec(epoch=SkyfieldConstants.timescale.J2000)
        )

        return moon_ra._degrees, moon_dec._degrees

    @property
    def moon_ra(self) -> float:
        moon_ra, _ = self.moon_ra_dec
        return moon_ra

    @property
    def moon_dec(self) -> float:
        _, moon_dec = self.moon_ra_dec
        return moon_dec

    @property
    def moon_fraction(self) -> float:
        return fraction_illuminated(SkyfieldConstants.ephemeris, "moon", self.time.sf)

    @property
    def moon_distance(self) -> float:
        return two_object_distance(
            self.sat_ra, self.sat_dec, self.moon_ra, self.moon_dec
        )

    @property
    def is_daytime(self) -> bool:
        return self.twilight_checker(self.time.sf) > 1

    @property
    def is_daytime_str(self) -> str:
        return "Y" if self.is_daytime else "N"

    @property
    def is_sat_sunlit(self) -> bool:
        return self.sat.at(self.time.sf).is_sunlit(SkyfieldConstants.ephemeris)

    @property
    def is_sat_sunlit_str(self) -> str:
        return "Y" if self.is_sat_sunlit else "N"

    @property
    def sat_geo_position(self) -> GeographicPosition:
        return iers2010.subpoint(self.sat.at(self.time.sf))

    @property
    def sat_lat_long(self) -> Tuple[float, float]:
        return (
            self.sat_geo_position.latitude._degrees,
            self.sat_geo_position.longitude._degrees,
        )

    @property
    def sat_lat(self) -> float:
        geo_lat, _ = self.sat_lat_long
        return geo_lat

    @property
    def sat_long(self) -> float:
        _, geo_long = self.sat_lat_long
        return geo_long

    @property
    def sat_obs_distance(self) -> float:
        return distance.distance(self.sat_lat_long, self.obs.lat_long).km
