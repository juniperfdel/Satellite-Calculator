import json
import shutil
import time
from enum import Enum
from pathlib import Path
from urllib.request import urlopen

import requests
from sgp4 import omm
from sgp4.api import Satrec
from skyfield.api import EarthSatellite
from skyfield.sgp4lib import EarthSatellite

from utils.skyfield_utils import SkyfieldConstants
from utils.time_utils import TimeObj, today


class FileGuess(Enum):
    csv = ["OBJECT_NAME", ",", "OBJECT_ID"]
    xml = ["<?xml"]
    tle = [""]


def guess_type(file: Path) -> FileGuess:
    with open(file) as f:
        first_line = f.readline()
    file_type = None
    for test in FileGuess:
        if all(map(lambda x: x in first_line, test.value)):
            file_type = test
            break
    return file_type


def lookup_name_sat_id(sat_ids: list[str]) -> dict[str, str]:
    itnames_path = Path(SkyfieldConstants.load.path_to("ids_to_names.txt"))
    itnames: dict[str, str] = {}
    if itnames_path.exists():
        itnames: dict[str, str] = dict(
            map(
                lambda x: x.split(","),
                itnames_path.read_text().splitlines(keepends=False),
            )
        )

    rv: dict[str, str] = {}
    for sat_id in sat_ids:
        sat_id_str = str(sat_id)
        if sat_id_str in itnames:
            sat_name = itnames[sat_id_str]
        else:
            print(f"Performing lookup of Catalog number {sat_id}")
            time.sleep(2.0)
            sat_name = json.loads(
                requests.get(
                    f"https://celestrak.org/satcat/records.php?CATNR={sat_id}"
                ).content.decode("ascii")
            )[0]["OBJECT_NAME"]
            print(f"Found as {sat_name}\n")
            itnames[sat_id_str] = sat_name
        rv[sat_id] = sat_name

    itnames_path.write_text(
        "\n".join(f"{_sat_id},{_sat_name}" for _sat_id, _sat_name in itnames.items())
    )
    return rv


def parse_hist_tle_file(in_file: str) -> list[EarthSatellite]:
    tle_lines = list(
        filter(
            Path(SkyfieldConstants.load.path_to(in_file))
            .read_text()
            .splitlines(keepends=False)
        )
    )
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
            raise ValueError(
                f"Unable to parse file as TLE file: {in_file} @ line {tle_curs}"
            )

        rv.append(EarthSatellite(line_one, line_two, name=sat_name))

    return rv


def csv_xml_handler(active: Path, ext: FileGuess) -> list[EarthSatellite]:
    omm_parsers = {FileGuess.csv: omm.parse_csv, FileGuess.xml: omm.parse_xml}

    sats = []
    with open(active) as f:
        for fields in omm_parsers[ext](f):
            sat = Satrec()
            omm.initialize(sat, fields)
            nsat = EarthSatellite.from_satrec(sat, SkyfieldConstants.timescale)
            nsat.name = fields["OBJECT_NAME"]
            sats.append(nsat)

    return sats


def handle_url(sat_url: str, reload_sat: bool) -> list[EarthSatellite]:

    is_url = sat_url.startswith("http")
    active = (
        Path(SkyfieldConstants.load.path_to("active.txt")) if is_url else Path(sat_url)
    )

    if reload_sat and is_url:
        data = urlopen(sat_url)
        active.write_bytes(data.read())

    guess = guess_type(active)

    return (
        SkyfieldConstants.load.tle_file(active, reload=reload_sat)
        if guess == FileGuess.tle
        else csv_xml_handler(active, guess)
    )


def make_locally_active_sats(
    reload_sat: bool,
    cache_sat: bool,
    use_all: bool,
    ignore_limit: bool,
    start_utc: TimeObj,
    sat_url: str = "",
) -> list[EarthSatellite]:
    loaded_satellites = handle_url(sat_url, reload_sat)

    is_url = sat_url.startswith("http")
    if (loaded_satellites is None) and (not is_url):
        print(
            f"Skyfield unable to load {sat_url}; attempting manual parsing as a TLE file"
        )
        loaded_satellites = parse_hist_tle_file(sat_url)

    if cache_sat and is_url:
        active_file = Path(SkyfieldConstants.load.path_to("active.txt"))
        cache_loc = Path(
            SkyfieldConstants.load.path_to(f"{today().get_file_format()}_sats.txt")
        )
        shutil.copy(active_file, cache_loc)

    lookup_sat_ids = [sat.model.satnum for sat in loaded_satellites if sat.name is None]
    sid_name = lookup_name_sat_id(lookup_sat_ids)
    for sat in loaded_satellites:
        if sat.name is None:
            sat.name = sid_name[sat.model.satnum]

    if not use_all:
        with open("config/sat_list.txt", "r") as fp:
            sat_names = [x.strip() for x in fp.readlines()]

        return [
            sat
            for sat in loaded_satellites
            if (ignore_limit or ((start_utc.sf - sat.epoch) < 14))
            and any(sat_name in sat.name for sat_name in sat_names)
        ]
    else:
        return loaded_satellites
