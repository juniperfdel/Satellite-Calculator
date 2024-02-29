from enum import Enum
from pathlib import Path
from urllib.request import urlopen

from sgp4 import omm
from sgp4.api import Satrec
from skyfield.api import EarthSatellite, load


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


omm_parsers = {FileGuess.csv: omm.parse_csv, FileGuess.xml: omm.parse_xml}

ts = load.timescale()
active = Path("active.temp")

if not active.exists():
    data = urlopen(
        "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=csv"
    )
    active.write_bytes(data.read())


ext = guess_type(active)

sats = []
with open(active) as f:
    for fields in omm_parsers[ext](f):
        sat = Satrec()
        omm.initialize(sat, fields)
        nsat = EarthSatellite.from_satrec(sat, ts)
        nsat.name = fields["OBJECT_NAME"]
        sats.append(nsat)

print(sats)
