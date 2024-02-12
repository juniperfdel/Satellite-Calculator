from pathlib import Path
from urllib.request import urlopen

from sgp4 import omm
from sgp4.api import Satrec
from skyfield.api import EarthSatellite, load

active = Path("active.csv")

if not active.exists():
    data = urlopen(
        "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=csv"
    )
    active.write_bytes(data.read())


ts = load.timescale()

sats = []
with open(active) as f:
	for fields in omm.parse_csv(f):
		sat = Satrec()
		omm.initialize(sat, fields)
		sats.append(EarthSatellite.from_satrec(sat, ts))

print(sats)
