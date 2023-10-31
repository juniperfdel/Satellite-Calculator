from skyfield.iokit import Loader


class SkyfieldConstants:
    load = Loader("skyfield_data")
    ephemeris = load("de421.bsp")
    timescale = load.timescale()
    earth = ephemeris["earth"]
    moon = ephemeris["moon"]
    sun = ephemeris["sun"]
