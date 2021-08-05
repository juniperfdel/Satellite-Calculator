# CALIPSO Calculator
A set of scripts whose goal is to calculate satellite positions given various conditions. This repository is a WIP, future scripts may be added. Pull Requests are welcomed.


## Find Satellite Events (find_satellite_events.py)

This is a script which calculates the various charateristics of a satellite when it rises, culminates, and falls . Outputs in a .csv file called `<sat name>_<obs name>_<today>-<final day>.csv`

### Inputs
#### obs_data.yaml
This file is in the `YAML` format, and therefore should be human-readable, specifies the list of observatories to calculate with respect to. The elements are (`<latitude in deg>`, `<longitude in deg>`, `<altitude in meters>`) with the keys being the observatory name.

#### sat_list.txt
This file lists the names of the satellites to look for on each line.

#### csv_config.ini
This file specifies the columns to activate for the CSV file and in what order the columns should be in

#### Command Line Inputs
`python3 find_calipso.py [<days> OPTIONS]`

All arguments are optional. The one positional argument, `<days>` determines how many days into the future to calculate the culminations for. By default it is 120.

Here is the list of non-positional arguments:
* `-h`, `--help` will print the help message
* `-tz`, `--set-timezone` will set the timezone to use in the output with you can find a list of possible inputs [here](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones). By default it uses the one your computer is set to.
* `-r`, `--reload` will re-download the TLE file from https://celestrak.com/NORAD/elements/active.txt as opposed to using a cached version.
* `-s`, `--solve` will calculate the satellite positions around culmination based on the next four parameters which are `<start time in seconds before culmination> <duration in seconds> <timestep in seconds> <alt/az flag>` and are all required. If `<alt/az flag>` is set to 1, then alt/az data will be added. Setting any other number will prevent this.
* `-a`, `--alt` will set the altitude at which the calculations for rising and falling occur; by default it will be 30.
* `-c`, `--csv` will specify the section in `csv_config.ini` to use. default is "DEFAULT"

## Get Satellite Positions (get_satellite_pos.py)
This script calculates the positions of the given satellites wrt the observatories and if those satellites go within a specified radius around a specified coordinate, then those coordinates and times will be written into a `.hdf5` file in `output/sat_pos`. 

### Note
This script accompishes it's task through calculating the positon of the satellite for every second of every day you ask of it. This means that it will run slow, so even if the progress bar isn't going up; that doesn't mean it stopped. Also before the file is written all the data is stored in RAM, so watch out for running out of RAM.

#### obs_data.yaml
This file is in the `YAML` format, and therefore should be human-readable, specifies the list of observatories to calculate with respect to. The elements are (`<latitude in deg>`, `<longitude in deg>`, `<altitude in meters>`) with the keys being the observatory name.

#### sat_list.txt
This file lists the names of the satellites to look for on each line.

#### Command Line Inputs
`python3 find_calipso.py <days> <Coordinate 1> <Coordinate 2> <Radius> <Flag> [OPTIONS]`

All positional arguments are required.

Here is the list of positional arguments:
* `days` - The number of days to calculate into the future
* `Coordinate 1` - The search area's center's first coordinate either RA or Alt depending on `<Flag>`
* `Coordinate 2` - The search area's center's second coordinate either Dec or Az depending on `<Flag>`
* `Radius` - The radius of the search area
* `Flag` - Indicates what type of positon is being specified where 0 = RA/Dec and 1 = Alt/Az

Here is the list of non-positional arguments:
* `-h`, `--help` will print the help message
* `-tz`, `--set-timezone` will set the timezone to use in the output with you can find a list of possible inputs [here](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones). By default it uses the one your computer is set to.
* `-r`, `--reload` will re-download the TLE file from https://celestrak.com/NORAD/elements/active.txt as opposed to using a cached version.

## Python Version
These scripts uses features from the `typing` module which was introduced in python 3.7, so that is the minimum python version required to run this script.

## HDF5 Note
If any script calculates the positions of a satellite, it will place those positions inside an individual table in an HDF5 file. Each table corresponds to a culmination or day (depending on the script). Each table is referenced by a key which is the time of the first element in [ISO-8601 format](https://en.wikipedia.org/wiki/ISO_8601) The [wikipedia page for HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) has a list of languages which have HDF5 parsers for them already written.

## Dependencies
* tqdm
* pytz
* tzlocal
* numpy
* skyfield
* geopy
* pyyaml
* h5py (If storing Satillite Positions)
