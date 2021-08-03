# CALIPSO Calculator
A script whose goal is to calculate the culmination times of CALIPSO given various conditions. This repository is a WIP, future scripts may be added. Pull Requests are welcomed.


## Find CALIPSO (find_calipso.py)

This is the basic script which finds CALIPSO. In order to edit the inputs, the script must be manually edited. Outputs in a .csv file called `<observatory name>.csv`

### Inputs
#### obs_data.yaml
This file is in the `YAML` format, and therefore should be human-readable, specifies the list of observatories to calculate with respect to. The elements are (`<latitude in deg>`, `<longitude in deg>`, `<altitude in meters>`) with the keys being the observatory name.

#### Command Line Inputs
`python3 find_calipso.py [<days> OPTIONS]`

All arguments are optional. The one positional argument, `<days>` determines how many days into the future to calculate the culminations for. By default it is 120.

Here is the list of non-positional arguments:
* `-h`, `--help` will print the help message
* `-tz`, `--set-timezone` will set the timezone to use in the output with you can find a list of possible inputs [here](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones). By default it uses the one your computer is set to.
* `-r`, `--reload` will re-download the TLE file from https://celestrak.com/NORAD/elements/active.txt as opposed to using a cached version.
* `-s`, `--solve` will calculate the satellite positions around culmination based on the next four parameters which are `<start time in seconds before culmination> <duration in seconds> <timestep in seconds> <alt/az flag>` and are all required. If `<alt/az flag>` is set to 1, then alt/az data will be added. Setting any other number will prevent this.

#### HDF5 Note
If you set `-s` or `--store` with the appropriate arguments, then script will calculate the positions of the satellite and place them inside an HDF5 file as individual tables, one for each culmination. Each table is referenced by a key which is the time of the first element in [ISO-8601 format](https://en.wikipedia.org/wiki/ISO_8601) The [wikipedia page for HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) has a list of languages which have HDF5 parsers for them already written.

### Python Version
This script uses features from the `typing` module which was introduced in python 3.7, so that is the minimum python version required to run this script.

### Dependencies
* tqdm
* pytz
* tzlocal
* numpy
* skyfield
* geopy
* pyyaml
* h5py (If storing Satillite Positions)
