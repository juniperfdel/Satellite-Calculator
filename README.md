Please go to the public version of this repo for more up-to-date code https://github.com/gregoryfdel/Satellite-Calculator

# Satellite Calculator (Private)
A set of scripts whose goal is to calculate satellite positions given various conditions. This repository is a WIP, future scripts may be added. Pull Requests are welcomed.


## Find Satellite Events (find_satellite_events.py)

This is a script which calculates the various characteristics of a satellite when it rises, culminates, and falls . Outputs in a .csv file called `culmination_output.csv`

### File Inputs
#### obs_data.yaml
This file is in the `YAML` format, and specifies the list of observatories to calculate with respect to. The elements are (`<latitude in deg>`, `<longitude in deg>`, `<altitude in meters>`) with the keys being the observatory name.

#### sat_list.txt
This file lists the names of the satellites to look for on each line.

#### csv_config.ini
This file specifies the columns to activate for the CSV file and in what order the columns should be in, just list the column headers in a new section.

### Command Line Inputs
```
usage: find_satellite_events.py [-h] [-tz TIMEZONE] [-r] [-ca] [-l] [-sd START_DATE] [-il] [--tles TLES [TLES ...]] [-o OUTPUT_FILE] [-a ALT [ALT ...]] [-c CSV] [-i] [-ir] [-is] [days]

Calculate Satellite rises, sets, and culminations for observatories found in config/obs_data.yaml

positional arguments:
  days                  The number of days to calculate for, default of 120 days

options:
  -h, --help            show this help message and exit
  -tz TIMEZONE, --timezone TIMEZONE
                        The timezone to calculate with respect to, default is your local timezone
  -r, --reload          Re-download the TLE file from the internet
  -ca, --cache          Automatically Cache TLE files
  -l, --all             Use all satellites from database
  -sd START_DATE, --start-date START_DATE
                        specify a starting date to begin calculations format is <year>-<month>-<day>
  -il, --ignore-limit   ignore the 14 day recommended TLE limit on satellite data
  --tles TLES [TLES ...]
                        Specify tle files to use instead of downloading them
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        name of the file to output
  -a ALT [ALT ...], --alt ALT [ALT ...]
                        Change the altitude at which the calculations for rising and falling occur
  -c CSV, --csv CSV     Specifies the section in csv_config.ini to use when building the csv files
  -i, --ignore-daytime  ignore the data points which occur during the day
  -ir, --ignore-rise-day
                        ignore the data points for which the rising occurs during the day
  -is, --ignore-set-day
                        ignore the data points for which the setting occurs during the day
```

## Get Satellite Positions (get_satellite_pos.py)
This script calculates the positions of the given satellites wrt the observatories and if those satellites go within a specified radius around a specified coordinate, then those coordinates and times will be written into a `.hdf5` file in `output/sat_pos`. 

### Note
This script accomplishes it's task through calculating the positon of the satellite for every second of every day you ask of it. This means that it will run slow, so even if the progress bar isn't going up; that doesn't mean it stopped. Also before the file is written all the data is stored in RAM, so watch out for running out of RAM.

### File Inputs

#### obs_data.yaml
This file is in the `YAML` format, and therefore should be human-readable, specifies the list of observatories to calculate with respect to. The elements are (`<latitude in deg>`, `<longitude in deg>`, `<altitude in meters>`) with the keys being the observatory name.

#### sat_list.txt
This file lists the names of the satellites to look for on each line.

### Command Line Inputs
```
usage: get_satellite_pos.py [-h] [-tz TIMEZONE] [-r] [-ca] [-l] [-sd START_DATE] [-il] [--tles TLES [TLES ...]] [-o OUTPUT_FILE] [-m] [-i] [--step STEP] coord1 coord2 radius flag [days]

search the inputted satellites track for how close each will get to your input.

positional arguments:
  coord1                RA or Alt depending on <flag>
  coord2                Dec or Az depending on <flag>
  radius                The radius around the coordinate
  flag                  0 = RA/Dec; 1 = Alt/Az; 2 = ignore
  days                  The number of days to calculate for, default of 120 days

options:
  -h, --help            show this help message and exit
  -tz TIMEZONE, --timezone TIMEZONE
                        The timezone to calculate with respect to, default is your local timezone
  -r, --reload          Re-download the TLE file from the internet
  -ca, --cache          Automatically Cache TLE files
  -l, --all             Use all satellites from database
  -sd START_DATE, --start-date START_DATE
                        specify a starting date to begin calculations format is <year>-<month>-<day>
  -il, --ignore-limit   ignore the 14 day recommended TLE limit on satellite data
  --tles TLES [TLES ...]
                        Specify tle files to use instead of downloading them
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        name of the file to output
  -m, --merge           write all data into a single file
  -i, --ignore-empty    If the dataset is empty, it will not be written to the file
  --step STEP           Step size in seconds
```
## Dependencies
### Python Version
These scripts uses features from the `typing` module which was introduced in python 3.7, so that is the minimum recommended python version required to run these scripts. 
### Formatting
All scripts have been formatted with [black](https://github.com/psf/black).
### Packages
This repository contains a `pipfile` which is used by the dependency management system [pipenv](https://pipenv.pypa.io/en/latest/). Otherwise, here is a list of the packages needed to run these scripts
```
pyyaml
tqdm
pytz
numpy
geopy
skyfield
tzlocal
h5py
cached-property
pandas
astropy
dateparser
h5py
```