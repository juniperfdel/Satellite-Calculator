# Satellite Calculator 
A set of scripts whose goal is to calculate satellite positions given various conditions.


## Find Satellite Events (find_satellite_events.py)
This is a script which calculates the various characteristics of a satellite when it rises, culminates, and falls.

### File Inputs
#### obs_data.yaml
This file is in the `YAML` format, and specifies the list of observatories to calculate with respect to. The elements are (`<latitude in deg>`, `<longitude in deg>`, `<altitude in meters>`) with the keys being the observatory name.

#### sat_list.txt
This file lists the names of the satellites to look for on each line.

#### csv_config.ini
This file specifies the columns to activate for the CSV file and in what order the columns should be in, just list the column headers in a new section.

### Script Usage
`find_satellite_events.py -h` for detailed help

## Get Satellite Positions (get_satellite_pos.py)
This script calculates the positions of the given satellites wrt the observatories and if those satellites go within a specified radius around a specified coordinate, then those coordinates and times will be written into a `.hdf5` file in `output/sat_pos`. 

### Note
This script accomplishes it's task through calculating the positon of the satellite for every second of every day you ask of it. This means that it will run slow, so even if the progress bar isn't going up; that doesn't mean it stopped. Also before the file is written all the data is stored in RAM, so watch out for running out of RAM.

### File Inputs

#### obs_data.yaml
This file is in the `YAML` format, and therefore should be human-readable, specifies the list of observatories to calculate with respect to. The elements are (`<latitude in deg>`, `<longitude in deg>`, `<altitude in meters>`) with the keys being the observatory name.

#### sat_list.txt
This file lists the names of the satellites to look for on each line.

### Script Usage
`get_satellite_pos.py -h` for detailed help

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
pandas
astropy
dateparser
h5py
```
