# CALIPSO Calculator
A script whose goal is to calculate the culmination times of CALIPSO given various conditions. This repository is a WIP, future scripts may be added. Pull Requests are welcomed.


### Find CALIPSO (find_calipso.py)

This is the basic script which finds calipso. In order to edit the inputs, the script must be manually edited. Outputs in a .csv file called `<observatory name>_CALIPSO_Info.csv`

#### Storing Satillite Position Data

In the script, if you set `solveForSatPos` to `True` then the script will by default calculate the positions of CALIPSO 10 seconds before and 10 seconds after culmination in RA and Dec inside an HDF5 file. Each table is referenced by a key which is the time of the first element in [ISO-8601 format](https://en.wikipedia.org/wiki/ISO_8601). Each row inside that table has format of `<Seconds since the start>, <Satillite RA>, <Satillite Dec>`. The file will be named after the observatory the calculations were done with respect to. The [wikipedia page for HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) has a list of languages which have HDF5 parsers for them already written.

#### Dependencies
* tqdm
* pytz
* tzlocal
* numpy
* skyfield
* geopy
* h5py (If storing Satillite Positions)
