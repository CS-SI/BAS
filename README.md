# BAS
[ Buffer Around Sections]

> Derive river width from a binary watermask and a set of sections

BAS is a simple Python-based toolbox designed to compute river widths from a space-borne watermask given a set of section (line geometries) set along the river.

## Download

Main version available on github : https://github.com/CS-SI/BAS

## Installation

### With conda and a yml file

The toolbox works in a specific environment that can be installed through conda using the environment_bas.yml file

```shell
conda env create -f environment_bas_py38.yml # for python 3.8
# or
conda env create -f environment_bas_py310.yml # for python 3.10

# then install BAS (make sure you are in the source directory)
pip install -e .

# then try to access the entry point
run_bas -h
```

### General

Using a virtual environment is strongly recommended.

#### Locally

If you have python-3.12 and `sw1Dto2D` is already installed, then running this command will install BAS requirements

```bash
pip install -e .

# then try to access the entry point
run_bas -h
```

This entrypoint is to compute BAS widths on Surfwater-like watermasks.

#### On TREX cluster

```bash
# On TREX cluster load a working version of python 3.12
module load otb/9.0.0-python3.12

# then create a dedicated virtual env
python3 -m venv ${YOUR_VENV_PATH}

# you will need sw1dto2d for python 3.12
cd ${SW1DTO2D_DIR_PATH}
${YOUR_VENV_PATH}/bin/pip install -e .

# then go back to your source directory to install BAS module
cd ${BAS_DIR_PATH}
${YOUR_VENV_PATH}/bin/pip install -e .

# then try to access the entry point
source ${YOUR_VENV_PATH}/bin/activate
run_bas -h
```

This entrypoint is to compute BAS widths on Surfwater-like watermasks.

## Examples

To run examples, use the entrypoint

```bash
run_examples
```

## Usage

To use BAS on Surfwater-like watermasks, here is the general command

```bash
run_bas -w /path/to/water_mask.tif -dt YYYYmmddThhmmss -r /path/to/eu_sword_reaches_hb23_v16.shp -n /path/to/eu_sword_nodes_hb23_v16.shp -o /path/to/output/directory
```

This will only save a csv file with widths data. You can had `--more_outputs` to also save the cleaned watermask and the shapefile.  


## Features

Using the BAS toolbox, you can perform two tasks :

- Given a watermask (as a GeoTiff file), a set on centerline reaches (as a shapefile of LineString) and a set of sections (as a shapefile of LineString),
derive the river widths along the sections
- If the sections lines are not available, using the centerline reaches (as a shapefile of LineString) and a set of segmentation points (as a shapefile of Point),
you can draw the sections yourself


## License

BAS is licensed by [CS GROUP](https://www.c-s.fr/) under
the [Apache License, version 2.0](http://www.apache.org/licenses/LICENSE-2.0.html).
A copy of this license is provided in the [LICENSE](LICENSE) file.