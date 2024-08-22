# BAS
[ Buffer Around Sections]

> Derive river width from a binary watermask and a set of sections

BAS is a simple Python-based toolbox designed to compute river widths from a space-borne watermask given a set of section (line geometries) set along the river.

## Download

Main version available on github : https://github.com/CS-SI/BAS

## Getting started

The toolbox works in a specific environment that can be installed through conda using the environment_bas.yml file

```shell
conda env create -f environment_bas_py38.yml
```

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