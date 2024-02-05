# Copyright (C) 2023-2024 CS GROUP France, https://csgroup.eu
#
# This file is part of BAS (Buffer Around Sections)
#
#     https://github.com/CS-SI/BAS
#
# Authors:
#     Charlotte Emery
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
tools.py
: module containing diverse python tools to follow execution
"""
import os
import time
import pyproj

os.environ['USE_PYGEOS'] = '0'


class FileExtensionError(TypeError):

    def __init__(self, message="File has to hage a different extension file"):
        """Class constructor
        """
        self.message = message
        super().__init__(self.message)


class DisjointBboxError(ValueError):

    def __init__(self, message="Compared scenes do not overlap."):
        """Class constructor
        """
        self.message = message
        super().__init__(self.message)


class Timer:

    def __init__(self):
        """Class constructor
        """

        self.start_time = 0.0
        self.tmp_time = 0.0
        self.stop_time = 0.0

    def start(self):
        """Initialize the time
        """

        self.start_time = time.time()

    def stop(self):
        """Stop the timer
        Returns
        -------
        message : str
        To display total duration time in seconds
        """

        # Calculate time since Timer initialization
        self.stop_time = time.time() - self.start_time

        # Output message
        message = "Total execution time in %s" % self.stop_time

        return message

    def info(self, flag_start=1):
        """Current time

        Parameters
        ----------
        flag_start : int
        set whether temporary duration is from the Timer initialiation (=1) or not

        Returns
        -------
        message : str
        To display total temporary time in seconds

        """

        # calulate duration time
        if (self.tmp_time == 0) or (flag_start == 1):
            current_time = time.time() - self.start_time
        else:
            current_time = time.time() - self.tmp_time

        # Update temporary time
        self.tmp_time = time.time()

        # Output message
        message = "Step executed in %s" % current_time

        return message


class RDFParser:

    def __init__(self, rdf_paramfile):
        """Class constructor

        Parameters
        ----------
        rdf_paramfile : str
        path to rdf file to parse
        """

        self.rdf_paramfile = rdf_paramfile
        self.rdf_dict = {}  # Dictionnary to store info from rdf file

        self.read_rdffile()

    def read_rdffile(self):
        """Fill rdf dictionnary
        """

        with open(self.rdf_paramfile, "r") as rdfparser:

            rdflines = rdfparser.readlines()

            for line in rdflines:
                if "=" in line:
                    line_key = line.split("=")[0].strip()

                    try:
                        line_value = line.split("=")[1].strip()
                        if line_value == "None":
                            line_value = None
                    except IndexError:
                        line_value = None

                    self.rdf_dict[line_key] = line_value

    @property
    def config(self):
        return self.rdf_dict


def project(lon, lat, proj="laea", ellps="WGS84", lon_0=None, lat_0=None, x_0=0.0, y_0=0.0, inverse=False, **proj_kwds):
    """Project/Reproject from WGS84-lon-lat into local "laea"

    Parameters
    ----------
    proj : str
        name of the output (input if inverse is True) projection
    ellps : str
        name eof the input (output if inverse is True) ellipsoid
    lon_0 :
    lat_0 :
    x_0 :
    y_0 :
    inverse : bool
        inverse transformation
    proj_kwds :
        other arguments for the projection function

    Returns
    -------

    """

    fct = pyproj.Proj(proj=proj,
                      lat_0=lat_0,
                      lon_0=lon_0,
                      x_0=x_0,
                      y_0=y_0,
                      ellps=ellps,
                      **proj_kwds)

    if not inverse:
        x, y = fct(lon, lat)
    else:
        x, y = fct(lon, lat, inverse=True)

    return x, y
