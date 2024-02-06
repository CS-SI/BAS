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
basprocessor.py
: From an external watermask + a (set of) river centerline(s) and associated node along it (stations or calculus points),
derive a width estimation at said-nodes observed within the mask
"""
import os
from datetime import datetime
import geopandas as gpd
from pyproj import CRS
import rasterio as rio
import shapely

from tools import DisjointBboxError, FileExtensionError
from watermask import WaterMask
from widths import compute_widths_from_single_watermask

os.environ['USE_PYGEOS'] = '0'


class BASProcessor:

    def __init__(self, str_watermask_tif=None, gdf_sections=None, gdf_reaches=None, attr_reachid=None, str_proj="proj",
                 str_provider=None,
                 str_datetime=None):
        """Class constructor
        """

        # Set watermask information
        if str_watermask_tif is not None:
            self.f_watermask_in = str_watermask_tif
            if not os.path.isfile(str_watermask_tif):
                raise FileExistsError("Input watermask GeoTiff does not exist")
        else:
            raise ValueError("Missing watermask GeoTiff input file")
        self.scene_name = os.path.basename(self.f_watermask_in).split(".")[0]
        self.proj = str_proj
        self.provider = str_provider
        self.watermask = None

        # Set sections information
        if gdf_sections is None:
            raise ValueError("Missing reaches geometries")
        else:
            if not isinstance(gdf_sections, gpd.GeoDataFrame):
                raise TypeError
        self.gdf_sections = gdf_sections

        # Set reaches information
        if gdf_reaches is None:
            raise ValueError("Missing reaches geometries")
        else:
            if not isinstance(gdf_reaches, gpd.GeoDataFrame):
                raise TypeError
        self.gdf_reaches = gdf_reaches
        self.attr_reachid = attr_reachid

        # Set datetime information if provided
        if str_datetime is not None:
            self.scene_datetime = str_datetime
            try:
                dt_scene = datetime.strptime(self.scene_datetime, "%Y%m%dT%H%M%S")
            except ValueError:
                raise ValueError("input datetime {} does not match format '%Y%m%dT%H%M%S'.".format(self.scene_datetime))
        else:
            self.scene_datetime = None

        # Processing default values
        self.dct_cfg = {"clean": {"bool_clean": True,
                                  "type_label": "base",
                                  "fpath_wrkdir": ".",
                                  "gdf_waterbodies": None
                                  },
                        "label": {"bool_label": False,
                                  "type_label": "base",
                                  "fpath_wrkdir": "."
                                  },
                        "widths": {"scenario": 0
                                   }
                        }

    def check_bbox_compatibility(self):
        """Check if sections and watermask are spatially compatible
        """

        print("----- WidthProcessing = CheckBboxCompatibility -----")
        print("")

        if self.gdf_sections.crs != CRS(4326):
            gdf_crs_sections = self.gdf_sections.to_crs(CRS(4326))
        else:
            gdf_crs_sections = self.gdf_sections.copy(deep=True)
        bbox_sections = gdf_crs_sections.total_bounds
        polygon_sections = shapely.geometry.box(*bbox_sections)  # Convert bbox to Polygon object
        # Shapely 2.0 : box(xmin, ymin, xmax, ymax, ccw=True, **kwargs)

        bbox_watermask = self.watermask.get_bbox()  # Load watermask bbox in lonlat
        polygon_watermask = shapely.geometry.box(*bbox_watermask)  # Convert watermask bbox to Polygon object

        if polygon_watermask.disjoint(polygon_sections):
            raise DisjointBboxError

        print("")
        print("----- CheckBboxCompatibility : Done -----")

    def preprocessing(self):
        """Preprocessing: load watermask, reproject sections et check bounding boxes intersections
        """

        print("----- WidthProcessing = Preprocessing -----")
        print("")

        # Load WaterMask object
        self.watermask = WaterMask.from_tif(self.f_watermask_in, self.provider, self.proj)

        # Reproject sections to watermask coordinate system
        self.gdf_reaches = self.gdf_reaches.to_crs(self.watermask.crs_epsg)
        self.gdf_sections = self.gdf_sections.to_crs(self.watermask.crs_epsg)

        # Check boundingbox compatibility
        self.check_bbox_compatibility()

        print("")
        print("----- Preprocessing : Done -----")

    def read_cfg(self, dct_cfg=None):
        """ Add default value to dct_cfg if keywords are missing

        Parameters
        ----------
        dct_cfg : dict

        Returns
        -------
        dct_cfg : dict

        """

        for key in ["clean", "label", "widths"]:
            if key not in dct_cfg.keys():
                dct_cfg[key] = self.dct_cfg[key]

            else:
                for subkey in self.dct_cfg[key].keys():
                    if subkey not in self.dct_cfg[key].keys():
                        dct_cfg[key][subkey] = self.dct_cfg[key][subkey]
        return dct_cfg

    def processing(self, dct_cfg=None):
        """Processing : extraction of widths from watermask

        Parameters
        ----------
        dct_cfg : set processing configuration
        { "clean" : { "bool_clean" : True/False,
                      "type_label" : base/waterbodies,
                      "fpath_wrkdir" : "."
                      "gdf_waterbodies" : gdf with polygon waterbodies to clean waterbodies [optionnal]
                    },
          "label" : { "bool_label" : True/False,
                      "type_label" : base,
                      "fpath_wrkdir" : "."
                    },
          "widths" : { scenario : 0/1/10/11
                     }
        }

        Returns
        -------
        gdf_widths : gpd.GeoDataFrame

        """

        print("----- WidthProcessing = Processing -----")
        print("")

        # Check cfg
        dct_cfg = self.read_cfg(dct_cfg)

        # Clean watermask
        if dct_cfg["clean"]["bool_clean"]:
            self.watermask.clean_watermask(gdf_reaches=self.gdf_reaches,
                                           out_dir=dct_cfg["clean"]["fpath_wrkdir"],
                                           scn_name=self.scene_name,
                                           gdf_waterbodies=dct_cfg["clean"]["gdf_waterbodies"])

        # Label watermask
        if dct_cfg["label"]["bool_label"]:
            self.watermask.label_watermask(gdf_reaches=self.gdf_reaches,
                                           attr_reachid=self.attr_reachid,
                                           out_dir=dct_cfg["label"]["fpath_wrkdir"],
                                           scn_name=self.scene_name)

        # Prepare sections
        gdf_wrk_sections = self.watermask.reduce_sections(gdf_reaches=self.gdf_reaches,
                                                          attr_reachid=self.attr_reachid,
                                                          gdf_sections_in=self.gdf_sections)

        # Process width
        print("---- Compute widths ----")
        with rio.open(self.watermask.rasterfile) as src:
            gdf_widths, _ = compute_widths_from_single_watermask(scenario=dct_cfg["widths"]["scenario"],
                                                                 watermask=src,
                                                                 sections=gdf_wrk_sections,
                                                                 buffer_length=8.0 * self.watermask.res)

        print("")
        print("----- Processing : Done -----")

        return gdf_widths
