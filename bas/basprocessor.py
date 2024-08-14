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
import numpy as np

from tools import DisjointBboxError, FileExtensionError
from watermask import WaterMask
from widths import compute_widths_from_single_watermask


# os.environ['USE_PYGEOS'] = '0'


class BASProcessor:

    def __init__(self,
                 str_watermask_tif=None,
                 gdf_sections=None,
                 gdf_reaches=None,
                 attr_reachid=None,
                 str_proj="proj",
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
        self.dct_cfg = {"clean": {"bool_toclean": True,
                                  "type_clean": "base",
                                  "fpath_wrkdir": ".",
                                  "gdf_waterbodies": None
                                  },
                        "label": {"bool_tolabel": False,
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
                      "type_clean" : base/waterbodies,
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
            self.clean_watermask(dct_cfg)
        self.watermask.save_wm(fmt="tif",
                               bool_clean=True,
                               bool_label=False,
                               str_fpath_dir_out="/home/cemery/Work/git/BAS/examples",
                               str_suffix="debug_example3")
        self.watermask.save_wm(fmt="pixc",
                               bool_clean=True,
                               bool_label=False,
                               str_fpath_dir_out="/home/cemery/Work/git/BAS/examples",
                               str_suffix="debug_example3")
        self.watermask.save_wm(fmt="shp",
                               bool_clean=True,
                               bool_label=False,
                               str_fpath_dir_out="/home/cemery/Work/git/BAS/examples",
                               str_suffix="debug_example3")

        # # Label watermask
        # if dct_cfg["label"]["bool_label"]:
        #     self.label_watermask(gdf_reaches=self.gdf_reaches,
        #                                    attr_reachid=self.attr_reachid,
        #                                    out_dir=dct_cfg["label"]["fpath_wrkdir"],
        #                                    scn_name=self.scene_name)

        # Prepare sections
        # gdf_wrk_sections = self.watermask.reduce_sections(gdf_reaches=self.gdf_reaches,
        #                                                   attr_reachid=self.attr_reachid,
        #                                                   gdf_sections_in=self.gdf_sections)
        gdf_wrk_sections = self.gdf_sections

        # Process width
        print("---- Compute widths ----")
        with rio.open(self.watermask.str_fpath_infile) as src:
            gdf_widths, _ = compute_widths_from_single_watermask(scenario=dct_cfg["widths"]["scenario"],
                                                                 watermask=src,
                                                                 sections=gdf_wrk_sections,
                                                                 buffer_length=8.0 * self.watermask.res)

        print("")
        print("----- Processing : Done -----")

        return gdf_widths

    def clean_watermask(self, dct_cfg=None):
        """Clean watermask from non-river waterbodies

        Parameters
        ----------
        Returns
        -------

        """

        print(" ---- Cleaning watermask ---- ")

         # Check config_dct
        if dct_cfg is None:
            dct_cfg=self.dct_cfg

        # Gather reaches and project them into the watermask coordinate system
        gdf_reaches_proj = self.gdf_reaches.to_crs(epsg=self.watermask.crs_epsg)

        # Gather wm as polygons
        gdf_wm_polygons = self.watermask.get_polygons(bool_clean=False,
                                                      bool_label=False,
                                                      bool_exterior_only=False,
                                                      bool_indices=True)

        # Apply regular cleaning
        gdf_join_wm_reaches = gpd.sjoin(left_df=gdf_wm_polygons,
                                        right_df=gdf_reaches_proj,
                                        how="inner",
                                        predicate="intersects")
        npar_idx_pol_notclean = np.setdiff1d(
            np.unique(gdf_wm_polygons.index),
            np.unique(gdf_join_wm_reaches.index)
        )

        # Apply waterbodies-type cleaning if activated
        if dct_cfg["clean"]["gdf_waterbodies"] is not None and dct_cfg["clean"]["type_clean"] == "waterbodies":

            if not isinstance(dct_cfg["clean"]["gdf_waterbodies"], gpd.GeoDataFrame):
                raise TypeError
            else:
                gdf_waterbodies_wrk = dct_cfg["clean"]["gdf_waterbodies"].to_crs(gdf_wm_polygons.crs)

            gdf_join_wm_waterbodies = gpd.sjoin(left_df=gdf_wm_polygons,
                                                right_df=gdf_waterbodies_wrk,
                                                how="inner",
                                                predicate="intersects")

            npar_idx_notclean_wb = np.setdiff1d(
                np.unique(gdf_wm_polygons.index),
                np.unique(gdf_join_wm_waterbodies.index)
            )

            npar_idx_pol_notclean = np.intersect1d(npar_idx_pol_notclean, npar_idx_notclean_wb)

        gdfsub_notclean_wm_polygons = gdf_wm_polygons.loc[npar_idx_pol_notclean,:].copy()
        gdfsub_notclean_wm_polygons.loc[:,["clean", "label", "geometry"]].to_file("/home/cemery/Work/git/BAS/examples/wm_notclean.shp")

        # Get pixc indexes from polygon indexes
        l_idx_pixc_notclean = [
            element for list_ in gdfsub_notclean_wm_polygons["indices"].values for element in list_
        ]

        # Update clean flag in the watermask
        self.watermask.update_clean_flag(mask=l_idx_pixc_notclean)


    # def label_watermask(self, gdf_reaches=None, attr_reachid=None, out_dir=".", scn_name="scn"):
    #     """Label watermask into individual regions associated to a unique reach each
    #
    #     Parameters
    #     ----------
    #     gdf_reaches : gpd.GeoDataFrame
    #         with shapely.geometry.LineString geometries
    #     attr_reachid : str
    #     scn_name : str
    #     out_dir : str
    #         path towards a directory where store labeled watermask as a GeoTiff
    #
    #     Returns
    #     -------
    #
    #     """
    #
    #     print(" ---- Label watermask ---- ")
    #
    #     if not self.bool_cleaned:
    #         raise Warning("Watermask not yet cleaned, should be done before..")
    #
    #     # Gather reaches and project them into the watermask coordinate system
    #     gser_reach_geom_proj = gdf_reaches["geometry"].to_crs(epsg=self.crs_epsg)
    #     gdf_reaches_proj = gpd.GeoDataFrame(
    #         pd.DataFrame({"reach_id": gdf_reaches[attr_reachid].tolist()}),
    #         geometry=gser_reach_geom_proj,
    #         crs=CRS(self.crs_epsg)
    #     )
    #
    #     with rio.open(self.rasterfile, 'r') as raster:
    #
    #         # Get raw watermask band
    #         band_clean = raster.read(1)
    #         mask_clean = raster.read_masks(1)
    #
    #         # Turn watermask into a point-cloud format
    #         band_clean_flat = band_clean.flatten()
    #         indices = np.where(band_clean_flat == 1)[0]
    #         l_coords = [raster.xy(i, j) for (i, j) in
    #                     zip(np.unravel_index(indices, band_clean.shape)[0],
    #                         np.unravel_index(indices, band_clean.shape)[1])]
    #         l_index = [t for t in np.unravel_index(indices, band_clean.shape)]
    #         gser_pixc = gpd.GeoSeries([Point(t[0], t[1]) for t in l_coords], crs=raster.crs)
    #         gdf_pixc = gpd.GeoDataFrame(
    #             pd.DataFrame({"i": [i for i in l_index[0]], "j": [j for j in l_index[1]], "indice": indices
    #                           },
    #                          index=pd.Index(range(len(l_coords))), dtype=np.int64),
    #             geometry=gser_pixc,
    #             crs=raster.crs
    #         )
    #
    #         # Segment pixc into reaches
    #         gdf_label = gpd.sjoin_nearest(left_df=gdf_pixc, right_df=gdf_reaches_proj, max_distance=3000., how="inner",
    #                                       distance_col="dist")
    #         if len(gdf_reaches_proj) < 254:
    #             band_label_flat = 255 * np.ones(band_clean.shape, dtype=np.uint8).flatten()
    #             band_dtype = rio.uint8
    #             raster_nodata = 255
    #         else:
    #             band_label_flat = -1 * np.ones(band_clean.shape, dtype=np.int16).flatten()
    #             band_dtype = rio.int16
    #             raster_nodata = -1
    #         band_label_flat[gdf_label["indice"].to_numpy()] = gdf_label["index_right"].to_numpy()
    #         band_label = band_label_flat.reshape(band_clean.shape)
    #
    #         # for width extraction later, can't have a region with value 0
    #         band_label = np.where(band_label!=raster_nodata, band_label+1, band_label)
    #
    #         # Save clean watermask into a new GeoTiff file
    #         self.rasterfile = os.path.join(out_dir, scn_name + "_label.tif")
    #         with rio.open(
    #                 self.rasterfile,
    #                 mode="w",
    #                 driver="GTiff",
    #                 height=band_label.shape[0],
    #                 width=band_label.shape[1],
    #                 count=1,
    #                 dtype=band_dtype,
    #                 crs=self.crs,
    #                 transform=raster.transform,
    #                 nodata=raster_nodata
    #         ) as new_dataset:
    #             new_dataset.write(band_label, 1)
    #
    #     self.bool_labelled = True
    #     return self.rasterfile
    #
    # def reduce_sections(self, gdf_reaches=None, attr_reachid=None, gdf_sections_in=None):
    #     """Reduce sections geometry to its associated region
    #
    #     Parameters
    #     ----------
    #     gdf_reaches : gpd.GeoDataFrame
    #         with shapely.geometry.LineString geometries
    #     attr_reachid : str
    #     gdf_sections_in : gpd.GeoDataFrame
    #         input sections geometry
    #
    #     Returns
    #     -------
    #     gdf_sections_out : gpd.GeoDataFrame
    #         reduced sections geometry
    #
    #     """
    #
    #     print(" ---- Reduce sections ---- ")
    #
    #     # Check if section geometries are in the right coordinate system
    #     if gdf_sections_in.crs != self.crs:
    #         raise ValueError("Sections geometries and watermask raster has to be in the same coordinate system")
    #
    #     with rio.open(self.rasterfile) as src:
    #
    #         band = src.read(1)
    #
    #         # Vectorize regions from watermask raster
    #         dct_regions = {}
    #         for shape, value in shapes(band, mask=(band > 0), transform=src.transform):
    #
    #             l_polygons = []
    #             for coords in shape["coordinates"]:
    #                 l_polygons.append(Polygon(coords))
    #             if len(l_polygons) > 1:
    #                 big_pol = l_polygons[0]
    #                 small_pol = MultiPolygon(l_polygons[1:])
    #                 polygon = big_pol.difference(small_pol)
    #             else:
    #                 polygon = l_polygons[0]
    #
    #             if value in dct_regions.keys():
    #                 dct_regions[int(value)]["all"].append(polygon)
    #             else:
    #                 dct_regions[int(value)] = {}
    #                 dct_regions[int(value)]["all"] = [polygon]
    #
    #         for key, value in dct_regions.items():
    #             dct_regions[key]["region"] = MultiPolygon(dct_regions[key]["all"])
    #
    #         # Reduce sections geometry to within the associated region
    #         if self.bool_labelled :
    #             l_gdfsub_sections = []
    #             for index in dct_regions.keys():
    #
    #                 if index != src.nodata:
    #
    #                     # Extract sections associated to unique current reach
    #                     reach_id = gdf_reaches.at[index-1, attr_reachid]
    #                     gdfsub_sections_byreach = gdf_sections_in[gdf_sections_in[attr_reachid] == reach_id].copy(
    #                         deep=True)
    #
    #                     # Get current region
    #                     pol_region = dct_regions[index]["region"]
    #
    #                     # In sections subset, keep only sections that intersect current region
    #                     ser_bool_intersects = gdfsub_sections_byreach["geometry"].intersects(pol_region)
    #                     gdfsub_sections_byreach_onregion = gdfsub_sections_byreach[ser_bool_intersects].copy(deep=True)
    #
    #                     # For remaining stations/sections, reduce their geometry to within the current region
    #                     if len(gdfsub_sections_byreach_onregion) > 0:
    #                         gdfsub_sections_byreach_onregion["geometry"] = gdfsub_sections_byreach_onregion[
    #                             "geometry"].apply( lambda line: reduce_section(line, pol_region) )
    #                     l_gdfsub_sections.append(gdfsub_sections_byreach_onregion)
    #
    #             # Gather all sections
    #             gdf_sections_wrk_out = pd.concat(l_gdfsub_sections)
    #
    #         else:
    #             gdf_sections_wrk_out = gdf_sections_in.copy(deep=True)
    #
    #             # If not label step, region is a unique polygon with value 1
    #             pol_region = dct_regions[1]["region"]
    #
    #             # Reduce section geometry
    #             gdf_sections_wrk_out["geometry"] = gdf_sections_in["geometry"].apply(lambda line: reduce_section(line, pol_region))
    #
    #         gdf_sections_out = gdf_sections_wrk_out.copy(deep=True)
    #
    #     return gdf_sections_out
