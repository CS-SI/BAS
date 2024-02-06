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
module watermask.py
: Contains classes to manipulate watermask raster
Created on Jan 18, 2023
"""
import os
import geopandas as gpd
import numpy as np
import osgeo
from osgeo import osr
import pandas as pd
from pyproj import CRS
import rasterio as rio
from rasterio.features import shapes, rasterize
from shapely.geometry import Point, Polygon, MultiPolygon, LineString, MultiLineString
from shapely.ops import linemerge

from tools import FileExtensionError

os.environ['USE_PYGEOS'] = '0'

def reduce_section(lin_long_in, pol_in):
    """Reduce linestring to the shortest linestring with the control polygon

    Parameters
    ----------
    lin_long_in : LineString
    pol_in : Polygon

    Returns
    -------
    lin_out : LineString

    """

    lin_cut = pol_in.intersection(lin_long_in)
    if isinstance(lin_cut, MultiLineString):
        l_xy = []
        for geom in lin_cut.geoms:
            l_xy += list(geom.coords)
        lin_out = LineString(l_xy)
    elif isinstance(lin_cut, LineString):
        lin_out = lin_cut
    else:
        raise NotImplementedError

    return lin_out


class WaterMask:

    def __init__(self):
        """Class constructor
        """

        self.origin = None
        self.rasterfile = None
        self.band = None

        self.bbox = None
        self.crs = None
        self.crs_epsg = None
        self.coordsyst = None  # "proj" or "lonlat"
        self.res = None

        self.bool_cleaned = False
        self.bool_labelled = False

    @classmethod
    def from_tif(cls, watermask_tif=None, str_origin=None, str_proj="proj"):
        """ Instanciate from any GeoTiff file

        Parameters
        ----------
        watermask_tif : str
            Path to GeoTiff file containting watermask
        str_origin : str
        str_proj : str
            ["proj", "lonlat"]

        """

        klass = WaterMask()

        # Set watermask origin (can be anything)
        klass.origin = str_origin

        # Set watermask rasterfile
        if not os.path.isfile(watermask_tif):
            raise FileExistsError("Input watermak_tif file does not exist..")
        else:
            if not watermask_tif.endswith(".tif"):
                raise FileExtensionError(message="Input file is not a .tif")
        klass.rasterfile = watermask_tif

        # Set raster coordinate system
        if str_proj not in ["proj", "lonlat"]:
            raise NotImplementedError("coordsyst available options are proj or lonlat")
        klass.coordsyst = str_proj

        # Set raster bounding box, crs and resolution
        with rio.open(watermask_tif, 'r') as src:
            klass.crs = src.crs
            klass.crs_epsg = src.crs.to_epsg()
            klass.res = src.transform.a
            klass.bbox = (src.bounds.left,
                          src.bounds.bottom,
                          src.bounds.right,
                          src.bounds.top)

        return klass

    def __str__(self):
        """ str method

        Returns
        -------
        message : str

        """
        if self.origin is not None:
            message = "Watermask product from {}\n".format(self.origin)
        else:
            message = "Watermask product from {}\n".format(self.rasterfile)

        return message

    def get_bbox(self):
        """Derive bounding box of current WaterMask product in lon-lat system

        Returns
        -------
        minlon : float
            minimum longitude
        minlat : float
            minimum latitude
        maxlon : float
            maximum longitude
        maxlat : float
            maximum latitude
        """

        if self.coordsyst == "lonlat":

            minlon = self.bbox[0]
            minlat = self.bbox[1]
            maxlon = self.bbox[2]
            maxlat = self.bbox[3]

        elif self.coordsyst == "proj":

            src = osr.SpatialReference()
            # src.ImportFromEPSG(self.crs.to_epsg(confidence_threshold=20))
            src.ImportFromProj4(self.crs.to_proj4())

            tgt = osr.SpatialReference()
            tgt.ImportFromEPSG(4326)
            if osgeo.__version__ >= "3.6.4":
                tgt.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

            osr_transform = osr.CoordinateTransformation(src, tgt)

            lonlat_bottomleft_edge = osr_transform.TransformPoint(self.bbox[0], self.bbox[1])
            lonlat_topleft_edge = osr_transform.TransformPoint(self.bbox[0], self.bbox[3])
            lonlat_bottomright_edge = osr_transform.TransformPoint(self.bbox[2], self.bbox[1])
            lonlat_topright_edge = osr_transform.TransformPoint(self.bbox[2], self.bbox[3])

            minlon = min([lonlat_bottomleft_edge[0], lonlat_topleft_edge[0]])
            maxlon = max([lonlat_bottomright_edge[0], lonlat_topright_edge[0]])
            minlat = min([lonlat_bottomleft_edge[1], lonlat_bottomright_edge[1]])
            maxlat = max([lonlat_topleft_edge[1], lonlat_topright_edge[1]])


        else:
            raise NotImplementedError

        return minlon, minlat, maxlon, maxlat

    def clean_watermask(self, gdf_reaches=None, out_dir=".", scn_name="scn", gdf_waterbodies=None):
        """Clean watermask from non-river waterbody

        Parameters
        ----------
        gdf_reaches : gpd.GeoDataFrame
            with shapely.geometry.LineString geometries
        scn_name : str
        out_dir : str
            path towards a directory where store cleaned watermask as a GeoTiff
        gdf_waterbodies : gpd.GeoDataFrame
            [optional] reference waterbodies : shapely.geometry.Polygon (or MultiPolygon) geometries

        Returns
        -------

        """

        print(" ---- Cleaning watermask ---- ")

        # Gather reaches and project them into the watermask coordinate system
        gser_reach_geom_proj = gdf_reaches["geometry"].to_crs(epsg=self.crs_epsg)

        with rio.open(self.rasterfile, 'r') as raster:

            # Get raw watermask band and mask of nodata
            band = raster.read(1)
            mask = raster.read_masks(1)

            # Clean watermask using waterbody database as reference (if provided)
            l_shapes_waterbodies = []
            if gdf_waterbodies is not None:

                # Check waterbodies coordinate system and reproject if necessary
                gdf_waterbodies_wrk = gdf_waterbodies.to_crs(self.crs)

                # Polygonize watermask and keep polygons that intersect waterbody database
                for shape, value in shapes(band, mask=(band == 1), transform=raster.transform):
                    l_polygons = []
                    for coords in shape["coordinates"]:
                        l_polygons.append(Polygon(coords))
                    if len(l_polygons) > 1:
                        big_pol = l_polygons[0]
                        small_pol = MultiPolygon(l_polygons[1:])
                        polygon = big_pol.difference(small_pol)
                    else:
                        polygon = l_polygons[0]
                    ser_intersects = gdf_waterbodies_wrk["geometry"].apply(
                        lambda pol: 1 if polygon.intersects(pol) else 0)
                    ser_count = ser_intersects.value_counts()
                    try:
                        if ser_count[1] > 0:
                            l_shapes_waterbodies.append(polygon)
                    except KeyError:
                        pass

            # Clean watermask using SWORD reaches
            l_shapes_reaches = []

            # Polygonize watermask and keep polygons that intersect SWORD reaches
            for shape, value in shapes(band, mask=(band == 1), transform=raster.transform):
                l_polygons = []
                for coords in shape["coordinates"]:
                    l_polygons.append(Polygon(coords))
                    if len(l_polygons) > 1:
                        big_pol = l_polygons[0]
                        small_pol = MultiPolygon(l_polygons[1:])
                        polygon = big_pol.difference(small_pol)
                    else:
                        polygon = l_polygons[0]
                    ser_intersects = gser_reach_geom_proj.apply(
                        lambda line: 1 if line.intersects(polygon) else 0)
                    npar_int_intersects = np.array(ser_intersects.tolist())
                    if np.sum(npar_int_intersects) > 0:
                        l_shapes_reaches.append(polygon)

            # Gather cleaning shapes
            l_shapes = l_shapes_reaches + l_shapes_waterbodies

            # Rasterize-back extracted shapes
            if len(l_shapes) > 0:
                band_new = rasterize(l_shapes,
                                     out_shape=band.shape,
                                     default_value=1,
                                     transform=raster.transform)
                band_new = band_new.astype(np.uint8)
                band_clean = np.where(mask == 0, raster.nodata, band_new)

            else:
                band_clean = band.copy()

            # Save clean watermask into a new GeoTiff file
            self.rasterfile = os.path.join(out_dir, scn_name + "_clean.tif")
            band_dtype = rio.uint8
            raster_nodata = 255
            with rio.open(
                    self.rasterfile,
                    mode="w",
                    driver="GTiff",
                    height=band_clean.shape[0],
                    width=band_clean.shape[1],
                    count=1,
                    dtype=band_dtype,
                    crs=self.crs,
                    transform=raster.transform,
                    nodata=raster_nodata
            ) as new_dataset:
                new_dataset.write(band_clean, 1)

            self.bool_cleaned = True
        return self.rasterfile

    def label_watermask(self, gdf_reaches=None, attr_reachid=None, out_dir=".", scn_name="scn"):
        """Label watermask into individual regions associated to a unique reach each

        Parameters
        ----------
        gdf_reaches : gpd.GeoDataFrame
            with shapely.geometry.LineString geometries
        attr_reachid : str
        scn_name : str
        out_dir : str
            path towards a directory where store labeled watermask as a GeoTiff

        Returns
        -------

        """

        print(" ---- Label watermask ---- ")

        if not self.bool_cleaned:
            raise Warning("Watermask not yet cleaned, should be done before..")

        # Gather reaches and project them into the watermask coordinate system
        gser_reach_geom_proj = gdf_reaches["geometry"].to_crs(epsg=self.crs_epsg)
        gdf_reaches_proj = gpd.GeoDataFrame(
            pd.DataFrame({"reach_id": gdf_reaches[attr_reachid].tolist()}),
            geometry=gser_reach_geom_proj,
            crs=CRS(self.crs_epsg)
        )

        with rio.open(self.rasterfile, 'r') as raster:

            # Get raw watermask band
            band_clean = raster.read(1)
            mask_clean = raster.read_masks(1)

            # Turn watermask into a point-cloud format
            band_clean_flat = band_clean.flatten()
            indices = np.where(band_clean_flat == 1)[0]
            l_coords = [raster.xy(i, j) for (i, j) in
                        zip(np.unravel_index(indices, band_clean.shape)[0],
                            np.unravel_index(indices, band_clean.shape)[1])]
            l_index = [t for t in np.unravel_index(indices, band_clean.shape)]
            gser_pixc = gpd.GeoSeries([Point(t[0], t[1]) for t in l_coords], crs=raster.crs)
            gdf_pixc = gpd.GeoDataFrame(
                pd.DataFrame({"i": [i for i in l_index[0]], "j": [j for j in l_index[1]], "indice": indices
                              },
                             index=pd.Index(range(len(l_coords))), dtype=np.int64),
                geometry=gser_pixc,
                crs=raster.crs
            )

            # Segment pixc into reaches
            gdf_label = gpd.sjoin_nearest(left_df=gdf_pixc, right_df=gdf_reaches_proj, max_distance=3000., how="inner",
                                          distance_col="dist")
            if len(gdf_reaches_proj) < 254:
                band_label_flat = 255 * np.ones(band_clean.shape, dtype=np.uint8).flatten()
                band_dtype = rio.uint8
                raster_nodata = 255
            else:
                band_label_flat = -1 * np.ones(band_clean.shape, dtype=np.int16).flatten()
                band_dtype = rio.int16
                raster_nodata = -1
            band_label_flat[gdf_label["indice"].to_numpy()] = gdf_label["index_right"].to_numpy()
            band_label = band_label_flat.reshape(band_clean.shape)

            # for width extraction later, can't have a region with value 0
            band_label = np.where(band_label!=raster_nodata, band_label+1, band_label)

            # Save clean watermask into a new GeoTiff file
            self.rasterfile = os.path.join(out_dir, scn_name + "_label.tif")
            with rio.open(
                    self.rasterfile,
                    mode="w",
                    driver="GTiff",
                    height=band_label.shape[0],
                    width=band_label.shape[1],
                    count=1,
                    dtype=band_dtype,
                    crs=self.crs,
                    transform=raster.transform,
                    nodata=raster_nodata
            ) as new_dataset:
                new_dataset.write(band_label, 1)

        self.bool_labelled = True
        return self.rasterfile

    def reduce_sections(self, gdf_reaches=None, attr_reachid=None, gdf_sections_in=None):
        """Reduce sections geometry to its associated region

        Parameters
        ----------
        gdf_reaches : gpd.GeoDataFrame
            with shapely.geometry.LineString geometries
        attr_reachid : str
        gdf_sections_in : gpd.GeoDataFrame
            input sections geometry

        Returns
        -------
        gdf_sections_out : gpd.GeoDataFrame
            reduced sections geometry

        """

        print(" ---- Reduce sections ---- ")

        # Check if section geometries are in the right coordinate system
        if gdf_sections_in.crs != self.crs:
            raise ValueError("Sections geometries and watermask raster has to be in the same coordinate system")

        with rio.open(self.rasterfile) as src:

            band = src.read(1)

            # Vectorize regions from watermask raster
            dct_regions = {}
            for shape, value in shapes(band, mask=(band > 0), transform=src.transform):

                l_polygons = []
                for coords in shape["coordinates"]:
                    l_polygons.append(Polygon(coords))
                if len(l_polygons) > 1:
                    big_pol = l_polygons[0]
                    small_pol = MultiPolygon(l_polygons[1:])
                    polygon = big_pol.difference(small_pol)
                else:
                    polygon = l_polygons[0]

                if value in dct_regions.keys():
                    dct_regions[int(value)]["all"].append(polygon)
                else:
                    dct_regions[int(value)] = {}
                    dct_regions[int(value)]["all"] = [polygon]

            for key, value in dct_regions.items():
                dct_regions[key]["region"] = MultiPolygon(dct_regions[key]["all"])

            # Reduce sections geometry to within the associated region
            if self.bool_labelled :
                l_gdfsub_sections = []
                for index in dct_regions.keys():

                    if index != src.nodata:

                        # Extract sections associated to unique current reach
                        reach_id = gdf_reaches.at[index-1, attr_reachid]
                        gdfsub_sections_byreach = gdf_sections_in[gdf_sections_in[attr_reachid] == reach_id].copy(
                            deep=True)

                        # Get current region
                        pol_region = dct_regions[index]["region"]

                        # In sections subset, keep only sections that intersect current region
                        ser_bool_intersects = gdfsub_sections_byreach["geometry"].intersects(pol_region)
                        gdfsub_sections_byreach_onregion = gdfsub_sections_byreach[ser_bool_intersects].copy(deep=True)

                        # For remaining stations/sections, reduce their geometry to within the current region
                        if len(gdfsub_sections_byreach_onregion) > 0:
                            gdfsub_sections_byreach_onregion["geometry"] = gdfsub_sections_byreach_onregion[
                                "geometry"].apply( lambda line: reduce_section(line, pol_region) )
                        l_gdfsub_sections.append(gdfsub_sections_byreach_onregion)

                # Gather all sections
                gdf_sections_wrk_out = pd.concat(l_gdfsub_sections)

            else:
                gdf_sections_wrk_out = gdf_sections_in.copy(deep=True)

                # If not label step, region is a unique polygon with value 1
                pol_region = dct_regions[1]["region"]

                # Reduce section geometry
                gdf_sections_wrk_out["geometry"] = gdf_sections_in["geometry"].apply(lambda line: reduce_section(line, pol_region))

            gdf_sections_out = gdf_sections_wrk_out.copy(deep=True)

        return gdf_sections_out
