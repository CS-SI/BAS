# Copyright (C) 2023-2024 CS GROUP France, https://csgroup.eu
# Copyright (C) 2024 CNES.
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
module rivergeomproduct.py
: Contains classes to manipulate river-related 1D-geometries
"""

import geopandas as gpd
import numpy as np
import os
import pandas as pd
from pyproj import CRS
from shapely.geometry import Point, LineString
import logging

from sw1dto2d.sw1dto2d import SW1Dto2D
from bas.tools import FileExtensionError
from bas.tools import project

_logger = logging.getLogger("bas.rivergeomproduct")


def default_reproject_reach(lin_reach, minlon, minlat, maxlon, maxlat):
    """Default reprojection of each reach geometry in laea system
    """

    arr_centerline_lon = [t[0] for t in lin_reach.coords]
    arr_centerline_lat = [t[1] for t in lin_reach.coords]
    arr_centerline_x, arr_centerline_y = project(arr_centerline_lon,
                                                 arr_centerline_lat,
                                                 lon_0=0.5 * (minlon + maxlon),
                                                 lat_0=0.5 * (minlat + maxlat))
    lin_laea_reach = LineString(
        [(x, y) for (x, y) in
         zip(arr_centerline_x, arr_centerline_y)])

    return lin_laea_reach


def get_linedge_pointwise_norm(npar_xycoord):
    """Calculate distance between each consecutive points

    Parameters
    ----------
    npar_xycoord : np.ndarray
        shape (nx, 2) with nx: number of edges

    Returns
    -------
    npar_norm : np.ndarray
        shape (nx-1,) with nx: number of edges

    """

    npar_dx = npar_xycoord[2:, 0] - npar_xycoord[0:-2, 0]
    npar_dy = npar_xycoord[2:, 1] - npar_xycoord[0:-2, 1]
    npar_dx = np.insert(npar_dx, 0, npar_xycoord[1, 0] - npar_xycoord[0, 0])
    npar_dx = np.append(npar_dx, npar_xycoord[-1, 0] - npar_xycoord[-2, 0])
    npar_dy = np.insert(npar_dy, 0, npar_xycoord[1, 1] - npar_xycoord[0, 1])
    npar_dy = np.append(npar_dy, npar_xycoord[-1, 1] - npar_xycoord[-2, 1])

    npar_norm = np.sqrt(npar_dx ** 2 + npar_dy ** 2)

    return npar_norm


def check_centerline_geometry(lin_centerline_in):
    """Remove double edges from the centerline geometry

    Parameters
    ----------
    lin_centerline_in : LineString

    Returns
    -------
    lin_centerline_out : LineString

    """

    # Get coordinates of each point within the line geometry
    npar_xycoord = np.array(lin_centerline_in.coords)

    # Get distance between consecutive edge
    npar_norm = get_linedge_pointwise_norm(npar_xycoord)

    # Check if two consecutive edge are identical or really close
    npar_idx_norm_chck = np.where(npar_norm < 0.000001)[0]

    if npar_idx_norm_chck.size == 0:  # Each point of the input line is unique
        lin_centerline_out = lin_centerline_in

    else:  # Some consecutive points are identical
        npar_base_xycoord = npar_xycoord.copy()
        npar_new_xycoord = npar_xycoord.copy()

        # While there are still non-unique points remaining in the line geometry
        while npar_idx_norm_chck.size > 0:
            # Point indexes to remove
            npar_idx_norm_chck = np.where(npar_idx_norm_chck == 0, 1, npar_idx_norm_chck)
            npar_idx_norm_chck = np.where(npar_idx_norm_chck == npar_idx_norm_chck.size - 1,
                                          npar_idx_norm_chck.size - 2, npar_idx_norm_chck)
            npar_idx_norm_chck = np.unique(npar_idx_norm_chck)

            # All point indexes
            npar_idx_norm_full = np.linspace(start=0,
                                             stop=npar_norm.size - 1,
                                             num=npar_norm.size,
                                             dtype=np.int32)

            # Remove redundant indexes
            npar_idx_norm_keep = np.setdiff1d(ar1=npar_idx_norm_full,
                                              ar2=npar_idx_norm_chck)
            npar_new_xycoord = npar_base_xycoord[npar_idx_norm_keep, :].copy()

            # Get distance between updated consecutive edge
            npar_norm = get_linedge_pointwise_norm(npar_new_xycoord)
            npar_idx_norm_chck = np.where(npar_norm < 0.000001)[0]
            npar_base_xycoord = npar_new_xycoord.copy()

        lin_centerline_out = LineString(npar_new_xycoord)

    return lin_centerline_out


def modified_compute_xs_parameters(obj_sw1dto2d, lin_centerline_in):
    """ Recalculate xs_normals attribute along a straight line

    Parameters
    ----------
    obj_sw1dto2d : SW1Dto2D object
    lin_centerline_in : LineString

    Returns
    -------
    npar_xs_normals_modified : np.array
    npar_angles : np.array

    """

    # Simplify centerline as a strainght line between both edges
    l_t_edge_s = list(lin_centerline_in.coords)[0]
    l_t_edge_e = list(lin_centerline_in.coords)[-1]

    # Get normal direction to straight line between edge of input line
    npar_dx = l_t_edge_e[0] - l_t_edge_s[0]
    npar_dy = l_t_edge_e[-1] - l_t_edge_s[-1]
    npar_norm = np.sqrt(npar_dx ** 2. + npar_dy ** 2.)
    npar_cl_nx = npar_dy / npar_norm
    npar_cl_ny = -npar_dx / npar_norm

    # Get original normals from full centerline
    npar_xs_normals_base = obj_sw1dto2d.xs_normals
    npar_xs_normals_modified = npar_xs_normals_base.copy()

    # Modify normals as all equal to the modified normal directions
    npar_xs_normals_modified[:, 0] = npar_cl_nx * np.ones(npar_xs_normals_base[:, 0].shape)
    npar_xs_normals_modified[:, -1] = npar_cl_ny * np.ones(npar_xs_normals_base[:, -1].shape)

    # calculate angle between ortho sections (original) and new sections with modified orientations
    npar_angles_base = np.arctan2(npar_xs_normals_base[:, 0], npar_xs_normals_base[:, 1])
    npar_angles = np.arctan2(npar_cl_nx, npar_cl_ny) * np.ones(npar_angles_base.shape) - npar_angles_base

    return npar_xs_normals_modified, npar_angles


class CloneSW1Dto2D(SW1Dto2D):
    """SW1Dto2D-derived class to access some hidden attributess
    """

    def __init__(self, model_output_1d: pd.DataFrame, curvilinear_abscissa_key: str, heights_key: str, widths_key: str,
                 centerline: LineString):
        """Class constructor : same as Parent class

        Parameters
        ----------
        model_output_1d :
        curvilinear_abscissa_key :
        heights_key :
        widths_key :
        centerline :
        """

        super().__init__(model_output_1d, curvilinear_abscissa_key, heights_key, widths_key, centerline)

    @property
    def xs_normals(self):
        """Get xs normals
        """
        return self._curv_absc_normals

    @xs_normals.setter
    def xs_normals(self, xs_normals_value):
        self._curv_absc_normals = xs_normals_value


class RiverGeomProduct:
    """Class to handle 1d-riverline geometries and sections
    """

    def __init__(self):
        """Class constructor
        """

        self.bool_isempty = True
        self.bool_edge = True

        # Optionnal projected system information
        self.bool_input_crs = False
        self.crs_in = None

        # 1D geometries bounding box - geographic system
        self.flt_minlon = None
        self.flt_minlat = None
        self.flt_maxlon = None
        self.flt_maxlat = None

        # 1D geometries counts
        self.int_node_dim = None  # Number of nodes within the product
        self.int_reach_dim = None  # Number of reaches within the product

        # Reach geometries information
        self.npar_int_reachgrp_reachid = None  # Reach ID
        self.dct_pcenterline = {}  # Reach LineString prior geometry defining

        # Node geometries information
        self.npar_int_nodegrp_nodeid = None  # Node id sequence
        self.npar_int_nodegrp_reachid = None  # Reach ID sequence for mapping
        self.npar_flt_nodegrp_plon = None  # Node prior lon
        self.npar_flt_nodegrp_plat = None  # Node prior lat
        self.npar_flt_nodegrp_pwidth = None  # Prior node width
        self.npar_flt_nodegrp_pwse = None  # Prior node wse

        # To produce
        self.npar_flt_nodegrp_px = None  # Node prior longitude projected
        self.npar_flt_nodegrp_py = None  # Node prior latitude projected
        self.dct_centerline = {}  # Reach-scale geometry/xs/p_wse/p_width -- use to derive sections from sw1dto2d

    @classmethod
    def from_gdf(cls, gdf_reaches=None, gdf_nodes=None, bool_edge=True, dct_attr=None, crs_in=None):
        """Instanciate object from geodataframe objects

        Parameters
        ----------
        gdf_reaches : gpd.GeoDataFrame
            Set of reaches LineString geometries
        gdf_nodes : gpd.GeoDataFrame
            Set of nodes Point geometries
        bool_edge : boolean
            Indicates if nodes segmentation includes reache edges (True) or not (False)
        dct_attr : dct
            { "reaches": { "reaches_id" : ""}, "nodes": {"reaches_id" : "", "nodes_id": "", "pwidth": "", "pwse": ""} }
            Dictionary matching input reaches/nodes attributes
        crs_in : CRS-like
            Projected CRS code in which reproject geometries to compute node xs along reaches
            If not provided, used default laea projection centered on reaches/nodes extent
        Returns
        -------
        klass : RiverGeomProduct object
        """

        klass = RiverGeomProduct()

        # Set projection system
        if crs_in is not None:
            klass.bool_input_crs = True

        # Count available geometries
        klass.int_reach_dim = len(gdf_reaches)
        klass.int_node_dim = len(gdf_nodes)
        if klass.int_reach_dim > 0 and klass.int_node_dim > 0:
            klass.isempty = False
        klass.bool_edge = bool_edge

        # Get reaches information
        klass.npar_int_reachgrp_reachid = gdf_reaches[dct_attr["reaches"]["reaches_id"]].to_numpy()
        klass.flt_minlon, klass.flt_minlat, klass.flt_maxlon, klass.flt_maxlat = gdf_reaches.total_bounds

        # Compute reach projected-geometry
        if klass.bool_input_crs:
            try:
                gser_projected_reaches = gdf_reaches["geometry"].to_crs(crs_in)
            except Exception as err:
                _logger.info("Error while trying to reproject reaches in input crs")
                _logger.error(err)
                _logger.info("Reproject in default laea system")
                gser_projected_reaches = gdf_reaches["geometry"].apply(
                    lambda lin_reach: default_reproject_reach(lin_reach,
                                                              klass.flt_minlon,
                                                              klass.flt_minlat,
                                                              klass.flt_maxlon,
                                                              klass.flt_maxlat))
        else:
            gser_projected_reaches = gdf_reaches["geometry"].apply(
                lambda lin_reach: default_reproject_reach(lin_reach,
                                                          klass.flt_minlon,
                                                          klass.flt_minlat,
                                                          klass.flt_maxlon,
                                                          klass.flt_maxlat))

        # Sort reach information
        for index, row in gdf_reaches.iterrows():
            # Set reach_id
            reach_id = row[dct_attr["reaches"]["reaches_id"]]

            # Initialize dictionnary element
            klass.dct_pcenterline[reach_id] = {}
            klass.dct_pcenterline[reach_id]["lonlat"] = row["geometry"]

            # Sort projected geometry
            klass.dct_pcenterline[reach_id]["xy"] = gser_projected_reaches.loc[index]

        # Get nodes information
        klass.npar_int_nodegrp_nodeid = gdf_nodes[dct_attr["nodes"]["nodes_id"]].to_numpy()
        klass.npar_int_nodegrp_reachid = gdf_nodes[dct_attr["nodes"]["reaches_id"]].to_numpy()

        klass.npar_flt_nodegrp_plon = gdf_nodes["geometry"].x.to_numpy()
        klass.npar_flt_nodegrp_plat = gdf_nodes["geometry"].y.to_numpy()

        try:
            klass.npar_flt_nodegrp_pwidth = gdf_nodes[dct_attr["nodes"]["pwidth"]].to_numpy()
        except KeyError:
            klass.npar_flt_nodegrp_pwidth = 600. * np.ones_like(klass.npar_flt_nodegrp_plon)
        try:
            klass.npar_flt_nodegrp_pwse = gdf_nodes[dct_attr["nodes"]["pwse"]].to_numpy()
        except KeyError:
            klass.npar_flt_nodegrp_pwse = 15. * np.ones_like(klass.npar_flt_nodegrp_plon)

        # Compute node projected-geometry
        if klass.bool_input_crs:
            try:
                gser_projected_nodes = gdf_nodes["geometry"].to_crs(crs_in)
                klass.npar_flt_nodegrp_px = gser_projected_nodes.x.to_numpy()
                klass.npar_flt_nodegrp_py = gser_projected_nodes.y.to_numpy()

            except Exception as err:
                _logger.info("Error while trying to reproject nodes in input crs")
                _logger.error(err)
                _logger.info("Reproject in default laea system")
                flt_mid_lon = 0.5 * (klass.flt_minlon + klass.flt_maxlon)
                flt_mid_lat = 0.5 * (klass.flt_minlat + klass.flt_maxlat)
                klass.npar_flt_nodegrp_px, klass.npar_flt_nodegrp_py = project(klass.npar_flt_nodegrp_plon,
                                                                               klass.npar_flt_nodegrp_plat,
                                                                               lon_0=flt_mid_lon,
                                                                               lat_0=flt_mid_lat)

        else:
            flt_mid_lon = 0.5 * (klass.flt_minlon + klass.flt_maxlon)
            flt_mid_lat = 0.5 * (klass.flt_minlat + klass.flt_maxlat)
            klass.npar_flt_nodegrp_px, klass.npar_flt_nodegrp_py = project(klass.npar_flt_nodegrp_plon,
                                                                           klass.npar_flt_nodegrp_plat,
                                                                           lon_0=flt_mid_lon,
                                                                           lat_0=flt_mid_lat)

        return klass

    @classmethod
    def from_shp(cls, reaches_shp=None, nodes_shp=None, bool_edge=True, dct_attr=None, crs_in=None):
        """Instanciate object from shapefiles

        Parameters
        ----------
        reaches_shp : str
            Full path towards reaches shapefile
        nodes_shp : str
            Full path towards nodes shapefile
        bool_edge : boolean
            Indicates if nodes segmentation includes reache edges (True) or not (False)
        dct_attr : dct
            { "reaches": { "reaches_id" : ""}, "nodes": {"reaches_id" : "", "nodes_id": "", "pwidth": "", "pwse": ""} }
            Dictionary matching input reaches/nodes attributes
        crs_in : CRS-like
            Projected CRS code in which reproject geometries to compute node xs along reaches
            If not provided, used default laea projection centered on reaches/nodes extent
        Returns
        -------
        klass : RiverGeomProduct object

        """

        # Check reaches_shp input
        if not os.path.isfile(reaches_shp):
            raise FileExistsError("Input reaches_shp file does not exist..")
        else:
            if not reaches_shp.endswith(".shp"):
                raise FileExtensionError(message="Input file is not a .shp")

        # Check reaches_shp input
        if not os.path.isfile(nodes_shp):
            raise FileExistsError("Input nodes_shp file does not exist..")
        else:
            if not nodes_shp.endswith(".shp"):
                raise FileExtensionError(message="Input file is not a .shp")

        # Load 1D geometries
        gdf_reaches = gpd.read_file(reaches_shp)
        gdf_nodes = gpd.read_file(nodes_shp)

        klass = RiverGeomProduct.from_gdf(gdf_reaches=gdf_reaches,
                                          gdf_nodes=gdf_nodes,
                                          bool_edge=bool_edge,
                                          dct_attr=dct_attr,
                                          crs_in=crs_in)

        return klass

    def __str__(self):
        """str method

        Returns
        -------
        message : str
        output for print function

        """

        message = "RiverGeom product:\n"
        message += "contains {} nodes covering {} reaches".format(self.int_node_dim,
                                                                  self.int_reach_dim)

        return message

    def get_bbox(self):
        """Return object bounding box
        """
        return self.flt_minlon, self.flt_minlat, self.flt_maxlon, self.flt_maxlat

    def draw_allreaches_centerline(self):
        """For each reach within the product, draw a centerline between the nodes
        """

        _logger.info(" ---- Draw all reaches centerlines ----")

        for int_reachid in self.npar_int_reachgrp_reachid:
            self.draw_singlereach_centerline(int_reachid)

    def sort_nodes_along_reach(self, reachid):
        """Use projected node information to sort - if necessary - nodes along reach
        as xs must be increasing for section computation

        Parameters
        ----------
        reachid : int or str
            reach_id

        Returns
        -------
        npar_int_xs_argsrt : np.array
            output from argsort over xs
        idx_first : int
            index of first node (with lowest xs)
        idx_last : int
            index of last node (with highest xs)

        """

        # Extract node projected information over current reach
        px_subset = self.npar_flt_nodegrp_px[np.where(self.npar_int_nodegrp_reachid == reachid)]
        py_subset = self.npar_flt_nodegrp_py[np.where(self.npar_int_nodegrp_reachid == reachid)]
        pxs_subset = np.array(
            [self.dct_pcenterline[reachid]["xy"].project(Point((x, y)), normalized=True) for (x, y) in
             zip(px_subset, py_subset)])

        # Sort node along reach
        npar_int_xs_argsrt = np.argsort(pxs_subset)

        # Set first (min xs) and last (max xs) node
        if npar_int_xs_argsrt[0] == 0:
            idx_first = 0
            idx_last = -1
        else:
            idx_first = -1
            idx_last = 0

        return npar_int_xs_argsrt, idx_first, idx_last

    def draw_singlereach_centerline(self, reachid):
        """Draw a centerline between the nodes along a single reach

        Parameters
        ----------
        reachid : int
            Current study reach ID

        """

        # Check input types
        if isinstance(reachid, str):
            if not isinstance(self.npar_int_nodegrp_reachid.dtype, (str, object)):
                raise TypeError(f"Input reachid of class {reachid.__class__} is different from attribute "
                                f"npar_int_nodegrp_reachid dtype {self.npar_int_nodegrp_reachid.dtype}")

        else:
            if reachid.__class__ != self.npar_int_nodegrp_reachid.dtype:
                raise TypeError(f"Input reachid of class {reachid.__class__} is different from attribute "
                                f"npar_int_nodegrp_reachid dtype {self.npar_int_nodegrp_reachid.dtype}")

        self.dct_centerline[reachid] = {}

        # Get geometries
        self.dct_centerline[reachid]["geom_xy"] = self.dct_pcenterline[reachid]["xy"]
        self.dct_centerline[reachid]["geom_lonlat"] = self.dct_pcenterline[reachid]["lonlat"]

        # Extract node projected information over current reach to sort nodes
        (npar_int_xs_argsrt, idx_first, idx_last) = self.sort_nodes_along_reach(reachid)

        # Extract other prior node information over current reach
        nodeid_subset = self.npar_int_nodegrp_nodeid[np.where(self.npar_int_nodegrp_reachid == reachid)]
        plon_subset = self.npar_flt_nodegrp_plon[np.where(self.npar_int_nodegrp_reachid == reachid)]
        plat_subset = self.npar_flt_nodegrp_plat[np.where(self.npar_int_nodegrp_reachid == reachid)]
        pwidth_subset = self.npar_flt_nodegrp_pwidth[np.where(self.npar_int_nodegrp_reachid == reachid)]
        pwse_subset = self.npar_flt_nodegrp_pwse[np.where(self.npar_int_nodegrp_reachid == reachid)]

        # Derive geographic-xs to be compliant with SW1Dto2D
        flt_centerline_length = self.dct_pcenterline[reachid]["xy"].length
        xs_subset = np.array(
            [self.dct_pcenterline[reachid]["lonlat"].project(Point((x, y)), normalized=True) for (x, y) in
             zip(plon_subset[npar_int_xs_argsrt], plat_subset[npar_int_xs_argsrt])]
        ) * flt_centerline_length

        # Prepare inputs for centerline drawing
        if self.bool_edge:
            npar_int_wrk_nodeid = nodeid_subset[npar_int_xs_argsrt]
            npar_flt_wrk_pwidth = pwidth_subset[npar_int_xs_argsrt]
            npar_flt_wrk_pwse = pwse_subset[npar_int_xs_argsrt]
            npar_flt_wrk_xs = xs_subset

        else:  # add reach edge to list of nodes
            npar_int_wrk_nodeid = np.zeros((nodeid_subset.size + 2,))
            npar_flt_wrk_pwidth = np.zeros((pwidth_subset.size + 2,))
            npar_flt_wrk_pwse = np.zeros((pwse_subset.size + 2,))
            npar_flt_wrk_xs = np.zeros((xs_subset.size + 2,))

            npar_flt_wrk_xs[1:-1] = xs_subset.copy()

            if "0" in nodeid_subset:
                _logger.warning(
                    "Warning: node_id set to '0' will match virtual edge node ids. issue will arise when dropping them.")
            npar_int_wrk_nodeid[0] = "0"
            npar_int_wrk_nodeid[1:-1] = nodeid_subset[npar_int_xs_argsrt]
            npar_int_wrk_nodeid[-1] = "0"

            npar_flt_wrk_pwidth[0] = np.mean(pwidth_subset)
            npar_flt_wrk_pwidth[1:-1] = pwidth_subset[npar_int_xs_argsrt]
            npar_flt_wrk_pwidth[-1] = np.mean(pwidth_subset)

            npar_flt_wrk_pwse[0] = np.mean(pwse_subset)
            npar_flt_wrk_pwse[1:-1] = pwse_subset[npar_int_xs_argsrt]
            npar_flt_wrk_pwse[-1] = np.mean(pwse_subset)

        # Force xs edge to match centerline edges (for use in sw1dto2d)
        npar_flt_wrk_xs[0] = 0.
        npar_flt_wrk_xs[-1] = flt_centerline_length

        # Node-level caracteristic along reach
        self.dct_centerline[reachid]["nodes_id"] = npar_int_wrk_nodeid
        self.dct_centerline[reachid]["xs"] = npar_flt_wrk_xs
        self.dct_centerline[reachid]["W"] = npar_flt_wrk_pwidth
        self.dct_centerline[reachid]["H"] = npar_flt_wrk_pwse

        # Check potential duplicate xs
        npar_idx_duplicates = np.where((npar_flt_wrk_xs[1:] - npar_flt_wrk_xs[:-1]) == 0.)[0]
        self.dct_centerline[reachid]["duplicates"] = npar_idx_duplicates

        if npar_idx_duplicates.size > 0:
            _logger.warning(f"Warning reach {reachid} : found duplicate xs")
            self.dct_centerline[reachid]["bool_duplicates"] = True

            # If last node is at edge, make sure removed duplicate is the virtual node at edge
            int_max_duplicate_idx = npar_flt_wrk_xs.size - 2
            if int_max_duplicate_idx in npar_idx_duplicates:
                _logger.warning(f"Warning reach {reachid}: last node is at edge => to debug")
                npar_idx_duplicates[-1] = npar_flt_wrk_xs.size - 1

        else:
            self.dct_centerline[reachid]["bool_duplicates"] = False

    def draw_allreaches_sections(self, type="ortho", flt_factor_width=10.):
        """ Draw all section geometries

        Parameters
        ----------
        flt_factor_width : float
            SWORD width multiplying factor to draw section
        type : str
            "ortho" : by default. sections orthogonal to centerlines
            "chck": over a given centerline: all sections are parallel and orthogonal to
            a straight line between the edges of the centerline

        Returns
        -------
        gdf_allreaches_sections : gpd.GeoDataFrame

        """

        # print(" ---- Draw all reaches sections : %s ----" % type)
        l_gdf_sections = []

        if type == "ortho":
            for reachid in self.npar_int_reachgrp_reachid:
                gdf_sections = self.draw_singlereach_sections_ortho(reachid, flt_factor_width)
                l_gdf_sections.append(gdf_sections)

        elif type == "chck":
            for reachid in self.npar_int_reachgrp_reachid:
                gdf_sections = self.draw_singlereach_sections_chck(reachid, flt_factor_width)
                l_gdf_sections.append(gdf_sections)

        else:
            raise ValueError

        gdf_allreaches_sections = pd.concat(l_gdf_sections).reset_index(drop=True)
        return gdf_allreaches_sections

    def _prepare_inputs_for_section_drawing(self, reachid, flt_factor_width=10.):
        """Format reach prior data to be compatible with class computing sections

        Parameters
        ----------
        reachid : int or str
            unique reach identifier

        Returns
        -------
        df_model1d : pd.DataFrame
            Prior reach data formatted for class SW1Dto2D instanciation
        lin_centerline_valid : LineString
            Rech line geometry in geographic coordinates

        """

        # Get attributes
        nid = self.dct_centerline[reachid]["nodes_id"]
        xs = self.dct_centerline[reachid]["xs"]
        w = self.dct_centerline[reachid]["W"] * flt_factor_width
        h = self.dct_centerline[reachid]["H"]
        df_model1d = pd.DataFrame(
            {"xs": xs,
             "id": nid,
             "H": h,
             "W": w}
        )

        # Remove potential duplicates
        if self.dct_centerline[reachid]["bool_duplicates"]:
            df_model1d.drop(labels=self.dct_centerline[reachid]["duplicates"],
                            axis=0,
                            inplace=True)
            _logger.warning(f"Warning reach {reachid} : Duplicate xs have been removed")

        # Check centerline geometry
        lin_centerline_valid = check_centerline_geometry(self.dct_pcenterline[reachid]["lonlat"])

        return df_model1d, lin_centerline_valid

    def _sort_sections(self, reachid, df_model1d, l_sections, type="ortho", npar_flt_theta=None):
        """Format section and prior reach information together while cleaning constuction node if necessary

        Parameters
        ----------
        reachid : int or str
            unique reach identifier
        df_model1d : pd.DataFrame
            Prior reach data formatted for class SW1Dto2D instanciation
        l_sections : list of LineStrings
            List of cross-section geometries sorted like df_model1d index
        type : str
            "ortho" or "chck"
        npar_flt_theta : np.array
            Array of angle in radian between ortho section and its check counterpart

        Returns
        -------
        gdf_sections : gpd.GeoDataFrame
            Clean and unique cross-section geometries with node prior information such as
                node_id : unique node identifier
                loc_xs : curvilinear abscissa or position of node along reach
        """

        dict_sections = {
            "reach_id": [reachid] * len(l_sections),
            "node_id": df_model1d["id"].to_numpy(),
            "loc_xs": df_model1d["xs"].to_numpy()
        }

        if type == "chck":
            if npar_flt_theta is None:
                _logger.error("missing theta inputs for section sorting.")
                raise ValueError
            dict_sections["theta"] = npar_flt_theta * 180. / np.pi
            dict_sections["sin_theta"] = np.sin(npar_flt_theta)

        gser_sections = gpd.GeoSeries(l_sections, crs=CRS(4326))
        df_sections = pd.DataFrame(dict_sections)
        gdf_sections = gpd.GeoDataFrame(
            df_sections,
            geometry=gser_sections,
            crs=CRS(4326)
        )

        if not self.bool_edge:
            _logger.info("Dropping virtual edge nodes")
            if isinstance(df_sections.at[0, "node_id"], str):
                idx_edge = df_sections[df_sections["node_id"] == "0"].index
            elif isinstance(df_sections.at[0, "node_id"], (int, float)):
                idx_edge = df_sections[df_sections["node_id"] == 0].index
            else:
                raise NotImplementedError("node_id must be a str or a numeric value")
            gdf_sections.drop(labels=idx_edge, axis=0, inplace=True)

        return gdf_sections

    def draw_singlereach_sections_ortho(self, reachid, flt_factor_width=10.):
        """Draw section of type "ortho" over a single reach

        Parameters
        ----------
        reachid : int or str
            unique reach identifier
        flt_factor_width : float
            SWORD width multiplying factor to draw section

        Returns
        -------
        gdf_sections : gpd.GeoDataFrame
            Clean and unique cross-section geometries ORTHOGONAL to riverline with node prior information such as
                node_id : unique node identifier
                loc_xs : curvilinear abscissa or position of node along reach
        """

        # Prepare inputs for the SW1Dto2D object
        df_model1d, lin_centerline_valid = self._prepare_inputs_for_section_drawing(reachid=reachid,
                                                                                    flt_factor_width=flt_factor_width)

        # Draw sections
        obj_sw1dto2d = SW1Dto2D(model_output_1d=df_model1d,
                                curvilinear_abscissa_key="xs",
                                heights_key="H",
                                widths_key="W",
                                centerline=lin_centerline_valid)
        obj_sw1dto2d.compute_xs_parameters(optimize_normals=False)
        l_sections = obj_sw1dto2d.compute_xs_cutlines()  # list of Linestring object

        # Format output sections
        gdf_sections = self._sort_sections(reachid=reachid,
                                           df_model1d=df_model1d,
                                           l_sections=l_sections)


        return gdf_sections

    def draw_singlereach_sections_chck(self, reachid, flt_factor_width=10.):
        """Draw sections of type "chck" over a single reach

        Parameters
        ----------
        reachid : int or str
            unique reach identifier
        flt_factor_width : float
            SWORD width multiplying factor to draw section

        Returns
        -------
        gdf_sections : gpd.GeoDataFrame
            Clean and unique cross-section geometries ALL PARALLEL with node prior information such as
                node_id : unique node identifier
                loc_xs : curvilinear abscissa or position of node along reach
        """

        # Prepare inputs for the SW1Dto2D object
        df_model1d, lin_centerline_valid = self._prepare_inputs_for_section_drawing(reachid=reachid,
                                                                                    flt_factor_width=flt_factor_width)

        # Instanciate object to derive sections
        obj_sw1dto2d = CloneSW1Dto2D(model_output_1d=df_model1d,
                                     curvilinear_abscissa_key="xs",
                                     heights_key="H",
                                     widths_key="W",
                                     centerline=lin_centerline_valid)
        obj_sw1dto2d.compute_xs_parameters(optimize_normals=False)

        # Update sections orientation
        modified_normals, angles = modified_compute_xs_parameters(obj_sw1dto2d, lin_centerline_valid)
        obj_sw1dto2d.xs_normals = modified_normals

        # Draw sections geometries
        l_sections = obj_sw1dto2d.compute_xs_cutlines()  # list of Linestring object

        # Format sections for outputs
        gdf_sections = self._sort_sections(reachid=reachid,
                                           df_model1d=df_model1d,
                                           l_sections=l_sections,
                                           type="chck",
                                           npar_flt_theta=angles)

        return gdf_sections
