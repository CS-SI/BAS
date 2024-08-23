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
module rivergeomproduct.py
: Contains classes to manipulate river-related 1D-geometries
"""

import geopandas as gpd
import numpy as np
import os
import pandas as pd
from pyproj import CRS
from shapely.geometry import Point, LineString
from sw1dto2d.sw1dto2d import SW1Dto2D
from tools import FileExtensionError
from tools import project

# os.environ['USE_PYGEOS'] = '0'


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

    def __init__(self):
        """Class constructor
        """

        self.bool_isempty = True
        self.bool_edge = True

        # 1D geometries bounding box
        self.flt_minlon = None
        self.flt_minlat = None
        self.flt_maxlon = None
        self.flt_maxlat = None

        # 1D geometries counts
        self.int_node_dim = None  # Number of nodes within the product
        self.int_reach_dim = None  # Number of reaches within the product

        # Reach geometries information
        self.npar_int_reachgrp_reachid = None  # Reach ID
        self.npar_flt_reachgrp_plon = None
        self.npar_flt_reachgrp_plat = None
        self.dct_pcenterline = {}  # Reach LineString prior geometry defining

        # Node geometries information
        self.npar_int_nodegrp_nodeid = None  # Node id sequence
        self.npar_int_nodegrp_reachid = None  # Reach ID sequence for mapping
        self.npar_flt_nodegrp_plon = None  # Node prior lon
        self.npar_flt_nodegrp_plat = None  # Node prior lat
        self.npar_flt_nodegrp_pwidth = None  # Prior node width
        self.npar_flt_nodegrp_pwse = None  # Prior node wse

        # To produce
        self.npar_flt_nodegrp_px = None  # Node prior longitude projected to "laea"
        self.npar_flt_nodegrp_py = None  # Node prior latitude projected to "laea"
        self.dct_centerline = {}  # Reach-scale geometry/xs/p_wse/p_width -- use to derive sections from sw1dto2d

    @classmethod
    def from_shp(cls, reaches_shp=None, nodes_shp=None, bool_edge=True, dct_attr=None):
        """Instanciate object from shapefiles

        Parameters
        ----------
        reaches_shp : str
        nodes_shp : str
        dct_attr : dct
        { "reaches": { "reaches_id" : ""}, "nodes": {"reaches_id" : "", "nodes_id": "", "pwidth": "", "pwse": ""} }
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

        klass = RiverGeomProduct()

        # Load 1D geometries
        gdf_reaches = gpd.read_file(reaches_shp)
        gdf_nodes = gpd.read_file(nodes_shp)

        # Count available geometries
        klass.int_reach_dim = len(gdf_reaches)
        klass.int_node_dim = len(gdf_nodes)
        if klass.int_reach_dim > 0 and klass.int_node_dim > 0:
            klass.isempty = False
        klass.bool_edge = bool_edge

        # Get reaches information
        klass.npar_int_reachgrp_reachid = gdf_reaches[dct_attr["reaches"]["reaches_id"]].to_numpy()
        klass.minlon, klass.minlat, klass.maxlon, klass.maxlat = gdf_reaches.total_bounds

        gser_reach_midpoint = gdf_reaches["geometry"].interpolate(0.5, normalized=True)
        klass.npar_flt_reachgrp_plon = gser_reach_midpoint.x.to_numpy()
        klass.npar_flt_reachgrp_plat = gser_reach_midpoint.y.to_numpy()
        for index, row in gdf_reaches.iterrows():
            reach_id = row[dct_attr["reaches"]["reaches_id"]]

            klass.dct_pcenterline[reach_id] = {}
            klass.dct_pcenterline[reach_id]["lonlat"] = row["geometry"]
            klass.dct_pcenterline[reach_id]["plon"] = klass.npar_flt_reachgrp_plon[index]
            klass.dct_pcenterline[reach_id]["plat"] = klass.npar_flt_reachgrp_plat[index]

            arr_centerline_lon = [t[0] for t in row["geometry"].coords]
            arr_centerline_lat = [t[1] for t in row["geometry"].coords]
            arr_centerline_x, arr_centerline_y = project(arr_centerline_lon,
                                                         arr_centerline_lat,
                                                         lon_0=0.5 * (klass.minlon + klass.maxlon),
                                                         lat_0=0.5 * (klass.minlat + klass.maxlat))
            klass.dct_pcenterline[reach_id]["xy"] = LineString(
                [(x, y) for (x, y) in
                 zip(arr_centerline_x, arr_centerline_y)])

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

        klass.npar_flt_nodegrp_px, klass.npar_flt_nodegrp_py = project(klass.npar_flt_nodegrp_plon,
                                                                       klass.npar_flt_nodegrp_plat,
                                                                       lon_0=0.5 * (klass.minlon + klass.maxlon),
                                                                       lat_0=0.5 * (klass.minlat + klass.maxlat))

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

        print(" ---- Draw all reaches centerlines ----")

        for int_reachid in self.npar_int_reachgrp_reachid:
            self.draw_singlereach_centerline(int_reachid)

    def draw_singlereach_centerline(self, reachid):
        """Draw a centerline between the nodes along a single reach

        Parameters
        ----------
        reachid : int
            Current study reach ID

        """

        # Check input types
        if isinstance(reachid, str) :
            if not isinstance(self.npar_int_nodegrp_reachid.dtype, (str,object)):
                raise TypeError(f"Input reachid of class {reachid.__class__} is different from attribute "
                                f"npar_int_nodegrp_reachid dtype {self.npar_int_nodegrp_reachid.dtype}")

        else:
            if reachid.__class__ != self.npar_int_nodegrp_reachid.dtype:
                raise TypeError(f"Input reachid of class {reachid.__class__} is different from attribute "
                            f"npar_int_nodegrp_reachid dtype {self.npar_int_nodegrp_reachid.dtype}")

        self.dct_centerline[reachid] = {}

        # Sort node along reach
        px_subset = self.npar_flt_nodegrp_px[np.where(self.npar_int_nodegrp_reachid == reachid)]
        py_subset = self.npar_flt_nodegrp_py[np.where(self.npar_int_nodegrp_reachid == reachid)]
        pxs_subset = np.array(
            [self.dct_pcenterline[reachid]["xy"].project(Point((x, y)), normalized=True) for (x, y) in
             zip(px_subset, py_subset)])
        npar_int_xs_argsrt = np.argsort(pxs_subset)
        if npar_int_xs_argsrt[0] == 0:
            idx_first = 0
            idx_last = -1
        else:
            idx_first = -1
            idx_last = 0

        # Extract basic node information
        nodeid_subset = self.npar_int_nodegrp_nodeid[np.where(self.npar_int_nodegrp_reachid == reachid)]
        plon_subset = self.npar_flt_nodegrp_plon[np.where(self.npar_int_nodegrp_reachid == reachid)]
        plat_subset = self.npar_flt_nodegrp_plat[np.where(self.npar_int_nodegrp_reachid == reachid)]
        pwidth_subset = self.npar_flt_nodegrp_pwidth[np.where(self.npar_int_nodegrp_reachid == reachid)]
        pwse_subset = self.npar_flt_nodegrp_pwse[np.where(self.npar_int_nodegrp_reachid == reachid)]

        # Prepare inputs for centerline drawing
        if self.bool_edge:

            npar_int_wrk_nodeid = nodeid_subset[npar_int_xs_argsrt]
            npar_flt_wrk_pwidth = pwidth_subset[npar_int_xs_argsrt]
            npar_flt_wrk_pwse = pwse_subset[npar_int_xs_argsrt]

            # Get prior node coords
            npar_flt_wrk_plon = plon_subset[npar_int_xs_argsrt]
            npar_flt_wrk_plat = plat_subset[npar_int_xs_argsrt]

        else:  # add reach edge to list of nodes
            npar_int_wrk_nodeid = np.zeros((nodeid_subset.size + 2,))
            npar_int_wrk_nodeid[0] = "0"
            npar_int_wrk_nodeid[1:-1] = nodeid_subset[npar_int_xs_argsrt]
            npar_int_wrk_nodeid[-1] = "0"

            pwidth_subset = self.npar_flt_nodegrp_pwidth[np.where(self.npar_int_nodegrp_reachid == reachid)]
            npar_flt_wrk_pwidth = np.zeros((pwidth_subset.size + 2,))
            npar_flt_wrk_pwidth[0] = np.mean(pwidth_subset)
            npar_flt_wrk_pwidth[1:-1] = pwidth_subset[npar_int_xs_argsrt]
            npar_flt_wrk_pwidth[-1] = np.mean(pwidth_subset)

            pwse_subset = self.npar_flt_nodegrp_pwse[np.where(self.npar_int_nodegrp_reachid == reachid)]
            npar_flt_wrk_pwse = np.zeros((pwse_subset.size + 2,))
            npar_flt_wrk_pwse[0] = np.mean(pwse_subset)
            npar_flt_wrk_pwse[1:-1] = pwse_subset[npar_int_xs_argsrt]
            npar_flt_wrk_pwse[-1] = np.mean(pwse_subset)

            # Get prior node coords
            npar_flt_wrk_plon = np.zeros((plon_subset.size + 2,))
            npar_flt_wrk_plon[0] = list(self.dct_pcenterline[reachid]["lonlat"].coords)[idx_first][0]
            npar_flt_wrk_plon[1:-1] = plon_subset[npar_int_xs_argsrt]
            npar_flt_wrk_plon[-1] = list(self.dct_pcenterline[reachid]["lonlat"].coords)[idx_last][0]

            npar_flt_wrk_plat = np.zeros((plat_subset.size + 2,))
            npar_flt_wrk_plat[0] = list(self.dct_pcenterline[reachid]["lonlat"].coords)[idx_first][1]
            npar_flt_wrk_plat[1:-1] = plat_subset[npar_int_xs_argsrt]
            npar_flt_wrk_plat[-1] = list(self.dct_pcenterline[reachid]["lonlat"].coords)[idx_last][1]

        # Get geometries
        self.dct_centerline[reachid]["geom_xy"] = self.dct_pcenterline[reachid]["xy"]
        self.dct_centerline[reachid]["geom_lonlat"] = self.dct_pcenterline[reachid]["lonlat"]

        # Derive xs
        flt_centerline_length = self.dct_pcenterline[reachid]["xy"].length
        subset_xs = np.array(
            [self.dct_pcenterline[reachid]["lonlat"].project(Point((x, y)), normalized=True) for (x, y) in
             zip(npar_flt_wrk_plon, npar_flt_wrk_plat)]
        ) * flt_centerline_length
        npar_flt_wrk_xs = np.copy(subset_xs)
        npar_flt_wrk_xs[0] = 0.  # Force first node xs at 0 for use of sw1dto2d
        npar_flt_wrk_xs[-1] = flt_centerline_length

        # Node-level caracteristic along reach
        self.dct_centerline[reachid]["nodes_id"] = npar_int_wrk_nodeid
        self.dct_centerline[reachid]["xs"] = npar_flt_wrk_xs
        self.dct_centerline[reachid]["W"] = npar_flt_wrk_pwidth
        self.dct_centerline[reachid]["H"] = npar_flt_wrk_pwse

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

    def draw_singlereach_sections_ortho(self, reachid, flt_factor_width=10.):
        """Draw section of type "ortho" over a single reach

        Parameters
        ----------
        flt_factor_width : float
            SWORD width multiplying factor to draw section
        reachid :

        Returns
        -------
        flt_factor_width : float
            SWORD width multiplying factor to draw section
        gdf_sections : gpd.GeoDataFrame

        """

        # Prepare inputs for the SW1Dto2D object
        xs = self.dct_centerline[reachid]["xs"]
        W = self.dct_centerline[reachid]["W"] * flt_factor_width
        H = self.dct_centerline[reachid]["H"]
        df_model1d = pd.DataFrame(
            {"xs": xs,
             "H": H,
             "W": W}
        )

        # Check centerline geometry
        lin_centerline_valid = check_centerline_geometry(self.dct_pcenterline[reachid]["lonlat"])

        # Draw sections
        obj_sw1dto2d = SW1Dto2D(model_output_1d=df_model1d,
                                curvilinear_abscissa_key="xs",
                                heights_key="H",
                                widths_key="W",
                                centerline=lin_centerline_valid)
        obj_sw1dto2d.compute_xs_parameters(optimize_normals=False)
        l_sections = obj_sw1dto2d.compute_xs_cutlines()  # list of Linestring object

        # Format output sections
        dict_sections = {
            "reach_id": [reachid] * (len(l_sections) - 2),
            "node_id": self.dct_centerline[reachid]["nodes_id"][1:-1],
            "loc_xs": xs[1:-1]
        }
        df_sections = pd.DataFrame(dict_sections)
        gdf_sections = gpd.GeoDataFrame(
            df_sections,
            geometry=gpd.GeoSeries(l_sections[1:-1], crs=CRS(4326)),
            crs=CRS(4326)
        )

        return gdf_sections

    def draw_singlereach_sections_chck(self, reachid, flt_factor_width=10.):
        """Draw sections of type "chck" over a single reach

        Parameters
        ----------
        flt_factor_width : float
        reachid :

        Returns
        -------
        gdf_sections : gpd.GeoDataFrame

        """

        # Prepare inputs for the SW1Dto2D object
        xs = self.dct_centerline[reachid]["xs"]
        W = self.dct_centerline[reachid]["W"] * flt_factor_width
        H = self.dct_centerline[reachid]["H"]
        df_model1d = pd.DataFrame(
            {"xs": xs,
             "H": H,
             "W": W}
        )

        # Check centerline geometry
        lin_centerline_valid = check_centerline_geometry(self.dct_pcenterline[reachid]["lonlat"])

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
        dict_sections = {
            "reach_id": [reachid] * (len(l_sections) - 2),
            "node_id": self.dct_centerline[reachid]["nodes_id"][1:-1],
            "loc_xs": xs[1:-1],
            "theta": angles[1:-1] * 180. / np.pi,
            "sin_theta": np.sin(angles)[1:-1]
        }
        df_sections = pd.DataFrame(dict_sections)
        gdf_sections = gpd.GeoDataFrame(
            df_sections,
            geometry=gpd.GeoSeries(l_sections[1:-1], crs=CRS(4326)),
            crs=CRS(4326)
        )

        return gdf_sections
