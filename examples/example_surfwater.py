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

# More complete example use of BAS tools:
# - How to add a new/specific type of watermask via child classes WaterMaskCE and BASProcessorCE
# - How to compute a width product as a combination of BAS-based widths products (from scenario 0 and 11)
# - A first uncertainty model for node-scale width

import os
os.environ['USE_PYGEOS'] = '0'

from argparse import ArgumentParser
from datetime import datetime
import geopandas as gpd
import logging
import numpy as np
import pandas as pd
import rasterio as rio
import sys

from bas.rivergeomproduct import RiverGeomProduct
from bas.basprocessor import BASProcessor
from bas.watermask import WaterMask
from bas.tools import FileExtensionError
from bas.widths import compute_widths_from_single_watermask

# Config BAS
DCT_CONFIG_O = {
    "clean": {"bool_clean": True,
              "type_clean": "base",  # "base"/"waterbodies"
              "fpath_wrkdir": ".",
              "gdf_waterbodies": None
              },
    "label": {"bool_label": True,
              "type_label": "base",  # "base"
              "fpath_wrkdir": "."
              },
    "widths": {"scenario": 11
               }
}


def compute_nodescale_width(gdf_widths_ortho=None, gdf_widths_chck=None):
    """ Compute node-scale widths

    :param gdf_widths_ortho:
    :param gdf_widths_chck:
    :param method:
    :return:
    """

    LOGGER.info("Get width_1-related info : width_1 + beta")
    gdf_widths_out = gdf_widths_ortho.copy(deep=True)
    # Width from ortho sections :: width_1
    gdf_widths_out.rename(mapper={"width": "width_1"}, axis=1, inplace=True)

    LOGGER.info("Get width_2-related info : width_2 + theta")
    # Width from parallel sections : width_2
    gdf_widths_out.insert(loc=len(gdf_widths_out.columns) - 1,
                          column="width_2",
                          value=gdf_widths_chck["width"])

    # Add theta and sin_theta columns to output width dataframe
    gdf_widths_out.insert(loc=len(gdf_widths_out.columns) - 1,
                          column="theta",
                          value=gdf_widths_chck["theta"])
    gdf_widths_out.insert(loc=len(gdf_widths_out.columns) - 1,
                          column="sin_theta",
                          value=gdf_widths_chck["sin_theta"])

    gdf_widths_out.insert(loc=len(gdf_widths_out.columns) - 1,
                          column="cos_theta",
                          value=np.nan)
    gdf_widths_out["cos_theta"] = gdf_widths_chck["theta"].apply(np.cos)

    # Add final width product column
    gdf_widths_out.insert(loc=len(gdf_widths_out.columns) - 1,
                          column="width",
                          value=np.nan)

    ser_cos_theta_tmp = gdf_widths_out["theta"].apply(lambda a: np.abs(np.cos(a)))
    gdf_widths_out["width"] = gdf_widths_out["width_1"] + gdf_widths_out["width_2"].mul(ser_cos_theta_tmp)
    gdf_widths_out["width"] = gdf_widths_out["width"].div(ser_cos_theta_tmp.apply(lambda ct: ct + 1.))

    del ser_cos_theta_tmp

    return gdf_widths_out


def compute_nodescale_widtherror(gdf_widths, flt_watermask_resol):
    """Compute node-scale uncertainty

    :param gdf_widths:
    :param flt_watermask_resol:
    :param method:
    :return:
    """

    # Observation/measurement error
    # related to the edges of the watermask
    ser_sigo = gdf_widths["nb_banks"].mul(flt_watermask_resol / np.sqrt(12.))

    # Representativeness error
    ser_sigr = pd.Series(flt_watermask_resol / np.sqrt(12.) * np.ones((len(gdf_widths),)), index=gdf_widths.index)

    # Structure error : simplified version
    ser_sigs = pd.Series(flt_watermask_resol / np.sqrt(12.) * np.ones((len(gdf_widths),)), index=gdf_widths.index)

    ser_errtot = np.sqrt(ser_sigo ** 2.0 + ser_sigr ** 2.0 + ser_sigs ** 2.0)

    return ser_errtot, ser_sigo, ser_sigr, ser_sigs


class WaterMaskCE(WaterMask):

    def __init__(self):

        # Init parent class
        super().__init__()

    @classmethod
    def from_surfwater(cls, surfwater_tif=None):
        """ Instanciate from a Surfwater watermask

        Parameters
        ----------
        surfwater_tif : str
            Full path to surfwater mask stored in a GeoTiff file
        """

        if not os.path.isfile(surfwater_tif):
            LOGGER.error("Watermask.from_surfwater: Input tif file does not exist..")
            raise FileExistsError("Input tif file does not exist..")

        else:
            if not surfwater_tif.endswith(".tif"):
                raise FileExtensionError

        klass = WaterMaskCE()
        klass.origin = "SurfWater"
        klass.rasterfile = surfwater_tif
        klass.coordsyst = "proj"

        with rio.open(surfwater_tif, 'r') as src:
            klass.crs = src.crs
            klass.crs_epsg = src.crs.to_epsg()
            klass.res = src.transform.a
            klass.bbox = (src.bounds.left,
                          src.bounds.bottom,
                          src.bounds.right,
                          src.bounds.top)

        return klass


class BASProcessorCE(BASProcessor):

    def __init__(self,
                 str_watermask_tif=None,
                 gdf_sections=None,
                 gdf_reaches=None,
                 attr_reachid=None,
                 str_proj="proj",
                 str_provider=None,
                 str_datetime=None):

        # Init parent class
        super().__init__(str_watermask_tif,
                         gdf_sections,
                         gdf_reaches,
                         attr_reachid,
                         str_proj,
                         str_provider,
                         str_datetime)

    def preprocessing(self):
        """Preprocessing: load watermask, reproject sections et check bounding boxes intersections
            """

        print("----- WidthProcessing = Preprocessing -----")
        print("")

        # Load WaterMask object
        if self.provider == "SW":  # Watermask is provided by Surfwater
            LOGGER.info("Load watermask from Surfwater")
            self.watermask = WaterMaskCE.from_surfwater(self.f_watermask_in)
        else:
            LOGGER.error(f"Provider {self.provider} not recognized")
            raise NotImplementedError
        LOGGER.info("Watermask loaded..")

        # Reproject sections to watermask coordinate system
        LOGGER.info("Reproject sections to watermask coordinate system")
        self.gdf_reaches = self.gdf_reaches.to_crs(self.watermask.crs_epsg)
        self.gdf_sections = self.gdf_sections.to_crs(self.watermask.crs_epsg)
        LOGGER.info("Reproject sections to watermask coordinate system done ..")

        # Check boundingbox compatibility
        self.check_bbox_compatibility()

        print("")
        print("----- Preprocessing : Done -----")


class WidthProcessor:

    def __init__(self, str_watermask_tif=None,
                 str_datetime=None,
                 str_reaches_shp=None,
                 str_nodes_shp=None):
        """Class constructor
        """

        LOGGER.info("Instanciate WidthProcessor")

        if str_watermask_tif is not None:
            self.f_watermask_in = str_watermask_tif
            if not os.path.isfile(str_watermask_tif):
                LOGGER.error("Input watermask GeoTiff does not exist")
                raise FileExistsError("Input watermask GeoTiff does not exist")
        else:
            LOGGER.error("Missing watermask GeoTiff input file")
            raise ValueError("Missing watermask GeoTiff input file")
        self.scene_name = os.path.basename(self.f_watermask_in).split(".")[0]

        if str_datetime is not None:
            self.scene_datetime = str_datetime
            try:
                dt_scene = datetime.strptime(self.scene_datetime, "%Y%m%dT%H%M%S")
            except ValueError:
                LOGGER.error("input datetime {} does not match format '%Y%m%dT%H%M%S'.".format(self.scene_datetime))
                raise ValueError("input datetime {} does not match format '%Y%m%dT%H%M%S'.".format(self.scene_datetime))
        else:
            LOGGER.error("Missing scene datetime information input")
            raise ValueError("Missing scene datetime information input")

        if str_reaches_shp is not None:
            if not os.path.isfile(str_reaches_shp):
                LOGGER.error("Input reaches shapefile does not exist")
                raise FileExistsError("Input reaches shapefile does not exist")
            self.reaches_shp = str_reaches_shp
            self.gdf_reaches = gpd.read_file(self.reaches_shp)
        else:
            LOGGER.error("Input reaches shapefile does not exist")
            raise FileExistsError("Input reaches shapefile does not exist")

        if str_nodes_shp is not None:
            if not os.path.isfile(str_nodes_shp):
                LOGGER.error("Input nodes shapefile does not exist")
                raise FileExistsError("Input nodes shapefile does not exist")
            self.nodes_shp = str_nodes_shp
        else:
            LOGGER.error("Input nodes shapefile does not exist")
            raise FileExistsError("Input nodes shapefile does not exist")

        self.gdf_sections_ortho = None
        self.gdf_sections_chck = None

        self.bas_processor_o = None
        self.bas_processor_c = None

        self.gdf_nodescale_widths = None

    def preprocessing(self, flt_factor_width=10.):
        """ Prepare SWOT-like product watermask to process and associated SWORD nodes and reaches
        """

        LOGGER = logging.getLogger('WidthProcessing.preprocessing')

        # Instanciate RiverGeom object
        LOGGER.info("Instanciate RiverGeomProduct object ..")
        try:
            dct_geom_attr = {"reaches": {"reaches_id": "reach_id"},
                             "nodes": {"reaches_id": "reach_id",
                                       "nodes_id": "node_id",
                                       "pwidth": "p_width",
                                       "pwse": "p_wse"}}
            obj_rivergeom = RiverGeomProduct.from_shp(reaches_shp=self.reaches_shp,
                                                      nodes_shp=self.nodes_shp,
                                                      bool_edge=False,
                                                      dct_attr=dct_geom_attr)
            LOGGER.info("Instanciation done..")
        except Exception as err:
            LOGGER.error(err)
            LOGGER.error("Instanciate RiverGeomProduct object KO ..")
            raise Exception

        # Set centerlines for section definition
        LOGGER.info("Set centerlines for section definition")
        try:
            obj_rivergeom.draw_allreaches_centerline()
        except Exception as err:
            LOGGER.error(err)
            LOGGER.error("Set centerlines for section definition KO ..")
            raise Exception


        # Traditionnal orthogonal sections
        LOGGER.info("Traditionnal orthogonal sections ..")
        try:
            self.gdf_sections_ortho = obj_rivergeom.draw_allreaches_sections(type="ortho",
                                                                             flt_factor_width=flt_factor_width)
        except Exception as err:
            LOGGER.error(err)
            LOGGER.error("Draw orthogonal sections KO ..")
            raise Exception

        # Checking parallel sections
        LOGGER.info("Parallel sections given the main direction of the reach")
        try:
            self.gdf_sections_chck = obj_rivergeom.draw_allreaches_sections(type="chck",
                                                                            flt_factor_width=flt_factor_width)
        except Exception as err:
            LOGGER.error(err)
            LOGGER.error("Draw parallel sections KO ..")
            raise Exception

        # Instanciate BASProcessorCalVal objects
        try:
            LOGGER.info("Instanciate BASProcessor object for sections_ortho")
            self.bas_processor_o = BASProcessorCE(
                str_watermask_tif=self.f_watermask_in,
                gdf_sections=self.gdf_sections_ortho,
                gdf_reaches=self.gdf_reaches,
                attr_reachid="reach_id",
                str_proj="proj",
                str_provider="SW"
            )
            self.bas_processor_o.preprocessing()
        except Exception as err:
            LOGGER.error(err)
            LOGGER.error("Instanciate width processor ortho KO ..")
            raise Exception

        try:
            LOGGER.info("Instanciate BASProcessor object for sections_chck")
            self.bas_processor_c = BASProcessorCE(
                str_watermask_tif=self.f_watermask_in,
                gdf_sections=self.gdf_sections_chck,
                gdf_reaches=self.gdf_reaches,
                attr_reachid="reach_id",
                str_proj="proj",
                str_provider="SW"
            )
            self.bas_processor_c.preprocessing()
        except Exception as err:
            LOGGER.error(err)
            LOGGER.error("Instanciate width processor check KO ..")
            raise Exception

    def processing(self,
                   out_dir=".",
                   str_pekel_shp=None,
                   str_type_clean=None,
                   str_type_label=None):
        """ Produce riverwidth derived from watermask
        """

        LOGGER.info("Processing based on orthogonal sections")
        try:
            dct_cfg_o = DCT_CONFIG_O

            # Set cleaning and labelling method
            LOGGER.info(f"Set cleaning method to : {str_type_clean}")
            if str_type_clean is not None:
                if str_type_clean not in ["base", "waterbodies"]:
                    LOGGER.info("Undefined cleaning type .. Ignore .. use default value")
                else:
                    dct_cfg_o["clean"]["type_clean"] = str_type_clean

            LOGGER.info(f"Set labelling method to : {str_type_label}")
            if str_type_label is not None:
                if str_type_label not in ["base"]:
                    LOGGER.info("Undefined labelling type .. Ignore .. use default value")
                else:
                    dct_cfg_o["label"]["type_label"] = str_type_label

            # Set waterbody mask and output directory
            if str_type_clean == "waterbodies":
                gdf_waterbodies = gpd.read_file(str_pekel_shp)
                dct_cfg_o["clean"]["gdf_waterbodies"] = gdf_waterbodies
            dct_cfg_o["clean"]["fpath_wrkdir"] = out_dir
            dct_cfg_o["label"]["fpath_wrkdir"] = out_dir

            gdf_widths_ortho = self.bas_processor_o.processing(dct_cfg_o)
        except Exception as err:
            LOGGER.error(err)
            LOGGER.error("Processing based on orthogonal section KO ..")
            raise Exception

        LOGGER.info("Processing based on paralell sections")
        try:
            dct_cfg_c = {"clean": {"bool_clean": False},
                         "label": {"bool_label": False},
                         "widths": {"scenario": 0
                                    }
                         }
            self.bas_processor_c.f_watermask_in = self.bas_processor_o.watermask.rasterfile
            self.bas_processor_c.watermask.rasterfile = self.bas_processor_o.watermask.rasterfile
            self.bas_processor_c.watermask.bool_labelled = True

            gdf_widths_chck = self.bas_processor_c.processing(dct_cfg_c)
        except Exception as err:
            LOGGER.error(err)
            LOGGER.error("Processing based on parallel section KO ..")
            raise Exception

        # Compute node-scale widths
        LOGGER.info("Compute widths at node scale")
        try:
            self.gdf_nodescale_widths = compute_nodescale_width(gdf_widths_ortho,
                                                                gdf_widths_chck)

            self.gdf_nodescale_widths.insert(loc=2, column="provider", value="SW")
            self.gdf_nodescale_widths.insert(loc=3, column="bool_ko", value=0)
            self.gdf_nodescale_widths["bool_ko"] = self.gdf_nodescale_widths["width"].apply(
                lambda w: np.logical_or(np.isnan(w), w == 0))

            LOGGER.info("Node-scale widths computed..")
        except Exception as err:
            LOGGER.error(err)
            LOGGER.info("Node-scale widths computed ko..")
            raise Exception

        # Compute node-scale width errors
        LOGGER.info("Compute width error at node scale")
        try:
            ser_errtot, ser_sigo, ser_sigr, ser_sigs = compute_nodescale_widtherror(self.gdf_nodescale_widths,
                                                                                    self.bas_processor_o.watermask.res)
            self.gdf_nodescale_widths.insert(loc=3, column="width_u", value=ser_errtot)
            self.gdf_nodescale_widths.insert(loc=4, column="sigo", value=ser_sigo)
            self.gdf_nodescale_widths.insert(loc=5, column="sigr", value=ser_sigr)
            self.gdf_nodescale_widths.insert(loc=6, column="sigs", value=ser_sigs)

            LOGGER.info("Node-scale width errors computed..")
        except Exception as err:
            LOGGER.error(err)
            LOGGER.info("Node-scale width errors computed ko..")
            raise Exception

    def postprocessing(self, output_dir="."):
        """ Save processed riverwidths into files

        """

        LOGGER = logging.getLogger('WidthProcessing.postprocessing')

        LOGGER.info("Save node-scale width as csv ")
        df_nodes_width = self.gdf_nodescale_widths.loc[:,
                         ["reach_id",
                          "node_id",
                          "provider",
                          "width",
                          "width_u"]].copy(deep=True)
        df_nodes_width = df_nodes_width.dropna()
        df_nodes_width["datetime"] = self.scene_datetime
        width_nodes_csv = self.scene_name + "_nodescale_BAS_widths.csv"
        df_nodes_width.to_csv(os.path.join(output_dir, width_nodes_csv))
        LOGGER.info("Node-scale width saved to csv..")

        # Save cross-sections with associated width as shp in epsg:4326
        LOGGER.info("Save cross-sections (node-scale) with associated width as shp in epsg:4326")
        width_nodes_shp = self.scene_name + "_nodescale_BAS_widths.shp"
        self.gdf_nodescale_widths.insert(loc=len(self.gdf_nodescale_widths.columns) - 1,
                                         column="datetime",
                                         value=self.scene_datetime)
        self.gdf_nodescale_widths["reach_id"] = self.gdf_nodescale_widths["reach_id"].astype(str)
        self.gdf_nodescale_widths["node_id"] = self.gdf_nodescale_widths["node_id"].astype(int).astype(str)
        self.gdf_nodescale_widths.to_file(os.path.join(output_dir, width_nodes_shp))
        LOGGER.info("Node-scale width saved to shp..")


def process_single_scene(str_watermask_tif=None,
                         str_scene_datetime=None,
                         str_reaches_shp=None,
                         str_nodes_shp=None,
                         flt_factor_width=10.,
                         str_outputdir=None,
                         str_type_clean=None,
                         str_type_label=None):
    LOGGER.info("=== Processing watermask: " + str_watermask_tif + " === : start\n")

    # watermask filename to process
    if not os.path.isfile(str_watermask_tif):
        LOGGER.error("Watermask file '{}' seems not to exist..".format(str_watermask_tif))
    str_scn_name = os.path.basename(str_watermask_tif).split(".")[0]

    # Width processing
    try:

        # Instanciate width processing class
        obj_widthprocessor = WidthProcessor(
            str_watermask_tif=str_watermask_tif,
            str_datetime=str_scene_datetime,
            str_reaches_shp=str_reaches_shp,
            str_nodes_shp=str_nodes_shp
        )

        obj_widthprocessor.preprocessing(flt_factor_width)

        # Processing
        obj_widthprocessor.processing(out_dir=str_outputdir,
                                      str_type_clean=str_type_clean,
                                      str_type_label=str_type_label)

        # PostProcessing
        obj_widthprocessor.postprocessing(str_outputdir)





    except:
        LOGGER.info("===> Fail during working with WidthProcessor object\n")
        pass


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Extract river width at SWORD nodes and reaches from a watermask")

    parser.add_argument(
        "-w", "--watermask_tif", type=str, default=None, help="Full path of the input watermask as GeoTiff to process")
    parser.add_argument(
        "-dt", "--scn_datetime", type=str, default="20200111T105853", help="Scene datetime")
    parser.add_argument(
        "-r", "--reaches_shp", type=str, default=None,
        help="Full path to shp with reach-line geometries")
    parser.add_argument(
        "-n", "--nodes_shp", type=str, default=None,
        help="Full path to shp with node-point geometries")
    parser.add_argument(
        "-tc", "--type_clean", type=str, default="base",
        help="Watermask cleaning procedure")
    parser.add_argument(
        "-tl", "--type_label", type=str, default="base",
        help="Watermask labelling procedure")
    parser.add_argument(
        "-fw", "--factor_width", type=float, default=10.,
        help="Multiplying factor applied to PRD width to draw section")
    parser.add_argument(
        "-o", "--outputdir", type=str, default=".", help="Full path to directory where to store logfiles and outputs.")
    parser.add_argument(
        "-l", "--loglevel", type=str.lower, default="info", help="Logging level: debug, info, warning, error")

    # Set input arguments
    args = parser.parse_args()

    # Set logs
    logformat = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    str_now = datetime.now().strftime("%Y%m%dT%H%M%S")
    logfile = os.path.join(args.outputdir, 'logfile_' + str_now + '.log')

    LOGGER = logging.getLogger('BAS PROCESSING')
    loglevel = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR}[args.loglevel]
    logging.basicConfig(
        level=loglevel,
        format=logformat,
        handlers=[
            logging.FileHandler(logfile),
            logging.StreamHandler(sys.stdout)
        ]
    )

    # Set output dir
    if args.outputdir == ".":
        LOGGER.warning("Output directory pass to argument does not exist, write log in current directory")

    if args.watermask_tif is None or args.reaches_shp is None or args.nodes_shp is None:
        str_err = "Missing one or more input arguments for single scene processing. watermask:{}, reaches_shp:{}, nodes_shp:{}".format(
            args.watermask_tif,
            args.reaches_shp,
            args.nodes_shp)
        LOGGER.error(str_err)
        raise ValueError(str_err)

    process_single_scene(str_watermask_tif=args.watermask_tif,
                         str_scene_datetime=args.scn_datetime,
                         str_reaches_shp=args.reaches_shp,
                         str_nodes_shp=args.nodes_shp,
                         str_outputdir=args.outputdir,
                         str_type_clean=args.type_clean,
                         str_type_label=args.type_label,
                         flt_factor_width=args.factor_width
                         )
