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

import sys
sys.path.append("/home/cemery/Work/git/BAS/bas")

import geopandas as gpd

from basprocessor import BASProcessor
from rivergeomproduct import RiverGeomProduct

# Input file
watermask_tif = "example_watermask.tif"
ref_watermask_tif = "example_ref_waterbodies.shp"

# Simple example : sections are ready to use
shp_reaches_smpl = "example_reaches_simple.shp"
shp_sections_smpl = "example_sections_simple.shp"

# Complex example : sections have to be derived
shp_reaches_cplx = "example_reaches_cplx.shp"
shp_nodes_cplx = "example_nodes_cplx.shp"


def example_1():
    """Example_1 :
        No watermask cleaning and no watermask labelling
        Sections available
        Width only
    """

    print("===== BASProcessing Example #1 = BEGIN =====")
    print("")

    # Set config #1
    dct_cfg_V1 = {"clean": {"bool_clean": False,
                            "type_clean": "base",
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

    # Instanciate basprocessor
    processor = BASProcessor(
        str_watermask_tif=watermask_tif,
        gdf_sections=gdf_sections,
        gdf_reaches=gdf_reaches,
        attr_reachid="id",
        str_proj="proj",
        str_provider="EO"
    )
    processor.preprocessing()
    gdf_widths, _ = processor.processing(dct_cfg_V1)
    print(gdf_widths)
    # gdf_widths.to_file("widths_example1.shp")

    print("")
    print("===== BASProcessing Example #1 = END =====")


def example_2():
    """Example_2 :
        Basic watermask cleaning without reference waterbodies and no watermask labelling
        Sections available
        Width + estimation of intersection width other sections
    """

    print("===== BASProcessing Example #2 = BEGIN =====")
    print("")

    # Set config #2
    dct_cfg_V2 = {"clean": {"bool_clean": True,
                            "type_clean": "base",
                            "fpath_wrkdir": ".",
                            "gdf_waterbodies": None
                            },
                  "label": {"bool_label": False,
                            "type_label": "base",
                            "fpath_wrkdir": "."
                            },
                  "widths": {"scenario": 11
                             }
                  }

    # Instanciate basprocessor
    processor = BASProcessor(
        str_watermask_tif=watermask_tif,
        gdf_sections=gdf_sections,
        gdf_reaches=gdf_reaches,
        attr_reachid="id",
        str_proj="proj",
        str_provider="EO"
    )
    processor.preprocessing()
    gdf_widths, _ = processor.processing(dct_cfg_V2)
    print(gdf_widths)
    gdf_widths.to_file("widths_example2.shp")

    print("")
    print("===== BASProcessing Example #2 = END =====")


def example_3():
    """Example_3 :
        Watermask cleaning with reference waterbodies and no watermask labelling
        Sections available
        Width only
    """

    print("===== BASProcessing Example #3 = BEGIN =====")
    print("")

    # Set config #3
    dct_cfg_V3 = {"clean": {"bool_clean": True,
                            "type_clean": "waterbodies",
                            "fpath_wrkdir": ".",
                            "gdf_waterbodies": gdf_waterbodies
                            },
                  "label": {"bool_label": False,
                            "type_label": "base",
                            "fpath_wrkdir": "."
                            },
                  "widths": {"scenario": 0
                             }
                  }

    # Instanciate basprocessor
    processor = BASProcessor(
        str_watermask_tif=watermask_tif,
        gdf_sections=gdf_sections,
        gdf_reaches=gdf_reaches,
        attr_reachid="id",
        str_proj="proj",
        str_provider="EO"
    )
    processor.preprocessing()
    gdf_widths, _ = processor.processing(dct_cfg_V3)
    print(gdf_widths)
    gdf_widths.to_file("widths_example3.shp")

    print("")
    print("===== BASProcessing Example #3 = END =====")


def example_4():
    """Example_4 :
        Watermask cleaning with reference waterbodies + watermask labelling
        Sections available
        Width only
    """

    print("===== BASProcessing Example #4 = BEGIN =====")
    print("")

    # Set config #4
    dct_cfg_V4 = {"clean": {"bool_clean": True,
                            "type_clean": "waterbodies",
                            "fpath_wrkdir": ".",
                            "gdf_waterbodies": gdf_waterbodies
                            },
                  "label": {"bool_label": True,
                            "type_label": "base",
                            "fpath_wrkdir": "."
                            },
                  "widths": {"scenario": 0
                             }
                  }

    # Instanciate basprocessor
    processor = BASProcessor(
        str_watermask_tif=watermask_tif,
        gdf_sections=gdf_sections,
        gdf_reaches=gdf_reaches,
        attr_reachid="id",
        str_proj="proj",
        str_provider="EO"
    )
    processor.preprocessing()
    gdf_widths, _ = processor.processing(dct_cfg_V4)
    print(gdf_widths)
    gdf_widths.to_file("widths_example4.shp")

    print("")
    print("===== BASProcessing Example #4 = END =====")


def example_5():
    """Example_5 :
        Watermask cleaning with reference waterbodies + watermask labelling
        Sections NOT available
        2 width products over the same mask
    """

    print("===== BASProcessing Example #5 = BEGIN =====")
    print("")

    # Load reaches
    gdf_reaches_cplx = gpd.read_file(shp_reaches_cplx)

    # Compute sections
    dct_geom_attr = {"reaches": {"reaches_id": "reach_id"},
                     "nodes": {"reaches_id": "reach_id",
                               "nodes_id": "node_id",
                               "pwidth": "p_width",
                               "pwse": "p_wse"}}
    obj_rivergeom = RiverGeomProduct.from_shp(reaches_shp=shp_reaches_cplx,
                                              nodes_shp=shp_nodes_cplx,
                                              bool_edge=False,
                                              dct_attr=dct_geom_attr)
    obj_rivergeom.draw_allreaches_centerline()
    gdf_sections_ortho = obj_rivergeom.draw_allreaches_sections(type="ortho")
    print(gdf_sections_ortho)
    print("")

    # Set configs #5
    dct_cfg_V5 = {"clean": {"bool_clean": True,
                            "type_clean": "waterbodies",
                            "fpath_wrkdir": "/home/charlotte/Work/AT-SWOT/cal-val/git/BAS/examples",
                            "gdf_waterbodies": gdf_waterbodies
                            },
                  "label": {"bool_label": True,
                            "type_label": "base",
                            "fpath_wrkdir": "/home/charlotte/Work/AT-SWOT/cal-val/git/BAS/examples"
                            },
                  "widths": {"scenario": 11
                             }
                  }

    # Instanciate basprocessor(s)
    processor_a = BASProcessor(
        str_watermask_tif=watermask_tif,
        gdf_sections=gdf_sections_ortho,
        gdf_reaches=gdf_reaches_cplx,
        attr_reachid="reach_id",
        str_proj="proj",
        str_provider="EO"
    )
    processor_a.preprocessing()
    gdf_widths_a, _ = processor_a.processing(dct_cfg_V5)
    gdf_widths_a["reach_id"] = gdf_widths_a["reach_id"].astype(str)
    gdf_widths_a["node_id"] = gdf_widths_a["node_id"].astype(int).astype(str)
    gdf_widths_a.to_file("widths_example5.shp")

    print("")
    print("===== BASProcessing Example #5 = END =====")


def example_6():
    """Example_6 :
        Watermask cleaning with reference waterbodies + watermask labelling
        Sections NOT available
        2 width products over the same mask
    """

    print("===== BASProcessing Example #6 = BEGIN =====")
    print("")

    # Load reaches
    gdf_reaches_cplx = gpd.read_file(shp_reaches_cplx)

    # Compute sections
    dct_geom_attr = {"reaches": {"reaches_id": "reach_id"},
                     "nodes": {"reaches_id": "reach_id",
                               "nodes_id": "node_id",
                               "pwidth": "p_width",
                               "pwse": "p_wse"}}
    obj_rivergeom = RiverGeomProduct.from_shp(reaches_shp=shp_reaches_cplx,
                                              nodes_shp=shp_nodes_cplx,
                                              bool_edge=False,
                                              dct_attr=dct_geom_attr)
    obj_rivergeom.draw_allreaches_centerline()
    gdf_sections_ortho = obj_rivergeom.draw_allreaches_sections(type="ortho")
    gdf_sections_chck = obj_rivergeom.draw_allreaches_sections(type="chck")

    # Set configs #6
    dct_cfg_V6a = {"clean": {"bool_clean": True,
                             "type_clean": "waterbodies",
                             "fpath_wrkdir": "/home/charlotte/Work/AT-SWOT/cal-val/git/BAS/examples",
                             "gdf_waterbodies": gdf_waterbodies
                             },
                   "label": {"bool_label": True,
                             "type_label": "base",
                             "fpath_wrkdir": "/home/charlotte/Work/AT-SWOT/cal-val/git/BAS/examples"
                             },
                   "widths": {"scenario": 11
                              }
                   }

    dct_cfg_V6b = {"clean": {"bool_clean": False,
                             "type_clean": "waterbodies",
                             "fpath_wrkdir": "/home/charlotte/Work/AT-SWOT/cal-val/git/BAS/examples",
                             "gdf_waterbodies": gdf_waterbodies
                             },
                   "label": {"bool_label": False,
                             "type_label": "base",
                             "fpath_wrkdir": "/home/charlotte/Work/AT-SWOT/cal-val/git/BAS/examples"
                             },
                   "widths": {"scenario": 0
                              }
                   }

    # Instanciate basprocessor(s)
    processor_a = BASProcessor(
        str_watermask_tif=watermask_tif,
        gdf_sections=gdf_sections_ortho,
        gdf_reaches=gdf_reaches_cplx,
        attr_reachid="reach_id",
        str_proj="proj",
        str_provider="EO"
    )
    processor_a.preprocessing()
    gdf_widths_a, str_fpath_updated_wm_tif = processor_a.processing(dct_cfg_V6a)
    gdf_widths_a["reach_id"] = gdf_widths_a["reach_id"].astype(str)
    gdf_widths_a["node_id"] = gdf_widths_a["node_id"].astype(int).astype(str)

    processor_b = BASProcessor(
        str_watermask_tif=str_fpath_updated_wm_tif,
        gdf_sections=gdf_sections_chck,
        gdf_reaches=gdf_reaches_cplx,
        attr_reachid="reach_id",
        str_proj="proj",
        str_provider="EO"
    )

    processor_b.preprocessing()
    processor_b.watermask.bool_labelled = True

    gdf_widths_b, _ = processor_b.processing(dct_cfg_V6b)
    gdf_widths_b["reach_id"] = gdf_widths_b["reach_id"].astype(str)
    gdf_widths_b["node_id"] = gdf_widths_b["node_id"].astype(int).astype(str)
    gdf_widths_b.to_file("widths_example6.shp")

    print("")
    print("===== BASProcessing Example #6 = END =====")


if __name__ == "__main__":
    # Load reference waterbodies - cfg 4-5
    gdf_waterbodies = gpd.read_file(ref_watermask_tif)

    # Load sections and reaches - cfg 1-4
    gdf_reaches = gpd.read_file(shp_reaches_smpl)
    gdf_sections = gpd.read_file(shp_sections_smpl)
    gdf_sections.rename(mapper={"segment": "id"}, inplace=True, axis=1)

    # Run example 1
    #example_1()

    # Run example 2
    #example_2()

    # Run example 3
    #example_3()

    # Run example 4
    #example_4()

    # # Run example 5
    # example_5()

    # Run example 6
    example_6()
