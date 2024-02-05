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

import os

os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd

from basprocessor import WidthProcessing

# Inputs for test_1
str_watermask_1_tif = "unit_watermask_1.tif"
str_reaches_1_shp = "unit_centerline_1.shp"
str_sections_1_shp = "unit_sections_1.shp"

# Inputs for test_2
str_watermask_2_tif = "unit_watermask_2.tif"
str_reaches_2_shp = "unit_centerline_2.shp"
str_sections_2_shp = "unit_sections_2.shp"

# Inputs for test_3
str_watermask_3_tif = "unit_watermask_3.tif"
str_reaches_3_shp = "unit_centerline_3.shp"
str_sections_3_shp = "unit_sections_3.shp"


def test_1():
    # Load sections and reaches
    gdf_reaches = gpd.read_file(str_reaches_1_shp)
    gdf_sections = gpd.read_file(str_sections_1_shp)

    # Set config #1
    dct_cfg_test1 = {"clean": {"bool_clean": False,
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

    # Instanciate basprocessor
    processor = WidthProcessing(
        str_watermask_tif=str_watermask_1_tif,
        gdf_sections=gdf_sections,
        gdf_reaches=gdf_reaches,
        attr_reachid="id",
        str_proj="proj"
    )
    processor.preprocessing()
    gdf_widths = processor.processing(dct_cfg_test1)
    print(gdf_widths)


def test_2():
    # Load sections and reaches
    gdf_reaches = gpd.read_file(str_reaches_2_shp)
    gdf_sections = gpd.read_file(str_sections_2_shp)

    # Set config #2
    dct_cfg_test2 = {"clean": {"bool_clean": False,
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

    # Instanciate basprocessor
    processor = WidthProcessing(
        str_watermask_tif=str_watermask_2_tif,
        gdf_sections=gdf_sections,
        gdf_reaches=gdf_reaches,
        attr_reachid="id",
        str_proj="proj"
    )
    processor.preprocessing()
    gdf_widths = processor.processing(dct_cfg_test2)
    print(gdf_widths)


def test_3():
    # Load sections and reaches
    gdf_reaches = gpd.read_file(str_reaches_3_shp)
    gdf_sections = gpd.read_file(str_sections_3_shp)

    # Set config #3
    dct_cfg_test3 = {"clean": {"bool_clean": False,
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

    # Instanciate basprocessor
    processor = WidthProcessing(
        str_watermask_tif=str_watermask_3_tif,
        gdf_sections=gdf_sections,
        gdf_reaches=gdf_reaches,
        attr_reachid="id",
        str_proj="proj"
    )
    processor.preprocessing()
    gdf_widths = processor.processing(dct_cfg_test3)
    print(gdf_widths)


if __name__ == "__main__":

    # Test 1
    test_1()

    # Test 2
    test_2()

    # Test 3
    test_3()
