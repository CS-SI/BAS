import os

os.environ['USE_PYGEOS'] = '0'
from argparse import ArgumentParser
import geopandas as gpd
import tarfile
from matplotlib import pyplot as plt

tar_gold_results = "widths_examples_gold.tar.gz"

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Parse and process data from ICESat-2"
    )

    parser.add_argument(
        "-id_test", dest="id_test", default=None,
        help="index(es) of tests to run")
    parser.add_argument(
        "-output_run", dest="output_run", default=None,
        help="Path(es) to shp files with width outputs to tests")

    # Read and check optionnal input arguments
    opt = parser.parse_args()
    if opt.id_test is None:
        l_id_test = [1, 2, 3, 4, 5, 6]

    else:
        l_id_test = [int(i) for i in opt.id_test.split(",")]
    print(l_id_test)
    if opt.output_run is None:
        raise ValueError("Missing widths file to validated")
    else:
        l_output_run = [s for s in opt.output_run.split(",")]
        if len(l_output_run) != len(l_id_test):
            raise ValueError("Number of input file different from number of input test id")

    # Loop over tests to realize
    for id_test in l_id_test:

        print("Check test {} ...".format(id_test))

        # Untar gold simulations
        str_gold = "widths_example{}.shp".format(id_test)
        tar = tarfile.open(tar_gold_results)
        for ext in ["shp", "shx", "prj", "dbf", "cpg"]:
            tar.extract("widths_example{}.{}".format(id_test, ext))
        tar.close()

        # Load gold width
        gdf_gold = gpd.read_file(str_gold)
        print(gdf_gold)

        # Load test width
        gdf_test = gpd.read_file(l_output_run[id_test-1])
        print(gdf_test)

        # merge width
        gdf_test.rename(mapper={"width": "width_ref"}, axis=1, inplace=True)
        gdf_test["width_gold"] = gdf_gold["width"]
        gdf_test["width_diff"] = gdf_test["width_ref"] - gdf_test["width_gold"]

        gdf_test.plot.scatter(x="width_ref", y="width_gold", c="width_diff", colormap="viridis")
        print("Correlation : {}".format(gdf_test["width_gold"].corr(gdf_gold["width"])))
        print(gdf_test["width_diff"].describe())

        if abs(gdf_test["width_diff"].mean()) > 0.1:
            print("Test failed ...")
        else:
            print("Test passed ...")
        plt.show()
