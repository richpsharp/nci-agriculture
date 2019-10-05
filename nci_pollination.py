"""
Pollination analysis for NCI project. This is based off the IPBES-Pollination
project so that Peter can run specific landcover maps with given price data.
"""
import argparse
import re
import collections
import csv
import shutil
import glob
import sys
import zipfile
import time
import os
import logging
import itertools
import tempfile

import ecoshard
import rtree
import shapely.wkb
import shapely.prepared
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import pandas
import numpy
import scipy.ndimage.morphology
import taskgraph
import pygeoprocessing

# format of the key pairs is [data suffix]: [landcover raster]
# these must be ESA landcover map type
LANDCOVER_DATA_MAP = {
    'data_suffix': 'landcover raster.tif',
}

# set a 1GB limit for the cache
gdal.SetCacheMax(2**30)

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger('nci_pollination')
logging.getLogger('taskgraph').setLevel(logging.ERROR)

_MULT_NODATA = -1
# the following are the globio landcover codes. A tuple (x, y) indicates the
# inclusive range from x to y. Pollinator habitat was defined as any natural
# land covers, as defined (GLOBIO land-cover classes 6, secondary vegetation,
# and  50-180, various types of primary vegetation). To test sensitivity to
# this definition we included "semi-natural" habitats (GLOBIO land-cover
# classes 3, 4, and 5; pasture, rangeland and forestry, respectively) in
# addition to "natural", and repeated all analyses with semi-natural  plus
# natural habitats, but this did not substantially alter the results  so we do
# not include it in our final analysis or code base.

GLOBIO_AG_CODES = [2, (10, 40), (230, 232)]
GLOBIO_NATURAL_CODES = [6, (50, 180)]

WORKING_DIR = './nci_pollination_workspace'
ECOSHARD_DIR = os.path.join(WORKING_DIR, 'ecoshard_dir')
CHURN_DIR = os.path.join(WORKING_DIR, 'churn')

NODATA = -9999
N_WORKERS = -1 #max(1, multiprocessing.cpu_count())

CROP_PRICES_URL = 'https://storage.googleapis.com/nci-ecoshards/prices_by_crop_and_country_md5_af6233d592d4a01d8413f50c8ccbc78d.csv'
COUNTRY_ISO_GPKG_URL = 'https://storage.googleapis.com/nci-ecoshards/country_shapefile-20191004T192454Z-001_md5_4eb621c6c74707f038d9ac86a4f2b662.gpkg'
YIELD_AND_HAREA_ZIP_URL = 'https://storage.googleapis.com/nci-ecoshards/ipbes_monfreda_2008_observed_yield_and_harea_md5_49b529f57071cc85abbd02b6e105089b.zip'


def calculate_for_landcover(landcover_path):
    """Entry point."""
    landcover_key = os.path.splitext(os.path.basename(landcover_path))[0]
    output_dir = os.path.join(WORKING_DIR, landcover_key)
    for dir_path in [output_dir, ECOSHARD_DIR, CHURN_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    task_graph = taskgraph.TaskGraph(
        CHURN_DIR, N_WORKERS, reporting_interval=5.0)

    # CROP dependent prices?
    crop_prices_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(CROP_PRICES_URL))
    crop_prices_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(CROP_PRICES_URL, crop_prices_path),
        target_path_list=[crop_prices_path],
        task_name='download crop prices')

    country_iso_gpkg_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(COUNTRY_ISO_GPKG_URL))
    country_iso_gpkg_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(COUNTRY_ISO_GPKG_URL, country_iso_gpkg_path),
        target_path_list=[country_iso_gpkg_path],
        task_name='download country iso gpkg')

    zip_touch_file_path = os.path.join(
        ECOSHARD_DIR, 'monfreda_2008_observed_yield_and_harea.txt')
    yield_and_harea_task = task_graph.add_task(
        func=download_and_unzip,
        args=(YIELD_AND_HAREA_ZIP_URL, ECOSHARD_DIR,
              zip_touch_file_path),
        target_path_list=[zip_touch_file_path],
        task_name='download and unzip yield and harea')

    crop_prices_task.join()

    country_crop_price_map = collections.defaultdict(
        lambda: collections.defaultdict(dict))
    LOGGER.debug('parse crop prices table')
    with open(crop_prices_path, 'r') as crop_prices_file:
        csv_reader = csv.DictReader(crop_prices_file)
        for row in csv_reader:
            price_list = [row[year] for year in [
                '2014', '2013', '2012', '2011', '2010'] if row[year] != '']
            if price_list:
                country_crop_price_map[
                    row['ISO3']][row['earthstat_filename_prefix']] = float(
                        price_list[0])

    yield_and_harea_raster_dir = os.path.join(
        ECOSHARD_DIR, 'monfreda_2008_observed_yield_and_harea')

    for yield_raster_path in glob.glob(os.path.join(
            yield_and_harea_raster_dir, '*_yield.tif')):
        crop_name = re.match(
            '([^_]+)_.*\.tif', os.path.basename(yield_raster_path))[1]
        LOGGER.debug(crop_name)
        crop_price_raster_path = os.path.join(
            CHURN_DIR, '%s_price.tif' % crop_name)
        price_raster_task = task_graph.add_task(
            func=create_price_raster,
            args=(
                yield_raster_path, country_iso_gpkg_path,
                country_crop_price_map, crop_name, crop_price_raster_path),
            dependent_task_list=[yield_and_harea_task, country_iso_gpkg_task],
            ignore_path_list=[country_iso_gpkg_path],
            target_path_list=[crop_price_raster_path],
            task_name='%s price raster' % crop_name)

    # Crop content of critical macro and micronutrients (KJ energy/100 g, IU
    #   Vitamin A/ 100 g and mcg Folate/100g) for the 115 crops were taken
    #   from USDA (2011) . The USDA (2011) data also provided estimated refuse
    #   of the food item (e.g., peels, seeds). The pollination-dependent yield
    #   was reduced by this refuse percent and then multiplied by nutrient
    #   content, and summed across all crops to derive pollination-dependent
    #   nutrient yields (KJ/ha, IU Vitamin A/ha, mcg Folate/ha) for each
    #   nutrient at 5 arc min. The full table used in this analysis can be
    # found at https://storage.googleapis.com/ecoshard-root/'
    # 'crop_nutrient_md5_d6e67fd79ef95ab2dd44ca3432e9bb4d.csv
    crop_nutrient_url = (
        'https://storage.googleapis.com/nci-ecoshards/'
        'crop_nutrient_md5_2fbe7455357f8008a12827fd88816fc1.csv')
    crop_nutrient_table_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(crop_nutrient_url))

    crop_nutrient_table_fetch_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(
            crop_nutrient_url, crop_nutrient_table_path),
        target_path_list=[crop_nutrient_table_path],
        task_name=f'fetch {os.path.basename(crop_nutrient_table_path)}')


    ######### DO crop value here
    target_10km_value_yield_path = os.path.join(
        CHURN_DIR, 'monfreda_2008_value_yield_rasters',
        f'monfreda_2008_value_yield_total_10km.tif')
    target_10s_value_yield_path = os.path.join(
        CHURN_DIR, 'monfreda_2008_value_yield_rasters',
        f'monfreda_2008_value_yield_total_10s.tif')
    target_10s_value_path = os.path.join(
        CHURN_DIR, 'monfreda_2008_value_rasters',
        f'monfreda_2008_value_total_10s.tif')
    value_total_task = task_graph.add_task(
        func=create_value_rasters,
        args=(
             crop_nutrient_table_path,
             yield_and_harea_raster_dir, CHURN_DIR, False, landcover_path,
             target_10km_value_yield_path, target_10s_value_yield_path,
             target_10s_value_path),
        target_path_list=[
            target_10km_value_yield_path, target_10s_value_yield_path,
            target_10s_value_path],
        dependent_task_list=[
            yield_and_harea_task,
            crop_nutrient_table_fetch_task],
        task_name=f"""create prod raster {
            os.path.basename(target_10s_value_path)}""")

    task_graph.join()
    task_graph.close()
    return

    ########## Now pollination



    # 1.2.3.  Crop production

    # Spatially-explicit global crop yields (tons/ha) at 5 arc min (~10 km)
    # were taken from Monfreda et al. (2008) for 115 crops (permanent link to
    # crop yield folder). These yields were multiplied by crop pollination
    # dependency to calculate the pollination-dependent crop yield for each 5
    # min grid cell. Note the monfreda maps are in units of per-hectare
    # yields

    prod_total_nut_10s_task_path_map = {}
    poll_dep_prod_nut_10s_task_path_map = {}
    for nut_id, nutrient_name in [
            ('en', 'Energy'), ('va', 'VitA'), ('fo', 'Folate')]:
        # total annual production of nutrient
        yield_total_nut_10km_path = os.path.join(
            CHURN_DIR, 'monfreda_2008_yield_nutrient_rasters',
            f'monfreda_2008_yield_total_{nut_id}_10km.tif')
        yield_total_nut_10s_path = os.path.join(
            CHURN_DIR, 'monfreda_2008_yield_nutrient_rasters',
            f'monfreda_2008_yield_total_{nut_id}_10s.tif')
        prod_total_nut_10s_path = os.path.join(
            CHURN_DIR, 'monfreda_2008_prod_nutrient_rasters',
            f'monfreda_2008_prod_total_{nut_id}_10s.tif')

        prod_total_task = task_graph.add_task(
            func=create_prod_nutrient_raster,
            args=(
                 crop_nutrient_table_path, nutrient_name,
                 yield_and_harea_raster_dir, False, landcover_path,
                 yield_total_nut_10km_path, yield_total_nut_10s_path,
                 prod_total_nut_10s_path),
            target_path_list=[
                yield_total_nut_10km_path, yield_total_nut_10s_path,
                prod_total_nut_10s_path],
            dependent_task_list=[
                yield_and_harea_task,
                crop_nutrient_table_fetch_task],
            task_name=f"""create prod raster {
                os.path.basename(prod_total_nut_10s_path)}""")
        prod_total_nut_10s_task_path_map[nut_id] = (
            prod_total_task, prod_total_nut_10s_path)

        # pollination-dependent annual production of nutrient
        poll_dep_yield_nut_10km_path = os.path.join(
            CHURN_DIR, 'monfreda_2008_yield_poll_dep_rasters',
            f'monfreda_2008_yield_poll_dep_{nut_id}_10km.tif')
        poll_dep_yield_nut_10s_path = os.path.join(
            CHURN_DIR, 'monfreda_2008_yield_poll_dep_rasters',
            f'monfreda_2008_yield_poll_dep_{nut_id}_10s.tif')
        poll_dep_prod_nut_10s_path = os.path.join(
            CHURN_DIR, 'monfreda_2008_prod_poll_dep_rasters',
            f'monfreda_2008_prod_poll_dep_{nut_id}_10s.tif')
        pol_dep_prod_task = task_graph.add_task(
            func=create_prod_nutrient_raster,
            args=(
                crop_nutrient_table_path, nutrient_name,
                yield_and_harea_raster_dir, True, landcover_path,
                poll_dep_yield_nut_10km_path, poll_dep_yield_nut_10s_path,
                poll_dep_prod_nut_10s_path),
            target_path_list=[
                poll_dep_yield_nut_10km_path, poll_dep_yield_nut_10s_path,
                poll_dep_prod_nut_10s_path],
            dependent_task_list=[yield_and_harea_task],
            task_name=f"""create poll dep production raster {
                os.path.basename(poll_dep_prod_nut_10s_path)}""")
        poll_dep_prod_nut_10s_task_path_map[nut_id] = (
            pol_dep_prod_task, poll_dep_prod_nut_10s_path)

    # The proportional area of natural within 2 km was calculated for every
    #  pixel of agricultural land (GLOBIO land-cover classes 2, 230, 231, and
    #  232) at 10 arc seconds (~300 m) resolution. This 2 km scale represents
    #  the distance most commonly found to be predictive of pollination
    #  services (Kennedy et al. 2013).
    kernel_raster_path = os.path.join(CHURN_DIR, 'radial_kernel.tif')
    kernel_task = task_graph.add_task(
        func=create_radial_convolution_mask,
        args=(0.00277778, 2000., kernel_raster_path),
        target_path_list=[kernel_raster_path],
        task_name='make convolution kernel')

    prod_poll_dep_realized_1d_task_path_map = {}
    prod_poll_dep_unrealized_1d_task_path_map = {}
    prod_total_realized_1d_task_path_map = {}
    # mask landcover into agriculture and pollinator habitat

    # This loop is so we don't duplicate code for 'ag' and 'hab' with the
    # only difference being the lulc codes and prefix
    for mask_prefix, globio_codes in [
            ('ag', GLOBIO_AG_CODES), ('hab', GLOBIO_NATURAL_CODES)]:
        mask_key = f'{landcover_key}_{mask_prefix}_mask'
        mask_target_path = os.path.join(
            CHURN_DIR, f'{mask_prefix}_mask', f'{mask_key}.tif')

        mask_task = task_graph.add_task(
            func=mask_raster,
            args=(landcover_path, globio_codes, mask_target_path),
            target_path_list=[mask_target_path],
            task_name=f'mask {mask_key}',)

        if mask_prefix == 'hab':
            hab_task_path_tuple = (mask_task, mask_target_path)
        elif mask_prefix == 'ag':
            ag_task_path_tuple = (mask_task, mask_target_path)

    pollhab_2km_prop_path = os.path.join(
        CHURN_DIR, 'pollhab_2km_prop',
        f'pollhab_2km_prop_{landcover_key}.tif')
    pollhab_2km_prop_task = task_graph.add_task(
        func=pygeoprocessing.convolve_2d,
        args=[
            (hab_task_path_tuple[1], 1), (kernel_raster_path, 1),
            pollhab_2km_prop_path],
        kwargs={
            'working_dir': CHURN_DIR,
            'ignore_nodata': True,
            'gtiff_creation_options': (
                'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
                'PREDICTOR=3', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
                'NUM_THREADS=2'),
            'n_threads': 4},
        dependent_task_list=[hab_task_path_tuple[0], kernel_task],
        target_path_list=[pollhab_2km_prop_path],
        task_name=(
            'calculate proportional'
            f' {os.path.basename(pollhab_2km_prop_path)}'))

    # calculate pollhab_2km_prop_on_ag_10s by multiplying pollhab_2km_prop
    # by the ag mask
    pollhab_2km_prop_on_ag_path = os.path.join(
        output_dir, f'''pollhab_2km_prop_on_ag_10s_{
            landcover_key}.tif''')
    pollhab_2km_prop_on_ag_task = task_graph.add_task(
        func=mult_rasters,
        args=(
            ag_task_path_tuple[1], pollhab_2km_prop_path,
            pollhab_2km_prop_on_ag_path),
        target_path_list=[pollhab_2km_prop_on_ag_path],
        dependent_task_list=[
            pollhab_2km_prop_task, ag_task_path_tuple[0]],
        task_name=(
            f'''pollhab 2km prop on ag {
                os.path.basename(pollhab_2km_prop_on_ag_path)}'''))

    #  1.1.4.  Sufficiency threshold A threshold of 0.3 was set to
    #  evaluate whether there was sufficient pollinator habitat in the 2
    #  km around farmland to provide pollination services, based on
    #  Kremen et al.'s (2005)  estimate of the area requirements for
    #  achieving full pollination. This produced a map of wild
    #  pollination sufficiency where every agricultural pixel was
    #  designated in a binary fashion: 0 if proportional area of habitat
    #  was less than 0.3; 1 if greater than 0.3. Maps of pollination
    #  sufficiency can be found at (permanent link to output), outputs
    #  "poll_suff_..." below.

    threshold_val = 0.3
    pollinator_suff_hab_path = os.path.join(
        CHURN_DIR, 'poll_suff_hab_ag_coverage_rasters',
        f'poll_suff_ag_coverage_prop_10s_{landcover_key}.tif')
    poll_suff_task = task_graph.add_task(
        func=threshold_select_raster,
        args=(
            pollhab_2km_prop_path,
            ag_task_path_tuple[1], threshold_val,
            pollinator_suff_hab_path),
        target_path_list=[pollinator_suff_hab_path],
        dependent_task_list=[
            pollhab_2km_prop_task, ag_task_path_tuple[0]],
        task_name=f"""poll_suff_ag_coverage_prop {
            os.path.basename(pollinator_suff_hab_path)}""")

    # tot_prod_en|va|fo_10s|1d_cur|ssp1|ssp3|ssp5
    # total annual production of energy (KJ/yr), vitamin A (IU/yr),
    # and folate (mg/yr)
    for nut_id in ('en', 'va', 'fo'):
        tot_prod_task, tot_prod_path = (
            prod_total_nut_10s_task_path_map[nut_id])

        prod_total_potential_path = os.path.join(
            output_dir, f'''prod_total_potential_{
                nut_id}_10s_{landcover_key}.tif''')
        prod_total_potential_task = task_graph.add_task(
            func=mult_rasters,
            args=(
                ag_task_path_tuple[1], tot_prod_path,
                prod_total_potential_path),
            target_path_list=[prod_total_potential_path],
            dependent_task_list=[tot_prod_task, ag_task_path_tuple[0]],
            task_name=(
                f'tot_prod_{nut_id}_10s_{landcover_key}'))
        schedule_aggregate_to_degree(
            task_graph, prod_total_potential_path, numpy.sum,
            prod_total_potential_task)

        poll_dep_prod_task, poll_dep_prod_path = (
            poll_dep_prod_nut_10s_task_path_map[nut_id])

        prod_poll_dep_potential_nut_scenario_path = os.path.join(
            output_dir,
            f'prod_poll_dep_potential_{nut_id}_10s_'
            f'{landcover_key}.tif')
        prod_poll_dep_potential_nut_scenario_task = task_graph.add_task(
            func=mult_rasters,
            args=(
                ag_task_path_tuple[1], poll_dep_prod_path,
                prod_poll_dep_potential_nut_scenario_path),
            target_path_list=[prod_poll_dep_potential_nut_scenario_path],
            dependent_task_list=[
                poll_dep_prod_task, ag_task_path_tuple[0]],
            task_name=(
                f'poll_dep_prod_{nut_id}_'
                f'10s_{landcover_key}'))
        (prod_poll_dep_potential_nut_scenario_1d_task,
         prod_poll_dep_potential_nut_scenario_1d_path) = (
            schedule_aggregate_to_degree(
                task_graph, prod_poll_dep_potential_nut_scenario_path,
                numpy.sum, prod_poll_dep_potential_nut_scenario_task))

        # pollination independent
        prod_poll_indep_nut_scenario_path = os.path.join(
            output_dir,
            f'prod_poll_indep_{nut_id}_10s_'
            f'{landcover_key}.tif')
        prod_poll_indep_nut_scenario_task = task_graph.add_task(
            func=subtract_2_rasters,
            args=(
                prod_total_potential_path,
                prod_poll_dep_potential_nut_scenario_path,
                prod_poll_indep_nut_scenario_path),
            target_path_list=[prod_poll_indep_nut_scenario_path],
            dependent_task_list=[
                prod_total_potential_task,
                prod_poll_dep_potential_nut_scenario_task],
            task_name=(
                f'prod_poll_indep_{nut_id}_'
                f'10s_{landcover_key}'))
        schedule_aggregate_to_degree(
            task_graph, prod_poll_indep_nut_scenario_path,
            numpy.sum, prod_poll_indep_nut_scenario_task)

        # prod_poll_dep_realized_en|va|fo_10s|1d_cur|ssp1|ssp3|ssp5:
        # pollination-dependent annual production of energy (KJ/yr),
        # vitamin A (IU/yr), and folate (mg/yr) that can be met by wild
        # pollinators due to the proximity of sufficient habitat.
        prod_poll_dep_realized_nut_scenario_path = os.path.join(
            output_dir,
            f'prod_poll_dep_realized_{nut_id}_10s_'
            f'{landcover_key}.tif')
        prod_poll_dep_realized_nut_scenario_task = task_graph.add_task(
            func=mult_rasters,
            args=(
                prod_poll_dep_potential_nut_scenario_path,
                pollinator_suff_hab_path,
                prod_poll_dep_realized_nut_scenario_path),
            target_path_list=[prod_poll_dep_realized_nut_scenario_path],
            dependent_task_list=[
                poll_suff_task, prod_poll_dep_potential_nut_scenario_task],
            task_name=(
                f'prod_poll_dep_realized_{nut_id}_'
                f'10s_{landcover_key}'))

        (prod_poll_dep_realized_1d_task,
         prod_poll_dep_realized_nut_scenario_1d_path) = (
            schedule_aggregate_to_degree(
                task_graph, prod_poll_dep_realized_nut_scenario_path,
                numpy.sum, prod_poll_dep_realized_nut_scenario_task))
        prod_poll_dep_realized_1d_task_path_map[
            (landcover_key, nut_id)] = (
                prod_poll_dep_realized_1d_task,
                prod_poll_dep_realized_nut_scenario_1d_path)

        # calculate prod_poll_dep_unrealized X1 as
        # prod_total - prod_poll_dep_realized
        prod_poll_dep_unrealized_nut_scenario_path = os.path.join(
            output_dir,
            f'prod_poll_dep_unrealized_{nut_id}_10s_'
            f'{landcover_key}.tif')
        prod_poll_dep_unrealized_nut_scenario_task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=([
                (prod_poll_dep_potential_nut_scenario_path, 1),
                (prod_poll_dep_realized_nut_scenario_path, 1),
                (_MULT_NODATA, 'raw'),
                (_MULT_NODATA, 'raw'),
                (_MULT_NODATA, 'raw')], sub_two_op,
                prod_poll_dep_unrealized_nut_scenario_path,
                gdal.GDT_Float32, -1),
            target_path_list=[prod_poll_dep_unrealized_nut_scenario_path],
            dependent_task_list=[
                prod_poll_dep_realized_nut_scenario_task,
                prod_poll_dep_potential_nut_scenario_task],
            task_name=f'''prod poll dep unrealized: {
                os.path.basename(
                    prod_poll_dep_unrealized_nut_scenario_path)}''')
        (prod_poll_dep_unrealized_nut_scenario_1d_task,
         prod_poll_dep_unrealized_nut_scenario_1d_path) = (
            schedule_aggregate_to_degree(
                task_graph, prod_poll_dep_unrealized_nut_scenario_path,
                numpy.sum, prod_poll_dep_unrealized_nut_scenario_task))
        prod_poll_dep_unrealized_1d_task_path_map[
            (landcover_key, nut_id)] = (
                prod_poll_dep_unrealized_nut_scenario_1d_task,
                prod_poll_dep_unrealized_nut_scenario_1d_path)

        # calculate prod_total_realized as
        #   prod_total_potential - prod_poll_dep_unrealized
        prod_total_realized_nut_scenario_path = os.path.join(
            output_dir,
            f'prod_total_realized_{nut_id}_10s_'
            f'{landcover_key}.tif')
        prod_total_realized_nut_scenario_task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=([
                (prod_total_potential_path, 1),
                (prod_poll_dep_unrealized_nut_scenario_path, 1),
                (_MULT_NODATA, 'raw'),
                (_MULT_NODATA, 'raw'),
                (_MULT_NODATA, 'raw')], sub_two_op,
                prod_total_realized_nut_scenario_path,
                gdal.GDT_Float32, _MULT_NODATA),
            target_path_list=[prod_total_realized_nut_scenario_path],
            dependent_task_list=[
                prod_poll_dep_unrealized_nut_scenario_task,
                prod_total_potential_task],
            task_name=f'''prod poll dep unrealized: {
                os.path.basename(
                    prod_total_realized_nut_scenario_path)}''')
        (prod_total_realized_nut_1d_scenario_task,
         prod_total_realized_nut_1d_scenario_path) = (
            schedule_aggregate_to_degree(
                task_graph, prod_total_realized_nut_scenario_path,
                numpy.sum, prod_total_realized_nut_scenario_task))
        prod_total_realized_1d_task_path_map[
            (landcover_key, nut_id)] = (
                prod_total_realized_nut_1d_scenario_task,
                prod_total_realized_nut_1d_scenario_path)

    task_graph.close()
    task_graph.join()


def build_spatial_index(vector_path):
    """Build an rtree/geom list tuple from ``vector_path``."""
    vector = gdal.OpenEx(vector_path)
    layer = vector.GetLayer()
    geom_index = rtree.index.Index()
    geom_list = []
    for index in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(index)
        geom = feature.GetGeometryRef()
        shapely_geom = shapely.wkb.loads(geom.ExportToWkb())
        shapely_prep_geom = shapely.prepared.prep(shapely_geom)
        geom_list.append(shapely_prep_geom)
        geom_index.insert(index, shapely_geom.bounds)

    return geom_index, geom_list


def calculate_total_requirements(
        pop_path_list, nut_need_list, target_path):
    """Calculate total nutrient requirements.

    Create a new raster by summing all rasters in `pop_path_list` multiplied
    by their corresponding scalar in `nut_need_list`.

    Parameters:
        pop_path_list (list of str): list of paths to population counts.
        nut_need_list (list): list of scalars that correspond in order to
            the per-count nutrient needs of `pop_path_list`.
        target_path (str): path to target file.

    Return:
        None.

    """
    nodata = -1
    pop_nodata = pygeoprocessing.get_raster_info(
        pop_path_list[0])['nodata'][0]

    def mult_and_sum(*arg_list):
        """Arg list is an (array0, scalar0, array1, scalar1,...) list.

        Returns:
            array0*scalar0 + array1*scalar1 + .... but ignore nodata.

        """
        result = numpy.empty(arg_list[0].shape, dtype=numpy.float32)
        result[:] = nodata
        array_stack = numpy.array(arg_list[0::2])
        scalar_list = numpy.array(arg_list[1::2])
        # make a valid mask as big as a single array
        valid_mask = numpy.logical_and.reduce(
            array_stack != pop_nodata, axis=0)

        # mask out all invalid elements but reshape so there's still the same
        # number of arrays
        valid_array_elements = (
            array_stack[numpy.broadcast_to(valid_mask, array_stack.shape)])
        array_stack = None

        # sometimes this array is empty, check first before reshaping
        if valid_array_elements.size != 0:
            valid_array_elements = valid_array_elements.reshape(
                -1, numpy.count_nonzero(valid_mask))
            # multiply each element of the scalar with each row of the valid
            # array stack, then sum along the 0 axis to get the result
            result[valid_mask] = numpy.sum(
                (valid_array_elements.T * scalar_list).T, axis=0)
        scalar_list = None
        valid_mask = None
        valid_array_elements = None
        return result

    pygeoprocessing.raster_calculator(list(itertools.chain(*[
        ((path, 1), (scalar, 'raw')) for path, scalar in zip(
            pop_path_list, nut_need_list)])), mult_and_sum, target_path,
        gdal.GDT_Float32, nodata)


def sub_two_op(a_array, b_array, a_nodata, b_nodata, target_nodata):
    """Subtract a from b and ignore nodata."""
    result = numpy.empty_like(a_array)
    result[:] = target_nodata
    valid_mask = (a_array != a_nodata) & (b_array != b_nodata)
    result[valid_mask] = a_array[valid_mask] - b_array[valid_mask]
    return result


def average_rasters(*raster_list, clamp=None):
    """Average rasters in raster list except write to the last one.

    Parameters:
        raster_list (list of string): list of rasters to average over.
        clamp (float): value to clamp the individual raster to before the
            average.

    Returns:
        None.

    """
    nodata_list = [
        pygeoprocessing.get_raster_info(path)['nodata'][0]
        for path in raster_list[:-1]]
    target_nodata = -1.

    def average_op(*array_list):
        result = numpy.empty_like(array_list[0])
        result[:] = target_nodata
        valid_mask = numpy.ones(result.shape, dtype=numpy.bool)
        clamped_list = []
        for array, nodata in zip(array_list, nodata_list):
            valid_mask &= array != nodata
            if clamp:
                clamped_list.append(
                    numpy.where(array > clamp, clamp, array))
            else:
                clamped_list.append(array)

        if valid_mask.any():
            array_stack = numpy.stack(clamped_list)
            result[valid_mask] = numpy.average(
                array_stack[numpy.broadcast_to(
                    valid_mask, array_stack.shape)].reshape(
                        len(array_list), -1), axis=0)
        return result

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in raster_list[:-1]], average_op,
        raster_list[-1], gdal.GDT_Float32, target_nodata)


def subtract_2_rasters(
        raster_path_a, raster_path_b, target_path):
    """Calculate target = a-b and ignore nodata."""
    a_nodata = pygeoprocessing.get_raster_info(raster_path_a)['nodata'][0]
    b_nodata = pygeoprocessing.get_raster_info(raster_path_b)['nodata'][0]
    target_nodata = -9999

    def sub_op(a_array, b_array):
        """Sub a-b-c as arrays."""
        result = numpy.empty(a_array.shape, dtype=numpy.float32)
        result[:] = target_nodata
        valid_mask = (
            (a_array != a_nodata) &
            (b_array != b_nodata))
        result[valid_mask] = (
            a_array[valid_mask] - b_array[valid_mask])
        return result

    pygeoprocessing.raster_calculator(
        [(raster_path_a, 1), (raster_path_b, 1)],
        sub_op, target_path, gdal.GDT_Float32, target_nodata)


def subtract_3_rasters(
        raster_path_a, raster_path_b, raster_path_c, target_path):
    """Calculate target = a-b-c and ignore nodata."""
    a_nodata = pygeoprocessing.get_raster_info(raster_path_a)['nodata'][0]
    b_nodata = pygeoprocessing.get_raster_info(raster_path_b)['nodata'][0]
    c_nodata = pygeoprocessing.get_raster_info(raster_path_c)['nodata'][0]
    target_nodata = -9999

    def sub_op(a_array, b_array, c_array):
        """Sub a-b-c as arrays."""
        result = numpy.empty(a_array.shape, dtype=numpy.float32)
        result[:] = target_nodata
        valid_mask = (
            (a_array != a_nodata) &
            (b_array != b_nodata) &
            (c_array != c_nodata))
        result[valid_mask] = (
            a_array[valid_mask] - b_array[valid_mask] - c_array[valid_mask])
        return result

    pygeoprocessing.raster_calculator(
        [(raster_path_a, 1), (raster_path_b, 1), (raster_path_c, 1)],
        sub_op, target_path, gdal.GDT_Float32, target_nodata)


def calculate_raster_ratio(
        raster_a_path, raster_b_path, target_raster_path, clamp_to=None):
    """Calculate the ratio of a:b and ignore nodata and divide by 0.

    Parameters:
        raster_a_path (string): path to numerator of ratio
        raster_b_path (string): path to denominator of ratio
        target_raster_path (string): path to desired target raster that will
            use a nodata value of -1.
        clamp_to (numeric): if not None, the ratio is capped to this value.

    Returns:
        None.

    """
    temp_dir = tempfile.mkdtemp()
    raster_a_aligned_path = os.path.join(
        temp_dir, os.path.basename(raster_a_path))
    raster_b_aligned_path = os.path.join(
        temp_dir, os.path.basename(raster_b_path))
    target_pixel_size = pygeoprocessing.get_raster_info(
        raster_a_path)['pixel_size']
    pygeoprocessing.align_and_resize_raster_stack(
        [raster_a_path, raster_b_path],
        [raster_a_aligned_path, raster_b_aligned_path], ['near']*2,
        target_pixel_size, 'intersection')

    nodata_a = pygeoprocessing.get_raster_info(
        raster_a_aligned_path)['nodata'][0]
    nodata_b = pygeoprocessing.get_raster_info(
        raster_b_aligned_path)['nodata'][0]
    target_nodata = -1.

    def ratio_op(array_a, array_b):
        result = numpy.empty(array_a.shape, dtype=numpy.float32)
        result[:] = target_nodata
        zero_mask = numpy.isclose(array_b, 0.)
        valid_mask = (
            ~numpy.isclose(array_a, nodata_a) &
            ~numpy.isclose(array_b, nodata_b) &
            ~zero_mask)
        result[valid_mask] = array_a[valid_mask] / array_b[valid_mask]
        if clamp_to:
            result[valid_mask & (result > clamp_to)] = clamp_to
        result[zero_mask] = 0.0
        return result

    pygeoprocessing.raster_calculator(
        [(raster_a_aligned_path, 1), (raster_b_aligned_path, 1)], ratio_op,
        target_raster_path, gdal.GDT_Float32, target_nodata)

    shutil.rmtree(temp_dir, ignore_errors=True)


def warp_to_raster(
        base_raster_path, canonical_raster_path, resample_method,
        target_raster_path):
    """Warp base raster to canonical example.

    Parameters:
        base_raster_path (str): path to base raster
        canonical_raster_path (str),
        resample_method (str): one of nearest, bilinear, cubic, cubic_spline,
            lanczos, average, mode, max, min, med, q1, q3.
        target_raster_path (str): path to target warped raster.

    Returns:
        None.

    """
    canonical_raster_info = pygeoprocessing.get_raster_info(
        canonical_raster_path)
    pygeoprocessing.warp_raster(
        base_raster_path, canonical_raster_info['pixel_size'],
        target_raster_path,
        resample_method, target_bb=canonical_raster_info['bounding_box'],
        n_threads=2)


def calculate_future_pop(
        cur_ssp_path, fut_ssp_path, gpw_tot_count_path, target_ssp_pop_path):
    """Calculate future population raster.

    Multiply `gpw_tot_count` by the fut_ssp/cur_ssp ratio.
    """
    target_nodata = -1
    ssp_cur_nodata = pygeoprocessing.get_raster_info(
        cur_ssp_path)['nodata'][0]
    ssp_fut_nodata = pygeoprocessing.get_raster_info(
        fut_ssp_path)['nodata'][0]
    count_nodata = pygeoprocessing.get_raster_info(
        gpw_tot_count_path)['nodata'][0]

    def _future_pop_op(cur_array, fut_array, count_array):
        """Calculate future pop by dividing fut/cur*cur_count."""
        result = numpy.empty(cur_array.shape, dtype=numpy.float32)
        result[:] = target_nodata
        zero_mask = cur_array == 0
        valid_mask = (
            (cur_array != ssp_cur_nodata) &
            (fut_array != ssp_fut_nodata) &
            (count_array != count_nodata) & ~zero_mask)
        result[valid_mask] = (
            (fut_array[valid_mask] / cur_array[valid_mask]).astype(
                numpy.float32) * count_array[valid_mask])
        # assume if denominator is 0 we don't want to mask anything out, just
        # get it working
        result[zero_mask] = 1.0
        return result

    pygeoprocessing.raster_calculator(
        [(cur_ssp_path, 1), (fut_ssp_path, 1), (gpw_tot_count_path, 1)],
        _future_pop_op, target_ssp_pop_path, gdal.GDT_Float32, target_nodata)


def calc_pop_count(gpw_dens_path, gpw_count_path):
    """Calculate population count from density.

    Parameters:
        gpw_dens_path (string): path to density raster in units of
            people / km^2.
        gpw_count_path (string): path to target raster to generate number of
            people / pixel.

    """
    gpw_dens_info = pygeoprocessing.get_raster_info(gpw_dens_path)
    y_lat_array = numpy.linspace(
        gpw_dens_info['geotransform'][3],
        gpw_dens_info['geotransform'][3] +
        gpw_dens_info['geotransform'][5] *
        gpw_dens_info['raster_size'][1],
        gpw_dens_info['raster_size'][1])

    # `area_of_pixel` is in m^2, convert to km^2
    y_km2_array = area_of_pixel(
        abs(gpw_dens_info['geotransform'][1]),
        y_lat_array) * 1e-6
    y_km2_column = y_km2_array.reshape((y_km2_array.size, 1))

    nodata = gpw_dens_info['nodata'][0]
    try:
        os.makedirs(os.path.dirname(gpw_count_path))
    except OSError:
        pass
    pygeoprocessing.raster_calculator(
        [(gpw_dens_path, 1), y_km2_column, (nodata, 'raw')],
        density_to_value_op, gpw_count_path, gdal.GDT_Float32,
        nodata)


def create_radial_convolution_mask(
        pixel_size_degree, radius_meters, kernel_filepath):
    """Create a radial mask to sample pixels in convolution filter.

    Parameters:
        pixel_size_degree (float): size of pixel in degrees.
        radius_meters (float): desired size of radial mask in meters.

    Returns:
        A 2D numpy array that can be used in a convolution to aggregate a
        raster while accounting for partial coverage of the circle on the
        edges of the pixel.

    """
    degree_len_0 = 110574  # length at 0 degrees
    degree_len_60 = 111412  # length at 60 degrees
    pixel_size_m = pixel_size_degree * (degree_len_0 + degree_len_60) / 2.0
    pixel_radius = numpy.ceil(radius_meters / pixel_size_m)
    n_pixels = (int(pixel_radius) * 2 + 1)
    sample_pixels = 200
    mask = numpy.ones((sample_pixels * n_pixels, sample_pixels * n_pixels))
    mask[mask.shape[0]//2, mask.shape[0]//2] = 0
    distance_transform = scipy.ndimage.morphology.distance_transform_edt(mask)
    mask = None
    stratified_distance = distance_transform * pixel_size_m / sample_pixels
    distance_transform = None
    in_circle = numpy.where(stratified_distance <= 2000.0, 1.0, 0.0)
    stratified_distance = None
    reshaped = in_circle.reshape(
        in_circle.shape[0] // sample_pixels, sample_pixels,
        in_circle.shape[1] // sample_pixels, sample_pixels)
    kernel_array = numpy.sum(reshaped, axis=(1, 3)) / sample_pixels**2
    normalized_kernel_array = kernel_array / numpy.sum(kernel_array)
    reshaped = None

    driver = gdal.GetDriverByName('GTiff')
    kernel_raster = driver.Create(
        kernel_filepath.encode('utf-8'), n_pixels, n_pixels, 1,
        gdal.GDT_Float32, options=[
            'BIGTIFF=IF_SAFER', 'TILED=YES', 'BLOCKXSIZE=256',
            'BLOCKYSIZE=256'])

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_raster.SetGeoTransform([-180, 1, 0, 90, 0, -1])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    kernel_raster.SetProjection(srs.ExportToWkt())
    kernel_band = kernel_raster.GetRasterBand(1)
    kernel_band.SetNoDataValue(NODATA)
    kernel_band.WriteArray(normalized_kernel_array)


def google_bucket_fetch_and_validate(url, json_key_path, target_path):
    """Create a function to download a Google Blob to a given path.

    Parameters:
        url (string): url to blob, matches the form
            '^https://storage.cloud.google.com/([^/]*)/(.*)$'
        json_key_path (string): path to Google Cloud private key generated by
        https://cloud.google.com/iam/docs/creating-managing-service-account-keys
        target_path (string): path to target file.

    Returns:
        a function with a single `path` argument to the target file. Invoking
            this function will download the Blob to `path`.

    """
    url_matcher = re.match(
        r'^https://[^/]*\.com/([^/]*)/(.*)$', url)
    LOGGER.debug(url)
    client = google.cloud.storage.client.Client.from_service_account_json(
        json_key_path)
    bucket_id = url_matcher.group(1)
    LOGGER.debug(f'parsing bucket {bucket_id} from {url}')
    bucket = client.get_bucket(bucket_id)
    blob_id = url_matcher.group(2)
    LOGGER.debug(f'loading blob {blob_id} from {url}')
    blob = google.cloud.storage.Blob(
        blob_id, bucket, chunk_size=2**24)
    LOGGER.info(f'downloading blob {target_path} from {url}')
    try:
        os.makedirs(os.path.dirname(target_path))
    except os.error:
        pass
    blob.download_to_filename(target_path)
    if not reproduce.valid_hash(target_path, 'embedded'):
        raise ValueError(f"{target_path}' does not match its expected hash")


def threshold_select_raster(
        base_raster_path, select_raster_path, threshold_val, target_path):
    """Select `select` if `base` >= `threshold_val`.

    Parameters:
        base_raster_path (string): path to single band raster that will be
            used to determine the threshold mask to select from
            `select_raster_path`.
        select_raster_path (string): path to single band raster to pass
            through to target if aligned `base` pixel is >= `threshold_val`
            0 otherwise, or nodata if base == nodata. Must be the same
            shape as `base_raster_path`.
        threshold_val (numeric): value to use as threshold cutoff
        target_path (string): path to desired output raster, raster is a
            byte type with same dimensions and projection as
            `base_raster_path`. A pixel in this raster will be `select` if
            the corresponding pixel in `base_raster_path` is >=
            `threshold_val`, 0 otherwise or nodata if `base` == nodata.

    Returns:
        None.

    """
    base_nodata = pygeoprocessing.get_raster_info(
        base_raster_path)['nodata'][0]
    target_nodata = -9999.

    def threshold_select_op(
            base_array, select_array, threshold_val, base_nodata,
            target_nodata):
        result = numpy.empty(select_array.shape, dtype=numpy.float32)
        result[:] = target_nodata
        valid_mask = (base_array != base_nodata) & (select_array == 1)
        result[valid_mask] = numpy.interp(
            base_array[valid_mask], [0, threshold_val], [0.0, 1.0], 0, 1)
        return result

    pygeoprocessing.raster_calculator(
        [(base_raster_path, 1), (select_raster_path, 1),
         (threshold_val, 'raw'), (base_nodata, 'raw'),
         (target_nodata, 'raw')], threshold_select_op,
        target_path, gdal.GDT_Float32, target_nodata, gtiff_creation_options=(
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
            'PREDICTOR=2', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
            'NUM_THREADS=2'))


def mask_raster(base_path, codes, target_path):
    """Mask `base_path` to 1 where values are in codes. 0 otherwise.

    Parameters:
        base_path (string): path to single band integer raster.
        codes (list): list of integer or tuple integer pairs. Membership in
            `codes` or within the inclusive range of a tuple in `codes`
            is sufficient to mask the corresponding raster integer value
            in `base_path` to 1 for `target_path`.
        target_path (string): path to desired mask raster. Any corresponding
            pixels in `base_path` that also match a value or range in
            `codes` will be masked to 1 in `target_path`. All other values
            are 0.

    Returns:
        None.

    """
    code_list = numpy.array([
        item for sublist in [
            range(x[0], x[1]+1) if isinstance(x, tuple) else [x]
            for x in codes] for item in sublist])
    LOGGER.debug(f'expanded code array {code_list}')

    base_nodata = pygeoprocessing.get_raster_info(base_path)['nodata'][0]
    mask_nodata = 2

    def mask_codes_op(base_array, codes_array):
        """Return a bool raster if value in base_array is in codes_array."""
        result = numpy.empty(base_array.shape, dtype=numpy.int8)
        result[:] = mask_nodata
        valid_mask = base_array != base_nodata
        result[valid_mask] = numpy.isin(
            base_array[valid_mask], codes_array)
        return result

    pygeoprocessing.raster_calculator(
        [(base_path, 1), (code_list, 'raw')], mask_codes_op, target_path,
        gdal.GDT_Byte, 2, gtiff_creation_options=(
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
            'PREDICTOR=2', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
            'NUM_THREADS=2'))


def unzip_file(zipfile_path, target_dir, touchfile_path):
    """Unzip contents of `zipfile_path`.

    Parameters:
        zipfile_path (string): path to a zipped file.
        target_dir (string): path to extract zip file to.
        touchfile_path (string): path to a file to create if unzipping is
            successful.

    Returns:
        None.

    """
    with zipfile.ZipFile(zipfile_path, 'r') as zip_ref:
        zip_ref.extractall(target_dir)

    with open(touchfile_path, 'w') as touchfile:
        touchfile.write(f'unzipped {zipfile_path}')


def _make_logger_callback(message):
    """Build a timed logger callback that prints `message` replaced.

    Parameters:
        message (string): a string that expects a %f placement variable,
            for % complete.

    Returns:
        Function with signature:
            logger_callback(df_complete, psz_message, p_progress_arg)

    """
    def logger_callback(df_complete, psz_message, p_progress_arg):
        """Log updates using GDAL API for callbacks."""
        try:
            current_time = time.time()
            if ((current_time - logger_callback.last_time) > 5.0 or
                    (df_complete == 1.0 and
                     logger_callback.total_time >= 5.0)):
                LOGGER.info(message, df_complete * 100)
                logger_callback.last_time = current_time
                logger_callback.total_time += current_time
        except AttributeError:
            logger_callback.last_time = time.time()
            logger_callback.total_time = 0.0

    return logger_callback


def total_yield_op(
        yield_nodata, pollination_yield_factor_list,
        *crop_yield_harea_array_list):
    """Calculate total yield.

    Parameters:
        yield_nodata (numeric): nodata value for the arrays in
            ``crop_yield_array_list```.
        pollination_yield_factor_list (list of float): list of non-refuse
            proportion of yield that is pollination dependent.
        crop_yield_harea_array_list (list of numpy.ndarray): list of length
            n of 2D arrays of n/2 yield (tons/Ha) for crops that correlate in
            order with the ``pollination_yield_factor_list`` followed by
            n/2 harea (proportional area) of those crops.

    """
    result = numpy.empty(
        crop_yield_harea_array_list[0].shape, dtype=numpy.float32)
    result[:] = 0.0
    all_valid = numpy.zeros(result.shape, dtype=numpy.bool)

    n_crops = len(crop_yield_harea_array_list) // 2

    for crop_index in range(n_crops):
        crop_array = crop_yield_harea_array_list[crop_index]
        harea_array = crop_yield_harea_array_list[crop_index + n_crops]
        valid_mask = crop_array != yield_nodata
        all_valid |= valid_mask
        result[valid_mask] += (
            crop_array[valid_mask] * harea_array[valid_mask] *
            pollination_yield_factor_list[crop_index])
    result[~all_valid] = yield_nodata
    return result


def total_price_yield_op(
        yield_nodata, pollination_yield_factor_list,
        *crop_yield_harea_price_array_list):
    """Calculate total yield.

    Parameters:
        yield_nodata (numeric): nodata value for the arrays in
            ``crop_yield_array_list```.
        pollination_yield_factor_list (list of float): list of non-refuse
            proportion of yield that is pollination dependent.
        crop_yield_harea_price_array_list (list of numpy.ndarray): list of
            length n of 2D arrays of n/2 yield (tons/Ha) for crops that
            correlate in order with the ``pollination_yield_factor_list``
            followed by n/2 harea (proportional area) of those crops.

    """
    result = numpy.empty(
        crop_yield_harea_price_array_list[0].shape, dtype=numpy.float32)
    result[:] = 0.0
    all_valid = numpy.zeros(result.shape, dtype=numpy.bool)

    n_crops = len(crop_yield_harea_price_array_list) // 3
    total_area_array = numpy.zeros(result.shape)
    total_area_if_price_array = numpy.zeros(result.shape)

    for crop_index in range(n_crops):
        crop_array = crop_yield_harea_price_array_list[crop_index]
        harea_array = crop_yield_harea_price_array_list[crop_index + n_crops]
        price_array = crop_yield_harea_price_array_list[crop_index + 2*n_crops]
        valid_mask = crop_array != yield_nodata
        all_valid |= valid_mask
        total_area_array[valid_mask] += harea_array[valid_mask]
        valid_if_prices = valid_mask & (price_array != 0)
        total_area_if_price_array[valid_if_prices] += (
            harea_array[valid_if_prices])
        result[valid_mask] += (
            crop_array[valid_mask] * harea_array[valid_mask] *
            pollination_yield_factor_list[crop_index])
    scaling_mask = total_area_array != 0
    scaling_factor = total_area_array[scaling_mask] / (
        total_area_if_price_array[scaling_mask])
    result[~all_valid] = yield_nodata
    # this line rescales the prices in case there are zero prices that should
    # otherwise take up harvested area
    result[all_valid] *= scaling_factor[all_valid]
    return result


def density_to_value_op(density_array, area_array, density_nodata):
    """Calculate production.

    Parameters:
        density_array (numpy.ndarray): array of densities / area in
            ``area_array``.
        area_array (numpy.ndarray): area of each cell that corresponds with
            ``density_array``.
        density_ndoata (numeric): nodata value of the ``density_array``.

    """
    result = numpy.empty(density_array.shape, dtype=numpy.float32)
    result[:] = density_nodata
    valid_mask = density_array != density_nodata
    result[valid_mask] = density_array[valid_mask] * area_array[valid_mask]
    return result


def create_value_rasters(
        crop_nutrient_df_path, yield_and_harea_raster_dir,
        price_raster_dir, consider_pollination, sample_target_raster_path,
        target_10km_value_yield_path, target_10s_value_yield_path,
        target_10s_value_path):
    """Create an economic value raster for all crops given a landcover raster.

    Parameters:
        crop_nutrient_df_path (str): path to CSV with at least the
            column `filenm` and `Pollination dependence crop`.
        yield_and_harea_raster_dir (str): path to a directory that has files
            of the format `[crop_name]_yield.tif` and
            `[crop_name]_harea.tif` where `crop_name` is a value
            in the `filenm` column of `crop_nutrient_df`.
        price_raster_dir (str): path to a directory containing files of the
            form form '[crop name]_price.tif'
        consider_pollination (bool): if True, multiply yields by pollinator
            dependence ratio.
        sample_target_raster_path (path): path to a file that has the raster
            pixel size and dimensions of the desired
            `target_10s_value_path`.
        sample_target_fetch_task (Task): must be complete before
            `sample_target_raster_path` is available.
        target_10km_value_yield_path (str): path to target raster that will
            contain total yield (tons/Ha)
        target_10s_value_yield_path (str): path to a resampled
            `target_10km_value_yield_path` at 10s resolution.
        target_10s_value_path (str): path to target raster that will
            contain a per-pixel amount of pollinator produced `nutrient_name`
            calculated as the sum(
                crop_yield_map * (100-Percent refuse crop) *
                (Pollination dependence crop) * nutrient) * (ha / pixel map))

    Returns:
        None.

    """
    for path in [
            target_10km_value_yield_path, target_10s_value_yield_path,
            target_10s_value_path]:
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass
    crop_nutrient_df = pandas.read_csv(crop_nutrient_df_path)
    yield_raster_path_list = []
    harea_raster_path_list = []
    price_raster_path_list = []
    pollination_yield_factor_list = []
    for _, row in crop_nutrient_df.iterrows():
        yield_raster_path = os.path.join(
            yield_and_harea_raster_dir, "%s_yield.tif" % row['filenm'])
        harea_raster_path = os.path.join(
            yield_and_harea_raster_dir, "%s_harea.tif" % row['filenm'])
        price_raster_path = os.path.join(
            price_raster_dir, '%s_price.tif' % row['filenm'])
        if os.path.exists(yield_raster_path):
            yield_raster_path_list.append(yield_raster_path)
            harea_raster_path_list.append(harea_raster_path)
            price_raster_path_list.append(price_raster_path)
            pollination_yield_factor_list.append(
                (1. - row['Percent refuse'] / 100.))
            if consider_pollination:
                pollination_yield_factor_list[-1] *= (
                    row['Pollination dependence'])
        else:
            raise ValueError(f"not found {yield_raster_path}")

    sample_target_raster_info = pygeoprocessing.get_raster_info(
        sample_target_raster_path)

    yield_raster_info = pygeoprocessing.get_raster_info(
        yield_raster_path_list[0])
    yield_nodata = yield_raster_info['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(yield_nodata, 'raw'), (pollination_yield_factor_list, 'raw')] +
        [(x, 1) for x in yield_raster_path_list + harea_raster_path_list +
         price_raster_path_list], total_price_yield_op, target_10km_value_yield_path,
        gdal.GDT_Float32, yield_nodata)

    y_lat_array = numpy.linspace(
        sample_target_raster_info['geotransform'][3],
        sample_target_raster_info['geotransform'][3] +
        sample_target_raster_info['geotransform'][5] *
        sample_target_raster_info['raster_size'][1],
        sample_target_raster_info['raster_size'][1])

    y_ha_array = area_of_pixel(
        abs(sample_target_raster_info['geotransform'][1]),
        y_lat_array) / 10000.0
    y_ha_column = y_ha_array.reshape((y_ha_array.size, 1))

    pygeoprocessing.warp_raster(
        target_10km_value_yield_path,
        sample_target_raster_info['pixel_size'], target_10s_value_yield_path,
        'cubicspline', target_bb=sample_target_raster_info['bounding_box'],
        n_threads=2)

    # multiplying the ha_array by 1e4 because the of yield are in
    # nutrient / 100g and yield is in Mg / ha.
    pygeoprocessing.raster_calculator(
        [(target_10s_value_yield_path, 1), y_ha_column * 1e4,
         (yield_nodata, 'raw')], density_to_value_op,
        target_10s_value_path, gdal.GDT_Float32, yield_nodata)


def create_prod_nutrient_raster(
        crop_nutrient_df_path, nutrient_name, yield_and_harea_raster_dir,
        consider_pollination, sample_target_raster_path,
        target_10km_yield_path, target_10s_yield_path,
        target_10s_production_path):
    """Create total production & yield for a nutrient for all crops.

    Parameters:
        crop_nutrient_df_path (str): path to CSV with at least the
            column `filenm`, `nutrient_name`, `Percent refuse crop`, and
            `Pollination dependence crop`.
        nutrient_name (str): nutrient name to use to index into the crop
            data frame.
        yield_and_harea_raster_dir (str): path to a directory that has files
            of the format `[crop_name]_yield.tif` and
            `[crop_name]_harea.tif` where `crop_name` is a value
            in the `filenm` column of `crop_nutrient_df`.
        consider_pollination (bool): if True, multiply yields by pollinator
            dependence ratio.
        sample_target_raster_path (path): path to a file that has the raster
            pixel size and dimensions of the desired
            `target_10s_production_path`.
        sample_target_fetch_task (Task): must be complete before
            `sample_target_raster_path` is available.
        target_10km_yield_path (str): path to target raster that will
            contain total yield (tons/Ha)
        target_10s_yield_path (str): path to a resampled
            `target_10km_yield_path` at 10s resolution.
        target_10s_production_path (str): path to target raster that will
            contain a per-pixel amount of pollinator produced `nutrient_name`
            calculated as the sum(
                crop_yield_map * (100-Percent refuse crop) *
                (Pollination dependence crop) * nutrient) * (ha / pixel map))

    Returns:
        None.

    """
    for path in [
            target_10km_yield_path, target_10s_yield_path,
            target_10s_production_path]:
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass
    crop_nutrient_df = pandas.read_csv(crop_nutrient_df_path)
    yield_raster_path_list = []
    harea_raster_path_list = []
    pollination_yield_factor_list = []
    for _, row in crop_nutrient_df.iterrows():
        yield_raster_path = os.path.join(
            yield_and_harea_raster_dir, f"{row['filenm']}_yield.tif")
        harea_raster_path = os.path.join(
            yield_and_harea_raster_dir, f"{row['filenm']}_harea.tif")
        if os.path.exists(yield_raster_path):
            yield_raster_path_list.append(yield_raster_path)
            harea_raster_path_list.append(harea_raster_path)
            pollination_yield_factor_list.append(
                (1. - row['Percent refuse'] / 100.) * row[nutrient_name])
            if consider_pollination:
                pollination_yield_factor_list[-1] *= (
                    row['Pollination dependence'])
        else:
            raise ValueError(f"not found {yield_raster_path}")

    sample_target_raster_info = pygeoprocessing.get_raster_info(
        sample_target_raster_path)

    yield_raster_info = pygeoprocessing.get_raster_info(
        yield_raster_path_list[0])
    yield_nodata = yield_raster_info['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(yield_nodata, 'raw'), (pollination_yield_factor_list, 'raw')] +
        [(x, 1) for x in yield_raster_path_list + harea_raster_path_list],
        total_yield_op, target_10km_yield_path, gdal.GDT_Float32,
        yield_nodata)

    y_lat_array = numpy.linspace(
        sample_target_raster_info['geotransform'][3],
        sample_target_raster_info['geotransform'][3] +
        sample_target_raster_info['geotransform'][5] *
        sample_target_raster_info['raster_size'][1],
        sample_target_raster_info['raster_size'][1])

    y_ha_array = area_of_pixel(
        abs(sample_target_raster_info['geotransform'][1]),
        y_lat_array) / 10000.0
    y_ha_column = y_ha_array.reshape((y_ha_array.size, 1))

    pygeoprocessing.warp_raster(
        target_10km_yield_path,
        sample_target_raster_info['pixel_size'], target_10s_yield_path,
        'cubicspline', target_bb=sample_target_raster_info['bounding_box'],
        n_threads=2)

    # multiplying the ha_array by 1e4 because the of yield are in
    # nutrient / 100g and yield is in Mg / ha.
    pygeoprocessing.raster_calculator(
        [(target_10s_yield_path, 1), y_ha_column * 1e4,
         (yield_nodata, 'raw')], density_to_value_op,
        target_10s_production_path, gdal.GDT_Float32, yield_nodata)


def area_of_pixel(pixel_size, center_lat):
    """Calculate m^2 area of a wgs84 square pixel.

    Adapted from: https://gis.stackexchange.com/a/127327/2397

    Parameters:
        pixel_size (float): length of side of pixel in degrees.
        center_lat (float): latitude of the center of the pixel. Note this
            value +/- half the `pixel-size` must not exceed 90/-90 degrees
            latitude or an invalid area will be calculated.

    Returns:
        Area of square pixel of side length `pixel_size` centered at
        `center_lat` in m^2.

    """
    a = 6378137  # meters
    b = 6356752.3142  # meters
    e = numpy.sqrt(1-(b/a)**2)
    area_list = []
    for f in [center_lat+pixel_size/2, center_lat-pixel_size/2]:
        zm = 1 - e*numpy.sin(numpy.radians(f))
        zp = 1 + e*numpy.sin(numpy.radians(f))
        area_list.append(
            numpy.pi * b**2 * (
                numpy.log(zp/zm) / (2*e) +
                numpy.sin(numpy.radians(f)) / (zp*zm)))
    return pixel_size / 360. * (area_list[0]-area_list[1])


def _mult_raster_op(array_a, array_b, nodata_a, nodata_b, target_nodata):
    """Multiply a by b and skip nodata."""
    result = numpy.empty(array_a.shape, dtype=numpy.float32)
    result[:] = target_nodata
    valid_mask = (array_a != nodata_a) & (array_b != nodata_b)
    result[valid_mask] = array_a[valid_mask] * array_b[valid_mask]
    return result


def mult_rasters(raster_a_path, raster_b_path, target_path):
    """Multiply a by b and skip nodata."""
    raster_info_a = pygeoprocessing.get_raster_info(raster_a_path)
    raster_info_b = pygeoprocessing.get_raster_info(raster_b_path)

    nodata_a = raster_info_a['nodata'][0]
    nodata_b = raster_info_b['nodata'][0]

    if raster_info_a['raster_size'] != raster_info_b['raster_size']:
        aligned_raster_a_path = (
            target_path + os.path.basename(raster_a_path) + '_aligned.tif')
        aligned_raster_b_path = (
            target_path + os.path.basename(raster_b_path) + '_aligned.tif')
        pygeoprocessing.align_and_resize_raster_stack(
            [raster_a_path, raster_b_path],
            [aligned_raster_a_path, aligned_raster_b_path],
            ['near'] * 2, raster_info_a['pixel_size'], 'intersection')
        raster_a_path = aligned_raster_a_path
        raster_b_path = aligned_raster_b_path

    pygeoprocessing.raster_calculator(
        [(raster_a_path, 1), (raster_b_path, 1), (nodata_a, 'raw'),
         (nodata_b, 'raw'), (_MULT_NODATA, 'raw')], _mult_raster_op,
        target_path, gdal.GDT_Float32, _MULT_NODATA)


def add_op(target_nodata, *array_list):
    """Add & return arrays in ``array_list`` but ignore ``target_nodata``."""
    result = numpy.zeros(array_list[0].shape, dtype=numpy.float32)
    valid_mask = numpy.zeros(result.shape, dtype=numpy.bool)
    for array in array_list:
        # nodata values will be < 0
        local_valid_mask = array >= 0
        valid_mask |= local_valid_mask
        result[local_valid_mask] += array[local_valid_mask]
    result[~valid_mask] = target_nodata
    return result


def schedule_sum_and_aggregate(
        task_graph, base_raster_path_list, aggregate_func,
        base_raster_task_list, target_10s_path):
    """Sum all rasters in `base_raster_path` and aggregate to degree."""
    path_list = [(path, 1) for path in base_raster_path_list]

    # this is the only way to get around that we need to check raster size
    for task in base_raster_task_list:
        task.join()

    raster_size_set = set(
        [pygeoprocessing.get_raster_info(path)['raster_size']
         for path in base_raster_path_list])

    dependent_task_list = list(base_raster_task_list)
    if len(raster_size_set) > 1:
        aligned_path_list = [
            target_10s_path + os.path.basename(path) + '_aligned.tif'
            for path in base_raster_path_list]
        align_task = task_graph.add_task(
            func=pygeoprocessing.align_and_resize_raster_stack,
            args=(
                base_raster_path_list,
                aligned_path_list,
                ['near'] * len(base_raster_path_list),
                pygeoprocessing.get_raster_info(
                    base_raster_task_list[0])['pixel_size'], 'intersection'),
            target_path_list=aligned_path_list,
            task_name='align base rasters for %s' % os.path.basename(
                target_10s_path))
        dependent_task_list.append(align_task)
        path_list = [(path, 1) for path in aligned_path_list]
    target_nodata = -9999.

    add_raster_task = task_graph.add_task(
        func=pygeoprocessing.raster_calculator,
        args=(
            [(target_nodata, 'raw')] + path_list, add_op, target_10s_path,
            gdal.GDT_Float32, target_nodata),
        target_path_list=[target_10s_path],
        dependent_task_list=dependent_task_list,
        task_name=f'add rasters {os.path.basename(target_10s_path)}')

    return schedule_aggregate_to_degree(
        task_graph, target_10s_path, aggregate_func, add_raster_task)


def schedule_aggregate_to_degree(
        task_graph, base_raster_path, aggregate_func, base_raster_task,
        target_raster_path_name=None):
    """Schedule an aggregate of 1D approximation."""
    if not target_raster_path_name:
        one_degree_raster_path = (
            base_raster_path.replace('_10s_', '_1d_').replace(
                '_30s_', '_1d_'))
    else:
        one_degree_raster_path = target_raster_path_name
    one_degree_task = task_graph.add_task(
        func=aggregate_to_degree,
        args=(base_raster_path, aggregate_func, one_degree_raster_path),
        target_path_list=[one_degree_raster_path],
        dependent_task_list=[base_raster_task],
        task_name=f'to degree {os.path.basename(one_degree_raster_path)}')
    return (one_degree_task, one_degree_raster_path)


def aggregate_to_degree(raster_path, aggregate_func, target_path):
    """Aggregate input raster to a degree.

    Parameters:
        base_raster_path (string): path to a WGS84 projected raster.
        target_path (string): path to desired target raster that will be
            an aggregated version of `base_raster_path` by the function
            `aggregate_func`.

    Returns:
        None.

    """
    base_raster = gdal.OpenEx(raster_path, gdal.OF_RASTER)
    base_gt = base_raster.GetGeoTransform()
    base_band = base_raster.GetRasterBand(1)
    base_nodata = base_band.GetNoDataValue()

    wgs84sr = osr.SpatialReference()
    wgs84sr.ImportFromEPSG(4326)

    driver = gdal.GetDriverByName('GTiff')
    n_rows = int(
        abs((base_gt[5] * base_band.YSize) / 1.0))
    n_cols = int(
        abs((base_gt[1] * base_band.XSize) / -1.0))
    target_raster = driver.Create(
        target_path, n_cols, n_rows, 1, gdal.GDT_Float32)
    target_raster.SetProjection(wgs84sr.ExportToWkt())
    degree_geotransform = [base_gt[0], 1., 0., base_gt[3], 0., -1.]
    target_raster.SetGeoTransform(degree_geotransform)
    target_band = target_raster.GetRasterBand(1)
    target_band.SetNoDataValue(base_nodata)
    target_band.Fill(base_nodata)

    base_y_winsize = int(round(abs(1. / base_gt[5])))
    base_x_winsize = int(round(abs(1. / base_gt[1])))

    last_time = time.time()
    for row_index in range(n_rows):
        lat_coord = (
            degree_geotransform[3] + degree_geotransform[5] * row_index)
        base_y_coord = int((lat_coord - base_gt[3]) / base_gt[5])
        target_y_coord = int(
            (lat_coord - degree_geotransform[3]) / degree_geotransform[5])
        for col_index in range(n_cols):
            long_coord = (
                degree_geotransform[0] + degree_geotransform[1] * col_index)
            base_x_coord = int((long_coord - base_gt[0]) / base_gt[1])
            target_x_coord = int(
                (long_coord - degree_geotransform[0]) / degree_geotransform[1])

            base_array = base_band.ReadAsArray(
                xoff=base_x_coord, yoff=base_y_coord,
                win_xsize=base_x_winsize, win_ysize=base_y_winsize)
            valid_array = ~numpy.isclose(base_array, base_nodata)
            if valid_array.any():
                target_band.WriteArray(
                    numpy.array([[aggregate_func(base_array[valid_array])]]),
                    xoff=target_x_coord, yoff=target_y_coord)

            current_time = time.time()
            if (current_time - last_time) > 5.0:
                LOGGER.info(
                    "%.2f%% complete", 100.0 * float(row_index+1) / n_rows)
                last_time = current_time
    target_band.FlushCache()
    target_band = None
    target_raster = None
    base_band = None
    base_raster = None
    LOGGER.info("100%% complete")


def sum_num_sum_denom(
        num1_array, num2_array, denom1_array, denom2_array, nodata):
    """Calculate sum of num divided by sum of denom."""
    result = numpy.empty_like(num1_array)
    result[:] = nodata
    valid_mask = (
        ~numpy.isclose(num1_array, nodata) &
        ~numpy.isclose(num2_array, nodata) &
        ~numpy.isclose(denom1_array, nodata) &
        ~numpy.isclose(denom2_array, nodata))
    result[valid_mask] = (
        num1_array[valid_mask] + num2_array[valid_mask]) / (
        denom1_array[valid_mask] + denom2_array[valid_mask] + 1e-9)
    return result


def avg_3_op(array_1, array_2, array_3, nodata):
    """Average 3 arrays. Skip nodata."""
    result = numpy.empty_like(array_1)
    result[:] = nodata
    valid_mask = (
        ~numpy.isclose(array_1, nodata) &
        ~numpy.isclose(array_2, nodata) &
        ~numpy.isclose(array_3, nodata))
    result[valid_mask] = (
        array_1[valid_mask] +
        array_2[valid_mask] +
        array_3[valid_mask]) / 3.
    return result


def weighted_avg_3_op(
        array_1, array_2, array_3,
        scalar_1, scalar_2, scalar_3,
        nodata):
    """Weighted average 3 arrays. Skip nodata."""
    result = numpy.empty_like(array_1)
    result[:] = nodata
    valid_mask = (
        ~numpy.isclose(array_1, nodata) &
        ~numpy.isclose(array_2, nodata) &
        ~numpy.isclose(array_3, nodata))
    result[valid_mask] = (
        array_1[valid_mask]/scalar_1 +
        array_2[valid_mask]/scalar_2 +
        array_3[valid_mask]/scalar_3) / 3.
    return result


def count_ge_one(array):
    """Return count of elements >= 1."""
    return numpy.count_nonzero(array >= 1)


def prop_diff_op(array_a, array_b, nodata):
    """Calculate prop change from a to b."""
    result = numpy.empty_like(array_a)
    result[:] = nodata
    valid_mask = (
        ~numpy.isclose(array_a, nodata) &
        ~numpy.isclose(array_b, nodata))
    # the 1e-12 is to prevent a divide by 0 error
    result[valid_mask] = (
        array_b[valid_mask] - array_a[valid_mask]) / (
            array_a[valid_mask] + 1e-12)
    return result


def calc_relevant_pop(
        base_pop_path,
        prod_poll_indep_en_path, nut_req_en_path,
        prod_poll_indep_fo_path, nut_req_fo_path,
        prod_poll_indep_va_path, nut_req_va_path,
        logic_ufunc, target_raster_path):
    """This function is meant to duplicate the SQL query that Becky wrote:

    ALTER TABLE relevant_population ADD relevant_pop_cur INTEGER;
    ALTER TABLE relevant_population ADD relevant_pop_ssp1 INTEGER;
    ALTER TABLE relevant_population ADD relevant_pop_ssp3 INTEGER;
    ALTER TABLE relevant_population ADD relevant_pop_ssp5 INTEGER;
    UPDATE relevant_population SET relevant_pop_cur = cur;
    UPDATE relevant_population SET relevant_pop_ssp1 = gpw_v4_e_atot_pop_30s_ssp1;
    UPDATE relevant_population SET relevant_pop_ssp3 = gpw_v4_e_atot_pop_30s_ssp3;
    UPDATE relevant_population SET relevant_pop_ssp5 = gpw_v4_e_atot_pop_30s_ssp5;
    UPDATE relevant_population SET relevant_pop_cur = 0 WHERE prod_poll_indep_en_cur > nut_req_en_1d_cur AND prod_poll_indep_fo_cur > nut_req_fo_1d_cur AND prod_poll_indep_va_cur > nut_req_va_1d_cur;
    UPDATE relevant_population SET relevant_pop_ssp1 = 0 WHERE prod_poll_indep_en_ssp1 > nut_req_en_1d_ssp1 AND prod_poll_indep_fo_ssp1 > nut_req_fo_1d_ssp1 AND prod_poll_indep_va_ssp1 > nut_req_va_1d_ssp1;
    UPDATE relevant_population SET relevant_pop_ssp3 = 0 WHERE prod_poll_indep_en_ssp3 > nut_req_en_1d_ssp3 AND prod_poll_indep_fo_ssp3 > nut_req_fo_1d_ssp3 AND prod_poll_indep_va_ssp3 > nut_req_va_1d_ssp3;
    UPDATE relevant_population SET relevant_pop_ssp5 = 0 WHERE prod_poll_indep_en_ssp5 > nut_req_en_1d_ssp5 AND prod_poll_indep_fo_ssp5 > nut_req_fo_1d_ssp5 AND prod_poll_indep_va_ssp5 > nut_req_va_1d_ssp5;

    """

    def relevant_pop_op(
            base_pop_array,
            prod_poll_indep_en_array, nut_req_en_array,
            prod_poll_indep_fo_array, nut_req_fo_array,
            prod_poll_indep_va_array, nut_req_va_array):
        result = numpy.copy(base_pop_array)
        result[
            logic_ufunc(
                (prod_poll_indep_en_array > nut_req_en_array),
                logic_ufunc(
                    (prod_poll_indep_fo_array > nut_req_fo_array),
                    (prod_poll_indep_va_array > nut_req_va_array)))] = 0.0
        return result

    temp_working_dir = tempfile.mkdtemp(
        dir=os.path.dirname(target_raster_path))

    base_raster_path_list = [
        base_pop_path,
        prod_poll_indep_en_path, nut_req_en_path,
        prod_poll_indep_fo_path, nut_req_fo_path,
        prod_poll_indep_va_path, nut_req_va_path]

    aligned_raster_path_list = [
        os.path.join(
            temp_working_dir, '%s_aligned%s' %
            os.path.splitext(os.path.basename(path)))
        for path in base_raster_path_list]

    base_info = pygeoprocessing.get_raster_info(base_pop_path)
    pygeoprocessing.align_and_resize_raster_stack(
        base_raster_path_list, aligned_raster_path_list,
        ['near']*len(base_raster_path_list), base_info['pixel_size'],
        'intersection')

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in aligned_raster_path_list],
        relevant_pop_op, target_raster_path,
        base_info['datatype'], base_info['nodata'][0])

    shutil.rmtree(temp_working_dir, ignore_errors=True)


def build_lookup_from_csv(
        table_path, key_field, to_lower=True, warn_if_missing=True):
    """Read a CSV table into a dictionary indexed by `key_field`.

    Creates a dictionary from a CSV whose keys are unique entries in the CSV
    table under the column named by `key_field` and values are dictionaries
    indexed by the other columns in `table_path` including `key_field` whose
    values are the values on that row of the CSV table.

    Parameters:
        table_path (string): path to a CSV file containing at
            least the header key_field
        key_field: (string): a column in the CSV file at `table_path` that
            can uniquely identify each row in the table.
        to_lower (bool): if True, converts all unicode in the CSV,
            including headers and values to lowercase, otherwise uses raw
            string values.
        warn_if_missing (bool): If True, warnings are logged if there are
            empty headers or value rows.

    Returns:
        lookup_dict (dict): a dictionary of the form {
                key_field_0: {csv_header_0: value0, csv_header_1: value1...},
                key_field_1: {csv_header_0: valuea, csv_header_1: valueb...}
            }

        if `to_lower` all strings including key_fields and values are
        converted to lowercase unicode.

    """
    # Check if the file encoding is UTF-8 BOM first, related to issue
    # https://bitbucket.org/natcap/invest/issues/3832/invest-table-parsing-does-not-support-utf
    encoding = None
    with open(table_path) as file_obj:
        first_line = file_obj.readline()
        if first_line.startswith('\xef\xbb\xbf'):
            encoding = 'utf-8-sig'
    table = pandas.read_csv(
        table_path, sep=None, engine='python', encoding=encoding)
    header_row = list(table)
    try:  # no unicode() in python 3
        key_field = unicode(key_field)
    except NameError:
        pass
    if to_lower:
        key_field = key_field.lower()
        header_row = [
            x if not isinstance(x, str) else x.lower()
            for x in header_row]

    if key_field not in header_row:
        raise ValueError(
            '%s expected in %s for the CSV file at %s' % (
                key_field, header_row, table_path))
    if warn_if_missing and '' in header_row:
        LOGGER.warn(
            "There are empty strings in the header row at %s", table_path)

    key_index = header_row.index(key_field)
    lookup_dict = {}
    for index, row in table.iterrows():
        if to_lower:
            row = pandas.Series([
                x if not isinstance(x, str) else x.lower()
                for x in row])
        # check if every single element in the row is null
        if row.isnull().values.all():
            LOGGER.warn(
                "Encountered an entirely blank row on line %d", index+2)
            continue
        if row.isnull().values.any():
            row = row.fillna('')
        lookup_dict[row[key_index]] = dict(zip(header_row, row))
    return lookup_dict


def download_and_unzip(url, target_dir, target_token_path):
    """Download `url` to `target_dir` and touch `target_token_path`."""
    zipfile_path = os.path.join(target_dir, os.path.basename(url))
    LOGGER.debug('url %s, zipfile_path: %s', url, zipfile_path)
    ecoshard.download_url(url, zipfile_path)

    with zipfile.ZipFile(zipfile_path, 'r') as zip_ref:
        zip_ref.extractall(target_dir)

    with open(target_token_path, 'w') as touchfile:
        touchfile.write(f'unzipped {zipfile_path}')


def create_price_raster(
        base_raster_path, country_vector_path, country_crop_price_map,
        crop_name, target_crop_price_raster_path):
    """Rasterize countries as prices.

    Parameters:
        base_raster_path (str): path to a raster that will be the base shape
            for `target_crop_price_ratser_path`.
        country_vector_path (str): path to country shapefile that has a
            field called `ISO3` that corresponds to the first key in
            `price_map`.
        country_crop_price_map (dict): map that indexes country ISO names to
            `crop_name`. If this crop is not in the subdictionary there is
            no price for that crop in that country.
        target_crop_price_raster_path (str): a raster with pixel values
            corresponding to the country in which the pixel resides and
            the price of that crop in the country.

    Returns:
        None.

    """
    pygeoprocessing.new_raster_from_base(
        base_raster_path, target_crop_price_raster_path, gdal.GDT_Float32,
        [-1], fill_value_list=[-1])
    memory_driver = ogr.GetDriverByName('MEMORY')
    country_vector = gdal.OpenEx(country_vector_path, gdal.OF_VECTOR)
    country_layer = country_vector.GetLayer()

    price_vector = memory_driver.CreateDataSource('price_vector')
    spat_ref = country_layer.GetSpatialRef()

    price_layer = price_vector.CreateLayer(
        'price_layer', spat_ref, ogr.wkbPolygon)
    price_layer.CreateField(ogr.FieldDefn('price', ogr.OFTReal))
    price_layer_defn = price_layer.GetLayerDefn()
    # add polygons to subset_layer
    price_layer.StartTransaction()
    for country_feature in country_layer:
        country_name = country_feature.GetField('ISO3')
        if crop_name in country_crop_price_map[country_name]:
            country_geom = country_feature.GetGeometryRef()
            new_feature = ogr.Feature(price_layer_defn)
            new_feature.SetGeometry(country_geom.Clone())
            new_feature.SetField(
                'price', country_crop_price_map[country_name][crop_name])
            price_layer.CreateFeature(new_feature)
    price_layer.CommitTransaction()

    target_crop_raster = gdal.OpenEx(
        target_crop_price_raster_path, gdal.OF_RASTER | gdal.GA_Update)
    gdal.RasterizeLayer(
        target_crop_raster, [1], price_layer,
        options=['ATTRIBUTE=price'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='NCI Pollination Analysis')
    parser.add_argument(
        'landcover rasters', type=str, nargs='+',
        help=(
            'Paths or patterns to landcover rasters that use ESA style '
            'encoding.'))
    args = parser.parse_args()
    landcover_raster_list = []
    for glob_pattern in vars(args)['landcover rasters']:
        for raster_path in glob.glob(glob_pattern):
            print(raster_path)
            # just try to open it in case it's not a raster, we won't see
            # anything otherwise
            r = gdal.OpenEx(raster_path, gdal.OF_RASTER)
            r = None
            landcover_raster_list.append(raster_path)

    for landcover_path in landcover_raster_list:
        calculate_for_landcover(landcover_path)
