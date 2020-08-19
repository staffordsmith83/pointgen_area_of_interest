#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VERSION 3.0
Created on Tuesday 18th of August 2020
Takes a project boundary and up to three exclusion layers.
Creates an Area of Interest by buffering all inputs.
Generates random points within the Area of Interest.
These random points align with the centre of a hypothetical
Sentinel raster pixel available at that location.

@author: Stafford Smith
"""

# %%###########################################################################
# IMPORTS
import datetime as dt
import fiona
import geopandas as gpd
import logging
import math
import matplotlib.pyplot as plt
import os
import random
import rasterio
import sys
from fiona import transform
from fiona.crs import from_epsg
from pyproj import CRS
from shapely.geometry import shape, Point, mapping, MultiPolygon
from shapely.ops import unary_union

# %%###########################################################################
# SET ENVIRONMENT VARIABLES - OPTIONAL
""" This is used if your python installation cannot find the proj.db library etc. 
Common with Anaconda installations on PC"""

# os.environ['GDAL_DATA'] = os.environ['CONDA_PREFIX'] + r'\Library\share\gdal'
# os.environ['PROJ_LIB'] = os.environ['CONDA_PREFIX'] + r'\Library\share'

# %%###########################################################################
# SET CONSTANTS
# Please provide inputs in GDA94 (lat and long).

PROJECT_SHAPEFILE = ''  # shp or gpkg

# The next three layers are optional. Leave as empty string '' if not required.
EXCLUSION1 = ''  # shp or gpkg

EXCLUSION2 = ''  # shp or gpkg

EXCLUSION3 = ''  # shp or gpkg

LAND_COVER = ''  # soil landscape mapping vector layer. 
FIELDS_TO_ADD = ['geometry', 'mu_name', 'mu_sum_des', 'mu_land_t1']  # list of the fields to join from LAND_COVER

# Please provide EITHER a PROJECT name to look up the classification raster path in the lookup table,
# a filepath to a raster classification, or a filepath to a vector classification (code will take the first specified in this order).

PROJECT_NAME = ''  # please do not specify unless you intend to use the classification file lookup table

RASTER_CLASSIFICATION = ''  # tif

VECTOR_CLASSIFICATION = ''  # shp or gpkg

# Set plots and buffers

NUM_POINTS = 1000  # int

BOUNDARY_BUFFER = -30  # int

EXCLUSION_BUFFER = 30  # int

# SEED variable takes any integer. So you can replace the line below for a manual seed.
SEED = int(dt.date.today().strftime('%Y%m%d'))  # uses today's date, expressed as yyyymmdd. Replace if necessary.

# %%###########################################################################
# PROJECTS Dictionary
# Key is PROJECT name, value is Class filename.
# Please don't use the reserved backslash character.
PROJECTS_DICT = {'Xxxxxxxx': 'C:/pathtofile',
                 'Yyyyyyyy': 'C:/pathtofile'
                 }


# %%###########################################################################
# FUNCTIONS:


def find_utm_zone(lon, lat):
    """ from https://stackoverflow.com/questions/40132542/get-a-cartesian-projection-accurate-around-a-lat-lng-pair """
    utm_band = str((math.floor((lon + 180) / 6) % 60) + 1)
    if len(utm_band) == 1:
        utm_band = '0' + utm_band
    if lat >= 0:
        epsg_code = '326' + utm_band
    else:
        epsg_code = '327' + utm_band
    crs_utm = fiona.crs.from_epsg(epsg_code)

    return crs_utm


def my_buffer(geoms, amount):
    # takes a list of geometries, buffers, returns list of geometries.
    # should be independent of geometry type...
    # checks for empty geometries due to negative buffering...
    geom_list = []

    for geom in geoms:
        geom = shape(geom)
        buffered = geom.buffer(amount)
        if not buffered.is_valid:
            print(f'Invalid geometry resulting from buffer operation using {amount}')
        if not buffered.is_empty:
            # plt.plot(*buffered.exterior.xy)
            geom_list.append(buffered)
    # plt.title(src)
    # plt.show()
    return geom_list


def roundTo20plus10(x):
    base = 20
    rounded = x // base * base
    return rounded + 10


def generate_random(number, polygon):
    points = []
    minx, miny, maxx, maxy = polygon.bounds
    while len(points) < number:
        pnt = Point(roundTo20plus10(random.uniform(minx, maxx)), roundTo20plus10(random.uniform(miny, maxy)))
        if polygon.contains(pnt):
            points.append(pnt)
    return points


def output_to_file(dst_file, geoms, src_crs, dst_crs):
    print('Writing point file...')

    with fiona.open(dst_file,
                    'w',
                    crs=dst_crs,
                    driver='ESRI Shapefile',
                    schema={'geometry': 'Point',
                            'properties': {'POINT_NUM': 'int',
                                           'X': 'float',
                                           'Y': 'float'}}) as dst:
        for counter, point in enumerate(geoms):
            # reproject if necessary
            new_point = shape(transform.transform_geom(src_crs, dst_crs, mapping(point)))

            prop = {'POINT_NUM': counter + 1,
                    'X': new_point.x,
                    'Y': new_point.y}
            dst.write({'geometry': mapping(new_point), 'properties': prop})


def output_poly_to_file(dst_file, polygons, src_crs, dst_crs):
    print('Writing polygon file...')

    with fiona.open(dst_file,
                    'w',
                    crs=dst_crs,
                    driver='ESRI Shapefile',
                    schema={'geometry': 'Polygon',
                            'properties': {'id': 'int'}}) as dst:
        for counter, polygon in enumerate(polygons):
            # reproject if necessary
            new_polygon = transform.transform_geom(src_crs, dst_crs, mapping(polygon))
            new_polygon = shape(new_polygon)

            prop = {'id': counter}
            dst.write({'geometry': mapping(new_polygon), 'properties': prop})


def get_utm_geoms(filepath):
    with fiona.open(filepath) as src:
        meta = src.meta

        bounds = src.bounds

        zone1 = find_utm_zone(bounds[0], bounds[1])
        zone2 = find_utm_zone(bounds[2], bounds[3])
        if not zone1 == zone2:
            print('Data is not contained within a single utm zone. Exiting program.')
            sys.exit()

        # reproject
        src_crs = meta['crs']
        dst_crs = zone1
        new_geoms = []
        for feature in src:
            geom = feature['geometry']
            if geom:  # this line should exclude empty geometries
                new_geom = transform.transform_geom(src_crs, dst_crs, geom)
                new_geoms.append(new_geom)
        return new_geoms


def get_raster_vals(points, raster):
    """Takes a geodataframe, and Returns a geodataframe"""
    print('Extracting raster values...')

    # Read points from shapefile
    pts = points
    pts.index = range(len(pts))
    coords = [(x, y) for x, y in zip(pts.X, pts.Y)]

    # Open the raster and store metadata
    src = rasterio.open(raster)

    # Sample the raster at every point location and store values in DataFrame
    pts['Class'] = [x[0] for x in src.sample(coords)]
    return pts


def spatial_join(points, polys, fields_to_add):
    """Takes geodataframe inputs, Returns a geodataframe
    NB: Will throw an error if the polygon file contains invalid geometries"""
    # TODO: make sure that the points and polys are in the same projection
    if len(fields_to_add) > 0:
        polys_subset = polys[fields_to_add]
    else:
        polys_subset = polys

    df = gpd.sjoin(points, polys_subset, op='within')
    df = df.drop(['index_right'], axis=1)

    return df


def check_projection(filepath, reference):
    """Takes a filepath and a projection name string as input"""
    file_crs = False
    filename = (os.path.basename(filepath))
    # Check if the file is a vector
    try:
        with fiona.open(filepath) as vector_file:
            v_crs = CRS.from_wkt(vector_file.crs_wkt)
            file_crs = v_crs.name
    except:
        # Check if the file is a raster
        try:
            with rasterio.open(filepath) as raster_file:
                r_crs = CRS.from_wkt(raster_file.crs.wkt)
                file_crs = r_crs.name
        except:
            print(f'{filepath} cannot be opened as either a raster or vector... Cannot check projection.')

    # Compare file crs to reference epsg
    if not file_crs == reference:
        print(f'*****WARNING: {filename} may not be in GDA94, '
              f'script will attempt to continue, but better to fix data and rerun!')


# %%###########################################################################
# MAIN:
def main():
    # setup outputs
    # global is required as this out of scope variable is modified within the main function sometimes:
    global RASTER_CLASSIFICATION
    dst_file = f'{os.path.splitext(PROJECT_SHAPEFILE)[0]}_plots_SEED-{SEED}{os.path.splitext(PROJECT_SHAPEFILE)[1]}'
    area_dst_file = f'{os.path.splitext(PROJECT_SHAPEFILE)[0]}_AOI{os.path.splitext(PROJECT_SHAPEFILE)[1]}'
    all_values_file = f'{os.path.splitext(PROJECT_SHAPEFILE)[0]}_random_plots_all{os.path.splitext(PROJECT_SHAPEFILE)[1]}'
    target_only_file = f'{os.path.splitext(PROJECT_SHAPEFILE)[0]}_plots_target{os.path.splitext(PROJECT_SHAPEFILE)[1]}'

    # setup logging
    dst_log_file = f'{os.path.splitext(PROJECT_SHAPEFILE)[0]}_plots_SEED-{SEED}.log'
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logging.basicConfig(filename=dst_log_file,
                        level=logging.INFO,
                        filemode='w',
                        format='%(asctime)s %(message)s')

    logger.info(f'\nPROJECT_SHAPEFILE:\t{PROJECT_SHAPEFILE}\n' + \
                f'SEED:\t\t\t\t{SEED}\n' + \
                f'PROJECT BOUNDARY BUFFER:\t{BOUNDARY_BUFFER}\n' + \
                f'EXCLUSION BUFFER:\t{EXCLUSION_BUFFER}\n')

    # set the random seed
    random.seed(SEED)

    # get some metadata from the PROJECT boundary file
    check_projection(PROJECT_SHAPEFILE, reference='GDA94')
    with fiona.open(PROJECT_SHAPEFILE) as project_boundary:
        meta = project_boundary.meta
        initial_crs = meta['crs']
        bounds = project_boundary.bounds
        utm_crs = find_utm_zone((bounds[2] + bounds[0]) / 2, (bounds[3] + bounds[1]) / 2)

    # Buffer the project boundary
    project_boundary = get_utm_geoms(PROJECT_SHAPEFILE)
    buffered_project = my_buffer(project_boundary, BOUNDARY_BUFFER)
    print(f'Buffered {len(buffered_project)} project polygons. Empty ones excluded.')

    # Buffer the exclusion
    if not EXCLUSION1 == '':
        check_projection(EXCLUSION1, reference='GDA94')
        exclu1 = get_utm_geoms(EXCLUSION1)
        buffered_exclu1 = my_buffer(exclu1, EXCLUSION_BUFFER)
        print(f'Buffered {len(buffered_exclu1)} exclu1.')
        logger.info(f'\nExcluded exclusion file {EXCLUSION1}\n')
    else:
        buffered_exclu1 = []

    if not EXCLUSION2 == '':
        check_projection(EXCLUSION2, reference='GDA94')
        exclu3 = get_utm_geoms(EXCLUSION2)
        buffered_exclu3 = my_buffer(exclu3, EXCLUSION_BUFFER)
        print(f'Buffered {len(buffered_exclu3)} exclu3.')
        logger.info(f'\nExcluded exclusion file {EXCLUSION2}\n')
    else:
        buffered_exclu3 = []

    if not EXCLUSION3 == '':
        check_projection(EXCLUSION3, reference='GDA94')
        exclu2 = get_utm_geoms(EXCLUSION3)
        buffered_exclu2 = my_buffer(exclu2, EXCLUSION_BUFFER)
        print(f'Buffered {len(buffered_exclu2)} exclu2.')
        logger.info(f'\nExcluded exclusion file {EXCLUSION3}\n')
    else:
        buffered_exclu2 = []

    # Combine the buffered exclusion layers
    all_cut_features = buffered_exclu1 + buffered_exclu2 + buffered_exclu3
    print(f'Total number of features to cut out of project is {len(all_cut_features)}.')

    # Clip our areas to make one 'area of interest':
    # NB: Make features Multipolygon objects to process the unarary union easier
    # NB: inputs to the MultiPolygon function should not already be multipolygons...

    # This conditional statement to deal with the case where no exclusion inputs
    if len(all_cut_features) == 0:
        area_of_interest = MultiPolygon(buffered_project)

    else:
        if all_cut_features[0].geom_type == 'Polygon':
            all_cut_features_multi = MultiPolygon(all_cut_features)

        all_cut_features_multi = unary_union(all_cut_features)
        buffered_project_multi = MultiPolygon(buffered_project)
        area_of_interest = unary_union(buffered_project_multi.difference(all_cut_features_multi))

    # output the AOI shapefile for records.
    output_poly_to_file(area_dst_file, area_of_interest, utm_crs, initial_crs)

    # generate the random points,
    # the 'generate_random' function will position them in the center of a hypothetical Sentinel pixel.
    # This is useful as these points will be used in conjunction with a landcover classification using Sentinel data.
    logger.info(f'\nStarted generation of random plots within area of interest...\n')
    points = generate_random(NUM_POINTS, area_of_interest)
    logger.info(f'\nGenerated {len(points)} points.\n')

    # Plot the output for review
    plt.title("Area of Interest and Points")
    for geom in area_of_interest:
        plt.plot(*geom.exterior.xy)
    plt.scatter([point.x for point in points], [point.y for point in points], s=1)
    plt.show()

    # output the points
    output_to_file(dst_file, points, utm_crs, initial_crs)

    # Join attributes from landcover
    # TODO: deal with bad landcover filename
    if not LAND_COVER == '':
        check_projection(LAND_COVER, reference='GDA94')
        points_for_join = gpd.read_file(dst_file)
        polys_for_join = gpd.read_file(LAND_COVER)
        print("Adding attributes from Landcover layer...")
        pts_landcover_added = spatial_join(points_for_join, polys_for_join,
                                           fields_to_add=FIELDS_TO_ADD)  # specified above as a constant

        if pts_landcover_added.empty:
            print('WARNING: landcover values not successfully added')
            pts_landcover_added = points_for_join
        else:
            logger.info(f'\nLandcover values added from {LAND_COVER}\n')
    else:
        pts_landcover_added = gpd.read_file(dst_file)

    ###################################
    # GET CLASSIFICATION VALUES SECTION

    # Check if PROJECT_NAME has been specified
    if not PROJECT_NAME == '':
        print('PROJECT_NAME detected... Looking up Classification file...')

        # give RASTER_CLASSIFICATION a value.
        # If it is already defined, it will be overwritten by the file from the lookup table
        try:
            RASTER_CLASSIFICATION = PROJECTS_DICT[PROJECT_NAME]

        except KeyError:
            print(f'({PROJECT_NAME}) does not exist in lookup table, script will continue and check for other options.')

    # now checks if RASTER_CLASSIFICATION has a value yet, either manually specified or looked up...
    if not RASTER_CLASSIFICATION == '':
        check_projection(RASTER_CLASSIFICATION, reference='GDA94')
        try:
            with open(RASTER_CLASSIFICATION, 'r') as raster_classification:
                print(f'{os.path.basename(RASTER_CLASSIFICATION)} successfully opened')

            pts_classification_added = get_raster_vals(pts_landcover_added, RASTER_CLASSIFICATION)
            # export file of all points, and one with just those in the TARGET area (classification ==1)
            print('Writing final output shapefiles - TARGET status from raster')
            pts_classification_added.to_file(all_values_file, encoding='utf-8')
            pts_target_only = pts_classification_added.loc[pts_classification_added['Class'] == 1]
            pts_target_only.to_file(target_only_file, encoding='utf-8')

            logger.info(f'\nClassification values added from {RASTER_CLASSIFICATION}\n')

        except FileNotFoundError:
            print('Raster Class filepath is invalid, script will continue and check for other options.')

    # VECTOR VERSION
    elif not VECTOR_CLASSIFICATION == '':
        check_projection(VECTOR_CLASSIFICATION, reference='GDA94')
        classification_polys = gpd.read_file(VECTOR_CLASSIFICATION)

        # find classification column
        classification_polys_columns = classification_polys.columns.to_list()
        classification_field = None
        while classification_field not in classification_polys_columns:
            print(f'The columns of the VECTOR_CLASSIFICATION file are {classification_polys_columns}')
            classification_field = input(
                'Please enter the exact name of the column that contains the relevant classification values:')

        # Add that field:
        pts_classification_added = spatial_join(pts_landcover_added, classification_polys,
                                                fields_to_add=['geometry', classification_field])
        logger.info(f'\nClassification values added from {VECTOR_CLASSIFICATION}\n')

        # find points in target and export target only
        classification_field_values = classification_polys[classification_field].unique().tolist()

        # remember data type of the TARGET field(to deal with user specifying integers as strings)
        classification_column_dtype = type(classification_field_values[0])

        # convert our list of the classification field values to strings, so that the user can choose between them
        classification_field_values = [str(i) for i in classification_field_values]

        # make a blank list, for the user to specify values that mean a polygon is TARGET
        # this will be built as a list of strings, and later converted back to a numeric type if they were initially.
        target_designator_list = []
        target_designator = None
        add_another_value = 'True'

        while add_another_value:
            print('The unique values of the classification column are:')
            for s in classification_field_values:
                print(s)

            # keeps asking until user provides a valid input, will also accept a blank string (return)
            while target_designator not in classification_field_values + ['']:
                target_designator = input(
                    'Please enter the exact name of the value that indicates a polygon is TARGET\n'
                    '(hit return to skip):')
            # allow a user to hit return in the last step and skip addition of values. The lists are not modified.
            if not target_designator == '':
                target_designator_list.append(target_designator)  # Add the value to the list
                classification_field_values.remove(
                    target_designator)  # Remove the value from the list of options to add

            # while loop limits the user responses to 'Y' or ''
            while not add_another_value in ['', 'Y']:
                print(f'Your specified TARGET values are {target_designator_list}')
                add_another_value = input('Please type Y if you would like to add another value to specify TARGET\n'
                                          '(hit return to skip):')
            if add_another_value == '':
                add_another_value = False
            elif add_another_value == 'Y':
                add_another_value = True

        # export file of all points, and one with just those in the TARGET area (classification == target_designator)
        print('Writing final output shapefiles - TARGET status from vector')
        pts_classification_added.to_file(all_values_file, encoding='utf-8')

        # convert the TARGET designator list back to its original datatype
        # (stored in classification_field_dtype variable)
        # this happens if the classification field datatype was numeric in the original shapefile
        target_designator_list = [classification_column_dtype(i) for i in target_designator_list]

        pts_target_only = pts_classification_added.loc[
            pts_classification_added[classification_field].isin(target_designator_list)]
        pts_target_only.to_file(target_only_file, encoding='utf-8')

    # NO CLASSIFICATION PROVIDED
    else:
        print("No Classification provided as raster or vector layer. TARGET subset not exported.")
        pts_landcover_added.to_file(all_values_file, encoding='utf-8')


# Run the main function
main()
