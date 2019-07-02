#!/usr/bin/python3

import pandas as pd
import geopandas as gpd
import numpy as np
import os
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon, Point, box
import datetime
import pickle
import json
import copy

from pyhdf.SD import SD, SDC
from pyhdf.HDF import *
from pyhdf.VS import *

from netCDF4 import Dataset
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.transform import from_bounds

from global_variables import data_path, default_crs

from global_functions_and_classes import specify_area_of_interest_EPSG3995
from global_functions_and_classes import radarPass, save_radarPass_object_pkl


def read_datasets(FILENAME):
    """Reads all datasets in an HDF-EOS file, replaces fill_values with NaNs, 
    and places those datasets in a dictionary.
    
    PARAMETERS: 
    -----------
    FILENAME: Name of HDF-EOS file. 

    OUTPUTS: 
    --------
    data_dict: Dictionary of datasets.  
    """

    # Open the file
    file = SD(FILENAME, SDC.READ)
    datasets_dict = file.datasets()
    dataset_names = [dname for dname in datasets_dict.keys()]

    # Extract a single dataset from the file, replace the fill values with NaN.
    data_dict = {}
    for ds in dataset_names:
        sds_obj = file.select(ds)
        data = sds_obj.get()
        fill_val = sds_obj.attributes().get('_FillValue')
        if fill_val:
            data = data.astype(np.float32)
            data[data == fill_val] = np.NaN
        data_dict[ds] = data

    return data_dict


def read_vdata(FILENAME):
    """Reads all vdata fields in an HDF-EOS file and places those fields 
    in a dictionary.
    
    PARAMETERS: 
    -----------
    FILENAME: Name of HDF-EOS file. 

    OUTPUTS: 
    --------
    vdata_dict: Dictionary of vdata fields.  
    """

    # Prepare to read the data.
    f = HDF(FILENAME, SDC.READ)        # Open the file
    vs = f.vstart()                    # Start the vdata interface
    data_info_list = vs.vdatainfo()    # List the vdata fields
    vdata_fieldnames = [a[0] for a in data_info_list]    # Get the names

    # Load the data, place in dictionary
    vdata_dict = {}
    for field in vdata_fieldnames:
        vdata_dict[field] = np.squeeze(np.asarray(vs.attach(field)[:]))

    # terminate the vdata interface, close the file.
    vs.end()
    f.close()

    return vdata_dict


def add_datetime(vdata_dict, FILENAME):
    """Adds a datetime vector (for the time of each measurement) to a 
    vdata_dict object. 
    
    PARAMETERS: 
    -----------
    vdata_dict: vdata_dict object. 

    FILENAME: Name of file that the vdata_dict object came from. 

    OUTPUTS: 
    --------
    vdata_dict: vdata_dict object with datetime vector. 
    """

    first_second = np.around(vdata_dict['UTC_start'], decimals=-2)
    first_dtime = np.asarray(datetime.datetime.strptime(FILENAME.split('/')[-1][:13], '%Y%j%H%M%S'))\
        .astype('datetime64[D]')+pd.Timedelta(str(first_second)+' seconds')
    tv = first_dtime + \
        (np.around(vdata_dict['Profile_time'],
                   decimals=2)*pd.Timedelta('1 seconds'))
    time_vec = np.asarray([pd.Timestamp(t) for t in tv])
    vdata_dict['datetime'] = time_vec
    return vdata_dict


def read_hdfeos_file(FILENAME):
    """Loads an entire HDF-EOS file and places the contents in a 
    more accessible format. 
    
    PARAMETERS: 
    -----------
    FILENAME: Name of the HDF-EOS file. 

    OUTPUTS: 
    --------
    data_dict: dictionary containing all the datasets in 
    the HDF-EOS file. 

    vdata_dict: dictionary containing all the vdata fields
    in the HDF-EOS file. 
    """

    data_dict = read_datasets(FILENAME)
    vdata_dict = read_vdata(FILENAME)
    vdata_dict = add_datetime(vdata_dict, FILENAME)

    return data_dict, vdata_dict


def get_all_downloaded_pass_ids(data_path):
    """Gets ID numbers for every CloudSat file that has been downloaded.
    
    PARAMETERS:
    ----------- 
    data_path: Path to main directory where all data is stored. Specified
    in global_variables.py

    OUTPUTS: 
    --------
    id_list: List of all ID numbers. 
    """
    id_list = sorted(list(set([a[:19] for a in os.listdir(
        os.path.join(data_path, 'cloudsat')) if a[0] != '.'])))
    return id_list


def create_radarPass_instance(data_path, pass_id):
    """Creates a radarPass instance using a GEOPROF CloudSat file.

    PARAMETERS: 
    -----------
    data_path: Path to main directory where all data is stored. Specified
    in global_variables.py

    pass_id: ID number for the given radar pass. Returned by 
    'get_all_downloaded_pass_ids' function. 

    OUTPUTS: 
    --------
    radar_pass: radarPass object. 
    """

    geoprof_fname = os.path.join(data_path, 'cloudsat',
                                 pass_id+'_CS_2B-'+'GEOPROF'+'_GRANULE_P1_R05_E06_F00.hdf')
    data_dict, vdata_dict = read_hdfeos_file(geoprof_fname)
    radar_pass = radarPass(vdata_dict['Longitude'], vdata_dict['Latitude'],
                           vdata_dict['datetime'], np.ravel(
                               np.nanmean(data_dict['Height'].T, 1)[::-1]),
                           np.flipud(data_dict['Radar_Reflectivity'].T),
                           np.flipud(data_dict['CPR_Cloud_mask'].T))

    return radar_pass


def add_cldclass_to_radarPass(radar_pass, data_path, pass_id):
    """Adds cloud classification data to a radarPass instance."""

    cldclass_fname = os.path.join(data_path, 'cloudsat',
                                  pass_id+'_CS_2B-'+'CLDCLASS'+'_GRANULE_P1_R05_E06_F00.hdf')
    data_dict, vdata_dict = read_hdfeos_file(cldclass_fname)
    radar_pass = radar_pass.add_cloudclass(data_dict, vdata_dict)
    
    return radar_pass


def add_cwc_to_radarPass(radar_pass, data_path, pass_id):
    """Adds cloud water content data to a radarPass instance."""

    cwc_fname = os.path.join(data_path, 'cloudsat',
                             pass_id+'_CS_2B-'+'CWC-RO'+'_GRANULE_P1_R05_E06_F00.hdf')
    data_dict, vdata_dict = read_hdfeos_file(cwc_fname)
    radar_pass = radar_pass.add_cwc(data_dict, vdata_dict)

    return radar_pass


if __name__ == '__main__': 
    
    area_of_interest = specify_area_of_interest_EPSG3995()
    pass_ids = get_all_downloaded_pass_ids(data_path)

    for pass_id in pass_ids[:5]:

        print("Working on pass ID: "+pass_id)

        # Create full radarPass object
        radar_pass_f = create_radarPass_instance(data_path, pass_id)
        radar_pass_f = add_cldclass_to_radarPass(radar_pass_f, data_path, pass_id)
        radar_pass_f = add_cwc_to_radarPass(radar_pass_f, data_path, pass_id)
        radar_pass_f = radar_pass_f.trim_pass(area_of_interest)
        radar_pass_f = radar_pass_f.add_era5_data()

        # Create smaller version: reduce resolution, remove less useful fields. 
        radar_pass = copy.deepcopy(radar_pass_f)
        radar_pass = radar_pass.reduce_size_cloudsat(reduction_factor=3)
        radar_pass = radar_pass.reduce_size_era5(reduction_factor=4)
        pass_timestamp = save_radarPass_object_pkl(radar_pass, data_path)



    


