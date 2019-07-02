#!/usr/bin/python3

import os 
import pickle 
import geopandas as gpd
import pandas as pd 
import numpy as np 
from shapely.geometry import Polygon, MultiPolygon, Point, box
import json
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.transform import from_bounds
import datetime
from netCDF4 import Dataset

from global_variables import data_path,default_crs


def save_radarPass_object_pkl(radar_pass, data_path):
    """Saves a radarPass instance to a pickle file with a standard filename format."""
    pass_timestamp = radar_pass.cloudsat['timestamp'][len(radar_pass.cloudsat['timestamp'])//2]\
                               .strftime('%Y%m%d_%H%M%S')
    fname = 'radarPass_plot_'+pass_timestamp+'.pkl'
    with open(os.path.join(data_path, 'radar_passes', fname), 'wb') as f:
        pickle.dump(radar_pass, f, 2)
    return pass_timestamp


def load_radarPass_object_pkl(pass_timestamp, data_path):
    """Loads a radarPass instance from a pickle file with a standard filename format."""
    fname = 'radarPass_plot_'+pass_timestamp+'.pkl'
    with open(os.path.join(data_path, 'radar_passes', fname), 'rb') as f:
        radar_pass = pickle.load(f)
    return radar_pass


def save_radarPass_object_json(radar_pass, data_path):
    """Saves a radarPass instance to a json file with a standard filename format."""
    pass_timestamp = radar_pass.cloudsat['timestamp'][len(radar_pass.cloudsat['timestamp'])//2]\
                               .strftime('%Y%m%d_%H%M%S')
    fname = 'radarPass_plot_'+pass_timestamp+'.json'
    radar_dict = radar_pass.get_json_serializable_obj()
    radar_dict['timestamp'] = [tstamp.strftime(
        '%Y%m%d_%H%M%S.%f') for tstamp in radar_dict['timestamp']]
    with open(os.path.join(data_path, 'radar_passes', fname), 'w') as f:
        json.dump(radar_dict, f)
    return pass_timestamp


def specify_area_of_interest_EPSG4326(bbox=(-180, 60, 180, 90)):
    """Specifies the area of interest for the plot in the lat/lon coordinate system"""
    area_of_int = gpd.GeoDataFrame({'geometry': [box(*bbox)]})
    area_of_int.crs = {'init': 'epsg:4326'}
    return area_of_int


def specify_area_of_interest_EPSG3995(bbox=(-3e6, -3e6, 3e6, 3e6)):
    """Specifies the area of interest for the plot in the Polar
    Stereographic coordinate system. 
    """
    area_of_int = gpd.GeoDataFrame({'geometry': [box(*bbox)]})
    area_of_int.crs = {'init': default_crs}
    return area_of_int


def load_country_geometries(area_of_interest):
    """Loads polygons for all of the world's countries, converts them to the default_crs."""
    countries = gpd.read_file(os.path.join(
        data_path, 'Countries_WGS84', 'Countries_WGS84.shp'))
    countries = countries.to_crs({'init': default_crs})
    countries = gpd.overlay(countries, area_of_interest, how='intersection')
    return countries


def prepare_polygon_coords_for_bokeh(countries):
    """Prepares the country polygons for plotting with Bokeh. 
    
    To plot series of polygons, Bokeh needs two lists of lists (one for x coordinates, and another
    for y coordinates). Each element in the outer list represents a single polygon, and each 
    element in the inner lists represents the coordinate for a single point in given polygon. 
    This function takes a GeoDataFrame with a given set of countries, and returns Bokeh-friendly
    lists of x coordinates and y coordinates for those countries. 

    PARAMETERS: 
    -----------
    countries: GeoDataFrame with a given set of countries. 

    OUTPUTS: 
    --------
    x_coords, y_coords: Bokeh-friendly lists of x and y coordinates for those countries. 
    """

    # Simplify shapes (to resolution of 10000 meters), convert polygons to multipolygons.
    list_of_polygons = []
    for raw_poly in countries['geometry']:
        raw_poly = raw_poly.simplify(10000, preserve_topology=False)
        if isinstance(raw_poly, Polygon):
            raw_poly = MultiPolygon([raw_poly])
        for poly in list(raw_poly):
            list_of_polygons.append(poly)

    # Create lists of lists.
    x_coords = [list(poly.exterior.coords.xy[0]) for poly in list_of_polygons]
    y_coords = [list(poly.exterior.coords.xy[1]) for poly in list_of_polygons]

    return x_coords, y_coords


class radarPass:
    """Class to hold all relevant data for a given radar pass."""

    def __init__(self, lon, lat, timestamp, height, radar_refl, cpr_cloud_mask):

        # Data for the full pass
        self.cloudsat = {}
        self.cloudsat['lon'] = lon.tolist()
        self.cloudsat['lat'] = lat.tolist()
        self.cloudsat['timestamp'] = timestamp.tolist()
        self.cloudsat['radar_refl'] = radar_refl.tolist()
        self.cloudsat['cpr_cloud_mask'] = cpr_cloud_mask.tolist()
        self.cloudsat['height'] = height.tolist()

        # Dictionaries for dealing with cloudsat cloud classification data
        self.cloudsat['cloud_dict'] = {
            '0000': (0, 'None'),
            '0001': (1, 'Ci'),
            '0010': (2, 'As'),
            '0011': (3, 'Ac'),
            '0100': (4, 'St'),
            '0101': (5, 'Sc'),
            '0110': (6, 'Cu'),
            '0111': (7, 'Ns'),
            '1000': (8, 'Deep')
        }

        self.cloudsat['precip_dict'] = {
            '00': (0, 'no precipitation'),
            '01': (1, 'liquid precipitation'),
            '10': (2, 'solid precipitation'),
            '11': (3, 'possible drizzle')
        }


    def add_cloudclass(self, cldclass_data, cldclass_vdata):
        """Adds CloudSat cloud classification data to the radarPass instance.
        
        PARAMETERS: 
        -----------
        cldclass_data: data dictionary output from 'read_cloudsat_file' function
        applied to a CloudSat cloud classification file. 

        cldclass_vdata: vdata dictionary output from 'read_cloudsat_file' function
        applied to a CloudSat cloud classification file.  

        OUTPUTS: 
        --------
        self (with cloud classification data)

        EXCEPTIONS: 
        -----------
        Exception: "Cloud classification file doesn't match with this radar pass."
        Raised when the first five latitude/longitude values in the 
        radarPass instance and cloud classification file don't match. 
        
        """

        if (sum([int(a[0] == a[1]) for a in zip(self.cloudsat['lon'][:5], cldclass_vdata['Longitude'][:5])]) == 5 and
                sum([int(a[0] == a[1]) for a in zip(self.cloudsat['lat'][:5], cldclass_vdata['Latitude'][:5])]) == 5):

            # All the cloud class and precipitation type information in stored in a binary format.
            # Need to decode the binary to retrieve the cloud class information.
            self.cloudsat['cloud_class'] = np.flipud(
                cldclass_data['cloud_scenario'].T).tolist()

            binary_cc = [["{0:b}".format(val) for val in line] for line in
                         self.cloudsat['cloud_class']]

            self.cloudsat['cloud_type'] = [[self.cloudsat['cloud_dict']
                                            [val[-5:-1]][0] for val in line]
                                           for line in binary_cc]

            self.cloudsat['precip_type'] = [[self.cloudsat['precip_dict']
                                             [val[-14:-12]
                                              ][0] if len(val) > 12 else -1
                                             for val in line] for line in binary_cc]

        else:

            # If the data from the cloud class file is incompatible with the CloudSat data (first five 
            # latitudes and longitudes don't match), set fields equal to None. 
            self.cloudsat['cloud_class'] = None
            self.cloudsat['cloud_type'] = None
            self.cloudsat['precip_type'] = None

            raise Exception(
                "Cloud classification file doesn't match with this radar pass.")

        return self

    def add_cwc(self, cwc_data, cwc_vdata):
        """Adds CloudSat cloud water content data to the radarPass instance.
        
        PARAMETERS: 
        -----------
        cwc_data:  data dictionary output from 'read_cloudsat_file' function
        applied to a CloudSat cloud water content file.

        cwc_vdata: vdata dictionary output from 'read_cloudsat_file' function
        applied to a CloudSat cloud water content file.  

        OUTPUTS: 
        --------
        self (with cloud water content data)

        EXCEPTIONS: 
        -----------
        Exception: "Cloud water content file doesn't match with this radar pass."
        Raised when the first five latitude/longitude values in the 
        radarPass instance and cloud classification file don't match. 
        
        """

        if (sum([int(a[0] == a[1]) for a in zip(self.cloudsat['lon'][:5], cwc_vdata["Longitude"][:5])]) == 5 and
                sum([int(a[0] == a[1]) for a in zip(self.cloudsat['lat'][:5], cwc_vdata["Latitude"][:5])]) == 5):

            self.cloudsat['lwc'] = np.flipud(
                cwc_data['RO_liq_water_content'].T)/1000
            self.cloudsat['lwc_unc'] = np.flipud(
                cwc_data['RO_liq_water_content_uncertainty'].T).astype(float)
            self.cloudsat['iwc'] = np.flipud(
                cwc_data['RO_ice_water_content'].T)/1000
            self.cloudsat['iwc_unc'] = np.flipud(
                cwc_data['RO_ice_water_content_uncertainty'].T).astype(float)
            self.cloudsat['lwc_unc'][self.cloudsat['lwc'] < 0] = np.nan
            self.cloudsat['iwc_unc'][self.cloudsat['iwc'] < 0] = np.nan
            self.cloudsat['lwc'][self.cloudsat['lwc'] < 0] = np.nan
            self.cloudsat['iwc'][self.cloudsat['iwc'] < 0] = np.nan

        else:

            self.cloudsat['lwc'] = None
            self.cloudsat['lwc_unc'] = None
            self.cloudsat['iwc'] = None
            self.cloudsat['iwc_unc'] = None

            raise Exception(
                "Cloud water content file doesn't match with this radar pass.")

        return self

    def trim_pass(self, area_of_interest):
        """Trim the radar pass so that it only contains profiles in the 'area of interest'. 

        PARAMETERS: 
        -----------
        area_of_interest: GeoDataFrame containing polygon that specifies 
        the "area of interest" for the dataframe. CRS must be equal to the 
        default_crs. 

        OUPTUTS: 
        --------
        self (with all CloudSat fields trimmed to the area of interest).
        """

        def determine_trim_indices(radar_pass, area_of_interest):
            """Determine which profile indices will be kept after the trimming."""
            points = [Point(lon, lat) for lon, lat in zip(radar_pass.cloudsat['lon'],
                                                          radar_pass.cloudsat['lat'])]
            traj_gdf = gpd.GeoDataFrame({'idx': list(range(len(radar_pass.cloudsat['lon']))),
                                         'geometry': points})
            traj_gdf.crs = {'init': 'epsg:4326'}
            traj_gdf = traj_gdf.to_crs({'init': default_crs})
            traj_gdf = gpd.sjoin(traj_gdf, area_of_interest, how='inner', op='within')\
                          .drop('index_right', axis=1)
            inds_to_keep = np.asarray(traj_gdf['idx'])
            return inds_to_keep

        inds_to_keep = determine_trim_indices(self, area_of_interest)

        # General radar pass data
        self.cloudsat['lon'] = (np.array(self.cloudsat['lon'])[
                                inds_to_keep]).tolist()
        self.cloudsat['lat'] = (np.array(self.cloudsat['lat'])[
                                inds_to_keep]).tolist()
        self.cloudsat['timestamp'] = (
            np.array(self.cloudsat['timestamp'])[inds_to_keep]).tolist()
        self.cloudsat['radar_refl'] = (np.array(self.cloudsat['radar_refl'])[
                                       :, inds_to_keep]).tolist()
        self.cloudsat['cpr_cloud_mask'] = (np.array(self.cloudsat['cpr_cloud_mask'])[
                                           :, inds_to_keep]).tolist()

        # Cloud mask data
        try:
            self.cloudsat['cloud_class'] = (np.array(self.cloudsat['cloud_class'])[
                                            :, inds_to_keep]).tolist()
            self.cloudsat['cloud_type'] = (np.array(self.cloudsat['cloud_type'])[
                                           :, inds_to_keep]).tolist()
            self.cloudsat['precip_type'] = (np.array(self.cloudsat['precip_type'])[
                                            :, inds_to_keep]).tolist()
        except:
            self.cloudsat['cloud_class'] = None
            self.cloudsat['cloud_type'] = None
            self.cloudsat['precip_type'] = None

        # Cloud water content data
        try:
            self.cloudsat['lwc'] = (np.array(self.cloudsat['lwc'])[
                                    :, inds_to_keep]).tolist()
            self.cloudsat['lwc_unc'] = (np.array(self.cloudsat['lwc_unc'])[
                                        :, inds_to_keep]).tolist()
            self.cloudsat['iwc'] = (np.array(self.cloudsat['iwc'])[
                                    :, inds_to_keep]).tolist()
            self.cloudsat['iwc_unc'] = (np.array(self.cloudsat['iwc_unc'])[
                                        :, inds_to_keep]).tolist()
        except:
            self.cloudsat['lwc'] = None
            self.cloudsat['lwc_unc'] = None
            self.cloudsat['iwc'] = None
            self.cloudsat['iwc_unc'] = None
        self.cloudsat['trim_inds'] = inds_to_keep

        return self

    def add_era5_data(self):
        """Adds ERA5 fields to the radar pass. The date/hour for the ERA5 
        fields is determined by rounding satellite pass' representative
        time to the nearest hour. 

        TERMINOLOGY NOTE: 'atmospheric variables' files are netCDF4 files that contain 
        ERA5 weather variables up in the atmosphere (eg. at 850 mb and 700 mb pressure 
        levels). Their filenames names have 'atmospheric_vbls' in them. Similarly, 
        'single level variables' files are files that contain ERA5 variables that 
        are only available near the surface. These files have 'singlelevel_vbls' 
        in the filename. 

        PARAMETERS: 
        -----------
        None
        
        OUTPUTS:
        --------
        self (with one additional attribute called 'era5'. 'era5' is a dictionary that contains 
        all the era5 data, including the x coordinates, y coordinates, and atmospheric variables)
        """

        def find_era_date_hour(radar_pass):
            """Find the appropriate ERA5 date/hour for the given radarPass instance.
            
            PARAMETERS: 
            -----------
            radar_pass: radarPass instance. 

            OUTPUTS: 
            --------
            era5_date: ERA5 date in string format (%Y%m%d).

            era5_hour: ERA5 hour in integer format. 
            
            """
            t = radar_pass.get_representative_time()
            if t.minute >= 30:
                rounded_time = t.replace(
                    second=0, microsecond=0, minute=0, hour=t.hour+1)
            else:
                rounded_time = t.replace(second=0, microsecond=0, minute=0)
            return datetime.datetime.strftime(rounded_time, '%Y%m%d'), rounded_time.hour

        def open_atmospheric_vbls_file(radar_pass):
            """Opens ERA5 'atmospheric variables' file corresponding to the 
            given radarPass instance.  

            PARAMETERS: 
            -----------
            radar_pass: radarPass instance. 

            OUTPUTS: 
            --------
            atm_dataset: netCDF4 Dataset object for atmospheric variables file. 

            hr: hour corresponding to atmospheric variables file (integer). 
            """

            dstr, hr = find_era_date_hour(radar_pass)
            file = 'ERA5_atmospheric_vbls_'+dstr+'.nc'
            atm_dataset = Dataset(os.path.join(data_path, 'ERA5', file))

            return atm_dataset, hr

        def open_singlelevel_vbls_file(radar_pass):
            """Opens ERA5 'single level variables' file corresponding to the 
            given radarPass instance.  

            PARAMETERS: 
            -----------
            radar_pass: radarPass instance. 

            OUTPUTS: 
            --------
            atm_dataset: netCDF4 Dataset object for single level variables file. 

            hr: hour corresponding to single level variables file (integer). 
            """

            dstr, hr = find_era_date_hour(radar_pass)
            file = 'ERA5_singlelevel_vbls_'+dstr+'.nc'
            slev_dataset = Dataset(os.path.join(data_path, 'ERA5', file))

            return slev_dataset, hr

        def reproject_from_epsg4326(lon, lat, field, dst_crs):
            """Reproject a data field from lat/lon coordinates to a specified crs.
            
            PARAMETERS: 
            -----------
            lon: array-like list/vector of longitude coordinates for the data field.  

            lat: array-like list/vector of latitude coordinates for the data field.  

            field: array-like data field of shape: len(lat), len(lon).  

            dst_crs: destination CRS. 

            OUTPUTS: 
            --------
            dst_x: vector of x coordinates for reprojected field

            dst_y: vector list/vector of y coordinates for reprojected field

            destination_array: array of reprojected field (of shape len(y), len(x)). 
            """

            # Parameters needed for the transformation.
            width = field.shape[1]
            height = field.shape[0]
            left, bottom, right, top = lon[0], lat[-1], lon[-1], lat[0]
            src_crs = {'init': 'epsg:4326'}

            # Calculate affine transformation matrix for the source field.
            src_transform = from_bounds(
                left, bottom, right, top, width, height)

            # Calculate affine transformation matrix, width, and height for the destination field.
            dst_transform, dst_width, dst_height = calculate_default_transform(src_crs, dst_crs, width,
                                                                               height, left=left,
                                                                               bottom=bottom,
                                                                               right=right, top=top)

            # Perform reprojection.
            destination_array = np.zeros((dst_height, dst_width))
            reproject(source=field, destination=destination_array, src_transform=src_transform,
                      src_crs=src_crs, dst_transform=dst_transform,
                      dst_crs=dst_crs, resampling=Resampling.nearest)

            # Replace fill values with NaN's.
            destination_array[destination_array == 1e20] = np.NaN

            # Get X and Y vectors for the transformed field. I could be doing this wrong,
            # but my method definitely works for this specific Arctic case!
            dst_x = np.linspace(dst_transform[2], dst_transform[5], dst_width)
            dst_y = np.linspace(dst_transform[2], dst_transform[5], dst_height)

            return dst_x, dst_y, destination_array

        def get_lat_lon_level_hours(atm_dataset):
            """Extracts latitude, longitude, time, and level vectors 
            from an ERA5 netCDF dataset.
            
            PARAMETERS: 
            -----------
            atm_dataset: netCDF4 Dataset object for a given ERA5 'atmospheric variables'
            or 'single level variables' file. 

            OUTPUTS: 
            --------
            lat: vector of latitude coordinates.  

            lon: vector of longitude coordinates.  

            hours: vector of hours.  

            levels: vector of levels. 'None' for single level variables files. 
            """

            lat = atm_dataset.variables['latitude'][:]
            lon = atm_dataset.variables['longitude'][:]
            hours = np.array([(pd.Timestamp('19000101')+(hrs_since*pd.Timedelta('1 hours'))).hour
                              for hrs_since in atm_dataset.variables['time'][:]])

            if 'level' in atm_dataset.variables.keys():
                levels = atm_dataset.variables['level'][:]
            else:
                levels = None

            return lat, lon, hours, levels

        def get_weather_variable_names(atm_dataset):
            """Retrieves the names of netCDF variables that correspond to weather variables
            (and are not latitude, longitude, level, or time).  

            PARAMETERS: 
            -----------
            atm_dataset: netCDF4 Dataset object for ERA5 'atmospheric variables' or 
            'single level variables' dataset. 

            OUTPUTS: 
            --------
            var_list = list of weather variable names. 
            """

            return [var for var in atm_dataset.variables.keys() if var not in
                    ['longitude', 'latitude', 'level', 'time']]

        def add_era5_atmospheric_vbls(atm_dataset, hr, era5={}):
            """Adds variables from an ERA5 'atmospheric variables' netCDF4 file 
            to an 'era5' dictionary. 

            PARAMETERS: 
            -----------
            atm_dataset: netCDF4 Dataset object corresponding to the 
            'atmospheric variables' file. 

            hr: integer hour corresponding to the radarPass object's 
            representative time. 

            era5: Dictionary object to add the ERA5 data to. 

            OUTPUTS: 
            --------
            era5: Dictionary with added ERA5 data. Contains keys for x 
            coordinates, y coordinates, and data fields for each variable/level. 
            """

            lat, lon, hours, levels = get_lat_lon_level_hours(atm_dataset)
            wx_vbls = get_weather_variable_names(atm_dataset)

            for vbl in wx_vbls:
                for level in levels:
                    i_time = np.where(hours == hr)[0][0]
                    i_level = np.where(levels == level)[0][0]
                    field = atm_dataset.variables[vbl][i_time, i_level, :, :]
                    reproj_x, reproj_y, reproj_field = reproject_from_epsg4326(lon, lat, field,
                                                                               {'init': default_crs})
                    if vbl in ['q', 'clwc', 'ciwc']:  # Change units
                        reproj_field = reproj_field*1000
                    elif vbl == 'z':
                        reproj_field = reproj_field/9.81

                    era5[vbl+'_'+str(level)] = np.flipud(reproj_field).tolist()

            # Convert vertical velocities to mm/s
            era5['w_700'] = (-np.array(era5['w_700'])*287 *
                             np.array(era5['t_700'])/(7e4*9.81)*1000).tolist()
            era5['w_850'] = (-np.array(era5['w_850'])*287 *
                             np.array(era5['t_850'])/(8.5e4*9.81)*1000).tolist()

            era5['x'], era5['y'] = reproj_x, reproj_y

            return era5

        def add_era5_singlelevel_vbls(slev_dataset, hr, era5={}):
            """Adds variables from an ERA5 'single level variables' netCDF4 file 
            to an 'era5' dictionary. 

            PARAMETERS: 
            -----------
            atm_dataset: netCDF4 Dataset object corresponding to the 
            'single level variables' file. 

            hr: integer hour corresponding to the radarPass object's 
            representative time. 

            era5: Dictionary object to add the ERA5 data to. 

            OUTPUTS: 
            --------
            era5: Dictionary with added ERA5 data. Contains keys for x 
            coordinates, y coordinates, and data fields for each added variable. 
            """

            lat, lon, hours, _ = get_lat_lon_level_hours(slev_dataset)
            wx_vbls = get_weather_variable_names(slev_dataset)

            for vbl in wx_vbls:
                i_time = np.where(hours == hr)[0][0]
                field = slev_dataset.variables[vbl][i_time, :, :]
                reproj_x, reproj_y, reproj_field = reproject_from_epsg4326(lon, lat, field,
                                                                           {'init': default_crs})
                if vbl in ['msl', 'sp']:  # Change units
                    reproj_field = reproj_field/100

                era5[vbl] = np.flipud(reproj_field).tolist()

            if 'x' not in era5.keys() and 'y' not in era5.keys():
                era5['x'], era5['y'] = reproj_x.tolist(), reproj_y.tolist()

            return era5

        atm_dataset, hr = open_atmospheric_vbls_file(self)
        era5 = add_era5_atmospheric_vbls(atm_dataset, hr, era5={})
        slev_dataset, hr = open_singlelevel_vbls_file(self)
        era5 = add_era5_singlelevel_vbls(slev_dataset, hr, era5=era5)
        self.era5 = era5

        return self

    def get_json_serializable_obj(self):
        """Converts structure of radarPass object to something that's 
        JSON serializable (contains only dictionaries and lists).
        
        PARAMETERS: 
        -----------
        None

        OUTPUTS: 
        --------
        radar_dict: dictionary containing all the information that the
        radarPass object has. 
        """

        radar_dict = self.__dict__
        for key, val in radar_dict['cloudsat'].items():
            if isinstance(val, np.ndarray):
                radar_dict['cloudsat'][key] = val.tolist()
        for key, val in radar_dict['era5'].items():
            if isinstance(val, np.ndarray):
                radar_dict['era5'][key] = val.tolist()
        return radar_dict

    def reduce_size_cloudsat(self, reduction_factor=3,
                             trim_vbls=['cpr_cloud_mask', 'cloud_class', 'precip_type', 'lwc_unc', 'iwc_unc']):
        """Reduces the size of the CloudSat data fields in two ways: 
        1. Reduces image sizes. 
        2. Removes fields that are less useful. 

        PARAMETERS: 
        -----------
        reduction_factor: Default 3. Image sizes are reduced by taking every n'th pixel in 
        the time dimension, where n is the reduction_factor. 

        trim_vbls: List of less-important variables that will be removed. 

        OUTPUTS: 
        --------
        self (with smaller and/or fewer data fields)
        """

        reduce_1d = ['lon', 'lat', 'timestamp']
        for var in reduce_1d:
            self.cloudsat[var] = (np.array(self.cloudsat[var])[
                                  ::reduction_factor]).tolist()

        reduce_2d = ['radar_refl', 'cpr_cloud_mask', 'cloud_class', 'cloud_type', 'precip_type',
                     'lwc', 'lwc_unc', 'iwc', 'iwc_unc']
        for var in reduce_2d:
            self.cloudsat[var] = (np.array(self.cloudsat[var])[
                                  :, ::reduction_factor]).tolist()

        for vbl in trim_vbls:
            del(self.cloudsat[vbl])

        return self

    def reduce_size_era5(self, reduction_factor=3, trim_vbls=['d2m', 'sp', 'z_850', 'z_700']):
        """Reduces the size of the ERA5 data fields in two ways: 
        1. Reduces image sizes. 
        2. Removes fields that are less useful. 

        PARAMETERS: 
        -----------
        reduction_factor: Default 3. Image sizes are reduced by taking every n'th pixel in 
        the x and y dimensions, where n is the reduction_factor. 

        trim_vbls: List of less-important variables that will be removed. 

        OUTPUTS: 
        --------
        self (with smaller and/or fewer data fields)
        """

        reduce_1d = ['x', 'y']
        for var in reduce_1d:
            self.era5[var] = (np.array(self.era5[var])[
                              ::reduction_factor]).tolist()

        reduce_2d = [var for var in self.era5.keys()
                     if var not in reduce_1d]
        for var in reduce_2d:
            self.era5[var] = (np.array(self.era5[var])[
                              ::reduction_factor, ::reduction_factor]).tolist()

        for vbl in trim_vbls:
            del(self.era5[vbl])

        return self

    def get_representative_time(self):
        """Retrieves a representative time for the radarPass instance. This 
        representative time is the median time for the radar pass. 

        PARAMETERS: 
        -----------
        None
        
        OUTPUTS: 
        --------
        rep_time: representative time (pd.Timestamp)
        """
        return self.cloudsat['timestamp'][len(self.cloudsat['timestamp'])//2]

