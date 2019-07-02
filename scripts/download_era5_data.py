#!/usr/bin/python3

import cdsapi
import os

def download_era5_data(year, month, days, hours, data_path):
    """Downloads ERA5 data for the specified days and hours in a given year/month.

    PARAMETERS: 
    -----------
    year:  Year of ERA5 data (four-character string)
    month: Month of ERA5 data (two-character string)
    days:  List of days in year/month for which ERA5 data is desired (list of 
           two-character strings)
    hours: List of hours for which ERA5 data is desired (list of two-character 
           strings)
    data_path: Path to main directory for data storage. Specified in 
               global_variables.py. 

    OUTPUTS: 
    --------
    None. 

    This function will alwayd download the same variables these variables can be 
    adjusted inside the function. 
    """

    c = cdsapi.Client()

    for day in days:

        atmos_fname = os.path.join(
            data_path, "ERA5", "ERA5_atmospheric_vbls_"+year+month+day+".nc")
        slvl_fname = os.path.join(
            data_path, "ERA5", "ERA5_singlelevel_vbls_"+year+month+day+".nc")

        c.retrieve("reanalysis-era5-pressure-levels", {
            "product_type":   "reanalysis",
            "format":         "netcdf",
            "area":           "90.00/-180.00/60.00/180.00",
            "variable":       ['geopotential', 'specific_cloud_ice_water_content',
                               'specific_cloud_liquid_water_content', 'specific_humidity',
                               'temperature', 'u_component_of_wind',
                               'v_component_of_wind', 'vertical_velocity'],
            "pressure_level": ["850", "700"],
            "year":           year,
            "month":          month,
            "day":            day,
            "time":           hours
        }, atmos_fname)

        c.retrieve('reanalysis-era5-single-levels', {
            'product_type':   'reanalysis',
            'format':         'netcdf',
            'area':           '90.00/-180.00/60.00/180.00',
            'variable':       ['2m_dewpoint_temperature', '2m_temperature',
                               'mean_sea_level_pressure', 'surface_pressure'],
            'year':           year,
            'month':          month,
            'day':            day,
            'time':           hours
        }, slvl_fname)

    return None


if __name__ == '__main__': 

    from global_variables import data_path

    year = "2015"
    month = "09"
    days = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11",
            "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
            "23", "24", "25", "26", "27", "28", "29", "30"]
    hours = ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11",
            "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
            "23"]

    download_era5_data(year, month, days, hours, data_path)
