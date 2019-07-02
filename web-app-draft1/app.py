#!/usr/bin/python3

from flask import Flask, render_template, request, redirect, send_file
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point, box
import pickle
import os

from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.io import show, output_file, output_notebook
from bokeh.models import ColumnDataSource, HoverTool, CustomJS, Panel, Tabs
from bokeh.models import LinearColorMapper, LogColorMapper, FixedTicker, ColorBar
from bokeh.plotting import figure, save
from bokeh.layouts import gridplot
from bokeh.models.widgets import Select
from bokeh.models.formatters import DatetimeTickFormatter

import sys
sys.path.insert(0, '../scripts')

from global_functions_and_classes import prepare_polygon_coords_for_bokeh
from global_functions_and_classes import specify_area_of_interest_EPSG3995, radarPass
from global_functions_and_classes import load_radarPass_object_pkl, load_country_geometries
from global_variables import data_path, default_crs


app = Flask(__name__)

default_crs = 'epsg:3995'

def get_radar_trajectory(radar_pass):
    """Returns a dataframe with the radar trajectory in the EPSG:3995 CRS"""
    points = [Point(lon,lat) for lon,lat in zip(radar_pass.cloudsat['lon'],radar_pass.cloudsat['lat'])]
    traj_gdf = gpd.GeoDataFrame({'t':radar_pass.cloudsat['timestamp'],'lat':radar_pass.cloudsat['lat'],
                                 'lon':radar_pass.cloudsat['lon'],'geometry':points})
    traj_gdf.crs = {'init':'epsg:4326'}
    traj_gdf = traj_gdf.to_crs({'init':default_crs})
    traj_gdf['x'] = [pt.coords.xy[0][0] for pt in traj_gdf['geometry']]
    traj_gdf['y'] = [pt.coords.xy[1][0] for pt in traj_gdf['geometry']]
    traj_df = pd.DataFrame(traj_gdf.drop('geometry',axis=1)).iloc[::12]
    traj_df['t_str'] = traj_df['t'].dt.strftime('%Y-%m-%d %H:%M:%S')
    return traj_df


def prepare_image_for_plot(image,min_val,fill_min,max_val,fill_max): 
    """Establishes bounds for the image values. Replaces values above and below those 
    bounds with the specified fill values.
    """
    image = np.array(image)
    image[image<=min_val] = fill_min
    image[image>=max_val] = fill_max
    image = image.tolist()
    return image


def get_hex_matplotlib_cmap(cmap_name,num_steps,reverse=False): 
    """Creates bokeh-compatible color maps (lists of hex colors) from 
    standard Matplotlib colormaps.
    
    PARAMETERS: 
    -----------
    cmap_name: Name of the Matplotlib colormap. 

    num_steps: Number of color shades in the colormap. 

    reverse: Boolean, default False. Reverse colormap if True. 

    OUTPUTS: 
    --------
    hex_map: The Bokeh-compatible colormap. List of hex colors. 
    """

    from matplotlib import cm
    from matplotlib.colors import rgb2hex

    cmap = cm.get_cmap(cmap_name, num_steps)
    hex_map = [rgb2hex(cmap(i)[:3]) for i in range(num_steps)]
    if reverse: 
        hex_map = hex_map[::-1]
    
    return hex_map


def find_corresponding_pass(entered_timestamp,data_path):  
    """Returns timestamp of radar pass that occurs closest to the entered time."""
    pass_times_all = sorted([pd.Timestamp(''.join(fname.split('_')[2:4])[:-4]) for fname in 
                             os.listdir(os.path.join(data_path,'radar_passes')) if fname[0]!='.'])
    pass_time = pass_times_all[np.argmin(np.abs(np.array(pass_times_all)-entered_timestamp))]
    pass_timestamp = pass_time.strftime('%Y%m%d_%H%M%S')
    return pass_timestamp


def create_CDS_land(area_of_interest): 
    """Creates ColumnDataSource for Bokeh land polygons in area of interest."""
    countries = load_country_geometries(area_of_interest)
    x_land_coords,y_land_coords = prepare_polygon_coords_for_bokeh(countries)
    land_src = ColumnDataSource({'x':x_land_coords,'y':y_land_coords})
    return land_src


def create_CDS_satellite_trajectory(radar_pass): 
    """Creates ColumnDataSource for satellite trajectory."""
    traj_df = get_radar_trajectory(radar_pass)
    traj_src = ColumnDataSource(traj_df)
    return traj_df,traj_src


def create_CDS_era5_field(radar_pass,era5_vbl): 
    """Creates ColumnDataSource for ERA5 field image."""
    era_img_src = ColumnDataSource({'image':[np.array(radar_pass.era5[era5_vbl]).tolist()]})
    return era_img_src


def create_map_panel(radar_pass, era5_vbl, era_img_src, land_src, traj_src, 
                     era5_color_maps, era5_color_bounds, era5_full_names): 
    """Creates top panel of data visualization figure."""

    # Create color mapper for ERA5 field. 
    color_mapper = LinearColorMapper(palette = era5_color_maps[era5_vbl], 
                                     low = era5_color_bounds[era5_vbl][0], 
                                     high = era5_color_bounds[era5_vbl][1])

    # Put the figure together. 
    pass_time = radar_pass.cloudsat['timestamp'][len(radar_pass.cloudsat['timestamp'])//2]\
                          .strftime('%Y-%m-%d %H:%M:%S')
    p1 = figure(title="ERA5 FIELD: "+era5_full_names[era5_vbl]+". TRAJECTORY TIME: "+
                pass_time+".", toolbar_location="right", 
                plot_width=800, plot_height=550,x_range=(-3e6, 3e6), 
                y_range=(-3e6, 3e6))
    era_im = p1.image('image', source=era_img_src, x=radar_pass.era5['x'][0], 
                      y=radar_pass.era5['y'][0], 
                      dw=radar_pass.era5['x'][-1]-radar_pass.era5['x'][0], 
                      dh=radar_pass.era5['y'][-1]-radar_pass.era5['y'][0], 
                      color_mapper=color_mapper)
    era_im.glyph.color_mapper.nan_color = (1, 1, 1, 0.1)
    pr_land = p1.patches('x', 'y', source=land_src, fill_color='black', fill_alpha=0.0, 
                         line_color="black", line_width=0.8)
    lr_traj = p1.line('x','y',source=traj_src,line_color='black',line_width=2)
    cr_traj = p1.circle('x','y',source=traj_src,fill_color='black',fill_alpha=0.3,
                        line_color='black',line_alpha=0,hover_fill_color='green',
                        hover_fill_alpha=1,hover_line_color='black',hover_line_alpha=1,
                        size=5)

    # Formatting. 
    p1.xaxis.major_label_text_font_size='9pt'
    p1.yaxis.major_label_text_font_size='9pt'

    # Add colorbar.  
    color_bar = ColorBar(color_mapper=color_mapper, bar_line_color='black',
                         major_tick_line_color='black',label_standoff=6, 
                         border_line_color=None, location=(0,0))
    p1.add_layout(color_bar, 'right')
    
    return p1, cr_traj


def create_CDS_satellite_position(traj_df): 
    """Creates ColumnDataSource for the line on the bottom panel that
    shows the satellite position.
    """
    sat_src = ColumnDataSource({'t':[min(traj_df['t']),min(traj_df['t'])],'y':[-5000,25000]})
    return sat_src


def create_CDS_satellite_image(radar_pass,cloudsat_vbl): 
    """Create ColumnDataSource for the satellite observation image."""
    if cloudsat_vbl == 'radar_refl': 
        img = prepare_image_for_plot(radar_pass.cloudsat['radar_refl'].copy(),-2300,np.nan,1500,1500)
    elif cloudsat_vbl == 'cloud_type':
        img = prepare_image_for_plot(radar_pass.cloudsat['cloud_type'].copy(),0,0,10,10)
    elif cloudsat_vbl == 'cpr_cloud_mask': 
        img = prepare_image_for_plot(radar_pass.cloudsat['cpr_cloud_mask'].copy(),0,0,100,100)
    elif cloudsat_vbl == 'iwc': 
        img = prepare_image_for_plot(radar_pass.cloudsat['iwc'].copy(),0,0,4000,4000)
    elif cloudsat_vbl == 'lwc': 
        img = prepare_image_for_plot(radar_pass.cloudsat['lwc'].copy(),0,0,4000,4000)
    elif cloudsat_vbl == 'iwc_unc': 
        img = prepare_image_for_plot(radar_pass.cloudsat['iwc_unc'].copy(),0,0,4000,4000)
    elif cloudsat_vbl == 'lwc_unc':
        img = prepare_image_for_plot(radar_pass.cloudsat['lwc_unc'].copy(),0,0,4000,4000)
    sat_img_src = ColumnDataSource({'image':[img]})
    return sat_img_src


def create_satellite_panel(radar_pass, cloudsat_vbl, sat_img_src, sat_src, 
                           traj_df, radar_full_names):
    """Creates bottom panel for the data visualization tool plot."""
    
    # Create color mapper. 
    if cloudsat_vbl == 'radar_refl': 
        color_mapper = LinearColorMapper(palette='Plasma256',low=-2300,high=1500)
    elif cloudsat_vbl == 'cloud_type':
        color_mapper = LinearColorMapper(palette='Paired9',low=-0.5,high=8.5)
    elif cloudsat_vbl == 'cpr_cloud_mask': 
        color_mapper = LinearColorMapper(palette='Paired11',low=-2,high=42)
    elif cloudsat_vbl == 'iwc': 
        cmp = get_hex_matplotlib_cmap('Blues',100,reverse=False)
        cmp.insert(0,'#FFFFFF')
        color_mapper = LinearColorMapper(palette=cmp,low=0,high=0.25)
    elif cloudsat_vbl == 'lwc': 
        cmp = get_hex_matplotlib_cmap('Greens',100,reverse=False)
        cmp.insert(0,'#FFFFFF')
        color_mapper = LinearColorMapper(palette=cmp,low=0,high=1)
    elif cloudsat_vbl == 'iwc_unc': 
        color_mapper = LinearColorMapper(palette='Viridis256',low=0,high=200)
    elif cloudsat_vbl == 'lwc_unc': 
        color_mapper = LinearColorMapper(palette='Viridis256',low=0,high=200)

    # Put the plot together. 
    p2 = figure(title="CloudSat: "+radar_full_names[cloudsat_vbl], toolbar_location="right",
                plot_width=800, plot_height=250, active_scroll = "wheel_zoom",
                x_range=(min(traj_df['t']), max(traj_df['t'])),y_range=(0, 15000))
    sat_im = p2.image('image', source=sat_img_src, x=min(traj_df['t']), 
                      y=radar_pass.cloudsat['height'][0], dw=max(traj_df['t'])-min(traj_df['t']), 
                      dh=(radar_pass.cloudsat['height'][-1]-radar_pass.cloudsat['height'][0]), 
                      color_mapper=color_mapper)
    sat_im.glyph.color_mapper.nan_color = (0, 0, 0, 1)
    lr_sat = p2.line('t','y',source=sat_src,line_color='gray',line_width=3)

    # Add colorbar. 
    #if cloudsat_vbl == 'radar_refl': 
    #    ticker = FixedTicker(ticks=np.linspace(-2300,1500,11))
    #elif cloudsat_vbl == 'cloud_type': 
    #    ticker = FixedTicker(ticks=[i for i in range(9)])
    #elif cloudsat_vbl == 'cpr_cloud_mask':
    #    ticker = FixedTicker(ticks=[4*i for i in range(11)])
    color_bar = ColorBar(color_mapper=color_mapper, bar_line_color='black',
                         major_tick_line_color='black',label_standoff=8, 
                         border_line_color=None, location=(0,0))
    p2.add_layout(color_bar, 'right')

    # Adjust some formatting. 
    p2.xaxis.axis_label = 'Time (UTC)'
    p2.xaxis.axis_label_text_font_size='11pt'
    p2.xaxis.axis_label_text_font_style='normal'
    p2.xaxis.major_label_text_font_size='9pt'
    p2.yaxis.axis_label = 'Height (m)'
    p2.yaxis.axis_label_text_font_size='11pt'
    p2.yaxis.axis_label_text_font_style='normal'
    p2.yaxis.major_label_text_font_size='9pt'

    # Create datetime x axis
    p2.xaxis.formatter = DatetimeTickFormatter(
        hours=["%H:%M:%S"],
        minutes=["%H:%M:%S"],
        seconds=["%H:%M:%S"],
    )

    return p2, sat_src


def add_hovertool(p1, cr_traj, traj_src, sat_src, traj_df): 
    """Adds a hovertool to the top panel of the data visualization tool plot."""

    # Create the JS callback for vertical line on radar plots. 
    callback_htool = CustomJS(args={'traj_src':traj_src,'sat_src':sat_src}, code="""
        const indices = cb_data.index["1d"].indices[0];

        var data_traj = traj_src.data
        var t_traj = data_traj['t']
        const t_val = t_traj[indices]

        var data_sat = sat_src.data;
        var t_sat = data_sat['t']
        t_sat[0] = t_val
        t_sat[1] = t_val
        sat_src.change.emit(); 
    """)

    # Add the hovertool for the satellite trajectory points on top panel, which are 
    # linked to the vertical line on the bottom panel. 
    htool_mode = ('vline' if max(traj_df['y'])-min(traj_df['y'])<=
                              (max(traj_df['x'])-min(traj_df['x'])) else 'hline')
    tooltips1 = [("lat", "@lat"),("lon", "@lon"),('time','@t_str')]
    p1.add_tools(HoverTool(renderers=[cr_traj],callback=callback_htool,
                           mode=htool_mode,tooltips=tooltips1))
    
    return p1


era5_color_bounds = {'q_700':[0,4],
                  'q_850':[0,6],
                  'z_700':[2500,3300],
                  'z_850':[1000,1700],
                  'ciwc_850':[0,0.1],
                  'ciwc_700':[0,0.1],
                  'clwc_850':[0,1],
                  'clwc_700':[0,1],
                  't_850':[230,316], 
                  't_700':[230,316],
                  'u_850':[-25,25],
                  'u_700':[-25,25],
                  'v_850':[-25,25],
                  'v_700':[-25,25],
                  'w_850':[-50,50],
                  'w_700':[-50,50],
                  'd2m':[230,316], 
                  't2m':[230,316], 
                  'msl':[960,1065], 
                  'sp':[650,1065]}

era5_color_maps = {'q_700':get_hex_matplotlib_cmap('BrBG',15,reverse=False), #r
                  'q_850':get_hex_matplotlib_cmap('BrBG',15,reverse=False), #r
                  'z_700':get_hex_matplotlib_cmap('PuOr',25,reverse=True),
                  'z_850':get_hex_matplotlib_cmap('PuOr',25,reverse=True),
                  'ciwc_850':get_hex_matplotlib_cmap('Blues',11,reverse=False),
                  'ciwc_700':get_hex_matplotlib_cmap('Blues',11,reverse=False),
                  'clwc_850':get_hex_matplotlib_cmap('Greens',11,reverse=False),
                  'clwc_700':get_hex_matplotlib_cmap('Greens',11,reverse=False),
                  't_850':get_hex_matplotlib_cmap('bwr',25,reverse=False), 
                  't_700':get_hex_matplotlib_cmap('bwr',25,reverse=False),
                  'u_850':get_hex_matplotlib_cmap('PRGn',15,reverse=False),
                  'u_700':get_hex_matplotlib_cmap('PRGn',15,reverse=False),
                  'v_850':get_hex_matplotlib_cmap('PRGn',15,reverse=False),
                  'v_700':get_hex_matplotlib_cmap('PRGn',15,reverse=False),
                  'w_850':get_hex_matplotlib_cmap('PRGn',15,reverse=False),
                  'w_700':get_hex_matplotlib_cmap('PRGn',15,reverse=False),
                  'd2m':get_hex_matplotlib_cmap('BrBG',25,reverse=False), 
                  't2m':get_hex_matplotlib_cmap('bwr',25,reverse=False), 
                  'msl':get_hex_matplotlib_cmap('PuOr',50,reverse=True), 
                  'sp':get_hex_matplotlib_cmap('PuOr',50,reverse=True)}

era5_color_maps['ciwc_850'].insert(0,'#FFFFFF')
era5_color_maps['clwc_850'].insert(0,'#FFFFFF')
era5_color_maps['ciwc_700'].insert(0,'#FFFFFF')
era5_color_maps['clwc_700'].insert(0,'#FFFFFF')

era5_full_names = {'q_700':'Specific humidity, 700 mb (g/kg)',
                  'q_850':'Specific humidity, 850 mb (g/kg)',
                  'z_700':'Geopotential height, 700 mb (m)', 
                  'z_850':'Geopotential height, 850 mb (m)',
                  'ciwc_850':'Cloud ice water content, 850 mb (g/kg)',
                  'ciwc_700':'Cloud ice water content, 700 mb (g/kg)',
                  'clwc_850':'Cloud liquid water content, 850 mb (g/kg)',
                  'clwc_700':'Cloud liquid water content, 700 mb (g/kg)',
                  't_850':'Temperature, 850 mb (K)', 
                  't_700':'Temperature, 700 mb (K)',
                  'u_850':'E/W Wind, 850 mb (m/s)',
                  'u_700':'E/W Wind, 700 mb (m/s)',
                  'v_850':'N/S Wind, 850 mb (m/s)',
                  'v_700':'N/S Wind, 700 mb (m/s)',
                  'w_850':'Vertical velocity, 850 mb (mm/s)',
                  'w_700':'Vertical velocity, 700 mb (mm/s)',
                  'd2m':'Dewpoint at surface (K)', 
                  't2m':'Temperature at 2 m (K)', 
                  'msl':'Mean sea level pressure (mb)', 
                  'sp':'Surface pressure (mb)'}

cloud_dict = {
    '0000': (0,'None'), 
    '0001': (1,'Ci'), 
    '0010': (2,'As'), 
    '0011': (3,'Ac'),
    '0100': (4,'St'),
    '0101': (5,'Sc'),
    '0110': (6,'Cu'),
    '0111': (7,'Ns'),
    '1000': (8,'Deep')
}

cloud_type_title_key = ', '.join([str(a[0])+'-'+a[1] for a in 
                                  cloud_dict.values()])
radar_full_names = {'radar_refl':'Radar Reflectivity',
                    'cloud_type':'Cloud Type, Key: '+cloud_type_title_key,
                    'cpr_cloud_mask':'Radar Cloud Mask',
                    'iwc':'Cloud Ice Water Content (g/m3)',
                    'lwc':'Cloud Liquid Water Content (g/m3)',
                    'iwc_unc':'Cloud Ice Water Content Uncertainty (%)',
                    'lwc_unc':'Cloud Liquid Water Content Uncertainty (%)'}


@app.route('/')
def index():
    
    # Input variables
    query_date = '2015-09-01'
    query_hour = '00'
    era5_vbl = 'q_850'
    cloudsat_vbl = 'radar_refl'
    entered_time = pd.Timestamp(query_date)+(int(query_hour)*pd.Timedelta('1 hours'))

    # Load the corrsponding radarPass instance. 
    pass_timestamp = find_corresponding_pass(entered_time,data_path)
    radar_pass = load_radarPass_object_pkl(pass_timestamp,data_path)

    # Top panel: map of Arctic with ERA5 field.
    area_of_interest = specify_area_of_interest_EPSG3995()
    land_src = create_CDS_land(area_of_interest)
    traj_df,traj_src = create_CDS_satellite_trajectory(radar_pass)
    era_img_src = create_CDS_era5_field(radar_pass,era5_vbl)
    p1, cr_traj = create_map_panel(radar_pass, era5_vbl, era_img_src, land_src, traj_src, 
                                   era5_color_maps, era5_color_bounds, era5_full_names)

    # Bottom panel: radar profiles. 
    sat_src = create_CDS_satellite_position(traj_df)
    sat_img_src = create_CDS_satellite_image(radar_pass,cloudsat_vbl)
    p2, sat_src = create_satellite_panel(radar_pass, cloudsat_vbl, sat_img_src, sat_src, 
                                         traj_df, radar_full_names) 
    p1 = add_hovertool(p1, cr_traj, traj_src, sat_src, traj_df)

    layout = gridplot([[p1],[p2]])
    
    script_plt, div_plt = components(layout)
    
    return render_template('index.html',bokeh_resources=CDN.render(),
                           date_val=query_date,hour_val=query_hour,
                           era5_val=era5_vbl,cloudsat_val=cloudsat_vbl,
                           div_plt=div_plt,script_plt=script_plt)


@app.route('/submit_query',methods=['POST'])
def submit_query(): 
    
    # Extract data from form
    query_date = request.form.get('query_date')
    query_hour = request.form.get('query_hour')
    era5_vbl = request.form.get('era5_field')
    cloudsat_vbl = request.form.get('sat_field')
    entered_time = pd.Timestamp(query_date)+(int(query_hour)*pd.Timedelta('1 hours'))
    
    # Load the corrsponding radarPass instance. 
    pass_timestamp = find_corresponding_pass(entered_time,data_path)
    radar_pass = load_radarPass_object_pkl(pass_timestamp,data_path)

    # Top panel: map of Arctic with ERA5 field. 
    area_of_interest = specify_area_of_interest_EPSG3995()
    land_src = create_CDS_land(area_of_interest)
    traj_df,traj_src = create_CDS_satellite_trajectory(radar_pass)
    era_img_src = create_CDS_era5_field(radar_pass,era5_vbl)
    p1, cr_traj = create_map_panel(radar_pass, era5_vbl, era_img_src, land_src, traj_src, 
                                   era5_color_maps, era5_color_bounds, era5_full_names)

    # Bottom panel: radar profiles. 
    sat_src = create_CDS_satellite_position(traj_df)
    sat_img_src = create_CDS_satellite_image(radar_pass,cloudsat_vbl)
    p2, sat_src = create_satellite_panel(radar_pass, cloudsat_vbl, sat_img_src, sat_src, 
                                         traj_df, radar_full_names) 
    p1 = add_hovertool(p1, cr_traj, traj_src, sat_src, traj_df)

    layout = gridplot([[p1],[p2]])
    
    script_plt, div_plt = components(layout)
    
    return render_template('index.html',bokeh_resources=CDN.render(),
                           date_val=query_date,hour_val=query_hour,
                           era5_val=era5_vbl,cloudsat_val=cloudsat_vbl,
                           div_plt=div_plt,script_plt=script_plt)


@app.route('/documentation')
def documentation():
    return render_template('documentation.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port,debug=True,threaded=True) 
