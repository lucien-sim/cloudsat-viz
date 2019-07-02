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
from bokeh.palettes import Reds6
from bokeh.plotting import figure, save
from bokeh.layouts import gridplot
from bokeh.models.widgets import Select
from bokeh.models.formatters import DatetimeTickFormatter

app = Flask(__name__)

data_path = os.path.join('.','data')
default_crs = 'epsg:3995'

def get_radar_trajectory(radar_pass):
    """Returns a dataframe with radar trajectory in EPSG:3995"""
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
    """Establishes bounds in the image values. Replaces values above and below those 
    bounds with specified fill values.
    """
    image = np.array(image)
    image[image<=min_val] = fill_min
    image[image>=max_val] = fill_max
    image = image.tolist()
    return image


def get_hex_matplotlib_cmap(cmap_name,num_steps,reverse=False): 
    """Creates bokeh-compatible color maps (lists of hex colors) from standard Matplotlib colormaps."""
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
    countries = load_country_geometries(area_of_interest)
    x_land_coords,y_land_coords = prepare_polygon_coords_for_bokeh(countries)
    land_src = ColumnDataSource({'x':x_land_coords,'y':y_land_coords})
    return land_src


def create_CDS_satellite_trajectory(radar_pass): 
    traj_df = get_radar_trajectory(radar_pass)
    traj_src = ColumnDataSource(traj_df)
    return traj_df,traj_src


def create_CDS_era5_field(radar_pass,era5_vbl): 
    era_img_src = ColumnDataSource({'image':[np.array(radar_pass.era5[era5_vbl]).tolist()]})
    return era_img_src


def create_map_panel(radar_pass, era5_vbl, era_img_src, land_src, traj_src, 
                     era5_color_maps, era5_color_bounds, era5_full_names): 

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
    sat_src = ColumnDataSource({'t':[min(traj_df['t']),min(traj_df['t'])],'y':[-5000,25000]})
    return sat_src


def create_CDS_satellite_image(radar_pass,cloudsat_vbl): 
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


#===================================================================================
# FUNCTIONS THAT WILL PROBABLY END UP IN OTHER FILES

def load_radarPass_object_pkl(pass_timestamp,data_path):
    fname = 'radarPass_plot_'+pass_timestamp+'.pkl'
    with open(os.path.join(data_path,'radar_passes',fname),'rb') as f:
        radar_pass = pickle.load(f)
    return radar_pass


def specify_area_of_interest_EPSG3995(bbox=(-3e6,-3e6,3e6,3e6)): 
    area_of_int = gpd.GeoDataFrame({'geometry':[box(*bbox)]})
    area_of_int.crs = {'init':default_crs}
    return area_of_int


def load_country_geometries(area_of_interest):
    countries = gpd.read_file(os.path.join(data_path,'Countries_WGS84','Countries_WGS84.shp'))
    countries = countries.to_crs({'init':default_crs})
    countries = gpd.overlay(countries,area_of_interest,how='intersection')
    return countries


def prepare_polygon_coords_for_bokeh(countries): 
    """Need to go from a list of points to two lists of lists (one for x and y coordinates). 
    For each list of lists, the inner lists contain the x or y coordinates for each point in
    a single polygon, while the outer list has one element for each polygon. 
    """
    # Simplify shapes (to resolution of 10000 meters), convert polygons to multipolygons. 
    list_of_polygons = []
    for raw_poly in countries['geometry']: 
        raw_poly = raw_poly.simplify(10000, preserve_topology=False)
        if isinstance(raw_poly,Polygon): 
            raw_poly = MultiPolygon([raw_poly])
        for poly in list(raw_poly): 
            list_of_polygons.append(poly)
            
    # Create lists of lists. 
    x_coords = [list(poly.exterior.coords.xy[0]) for poly in list_of_polygons]
    y_coords = [list(poly.exterior.coords.xy[1]) for poly in list_of_polygons]
    
    return x_coords,y_coords

class radarPass: 
    """Class for the radar pass... will hold all the relevant data!"""
    
    def __init__(self,lon,lat,timestamp,height,radar_refl,cpr_cloud_mask):
        
        # Data for the full pass
        self.cloudsat = {}
        self.cloudsat['lon'] = lon.tolist()
        self.cloudsat['lat'] = lat.tolist()
        self.cloudsat['timestamp'] = timestamp.tolist()
        self.cloudsat['radar_refl'] = radar_refl.tolist()
        self.cloudsat['cpr_cloud_mask'] = cpr_cloud_mask.tolist()
        self.cloudsat['height'] = height.tolist()
        
        # For dealing with cloud classification data
        self.cloudsat['cloud_dict'] = {
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
        self.cloudsat['precip_dict'] = {
            '00': (0,'no precipitation'),
            '01': (1,'liquid precipitation'),
            '10': (2,'solid precipitation'), 
            '11': (3,'possible drizzle')
        }
    
    def add_cloudclass(self, cldclass_data, cldclass_vdata): 
        if (sum([int(a[0]==a[1]) for a in zip(self.cloudsat['lon'][:5],cldclass_vdata['Longitude'][:5])])==5 and 
            sum([int(a[0]==a[1]) for a in zip(self.cloudsat['lat'][:5],cldclass_vdata['Latitude'][:5])])==5): 
            
            # All the cloud class and precipitation type information in stored in a binary format. 
            # Need to decode the binary to retrieve the information. 
            self.cloudsat['cloud_class'] = np.flipud(cldclass_data['cloud_scenario'].T).tolist()
            binary_cc = [["{0:b}".format(val) for val in line]  for line in 
                         self.cloudsat['cloud_class']]
            self.cloudsat['cloud_type'] = [[self.cloudsat['cloud_dict']
                                            [val[-5:-1]][0] for val in line] 
                                             for line in binary_cc]
            self.cloudsat['precip_type'] = [[self.cloudsat['precip_dict']
                                            [val[-14:-12]][0] if len(val)>12 else -1 
                                             for val in line] for line in binary_cc]
            
        else: 
            # Set fields equal to None if the lat/lon data from the two files doesn't match. 
            # Means that the two files are incompatible. 
            self.cloudsat['cloud_class'] = None
            self.cloudsat['cloud_type']  = None
            self.cloudsat['precip_type'] = None
            print("Cloud classification file doesn't match up with this radar pass.") 
            
        return self
    
    def add_cwc(self, cwc_data, cwc_vdata): 
        if (sum([int(a[0]==a[1]) for a in zip(self.cloudsat['lon'][:5],cwc_vdata["Longitude"][:5])])==5 and 
            sum([int(a[0]==a[1]) for a in zip(self.cloudsat['lat'][:5],cwc_vdata["Latitude"][:5])])==5): 
            
            self.cloudsat['lwc'] = np.flipud(cwc_data['RO_liq_water_content'].T)/1000
            self.cloudsat['lwc_unc'] = np.flipud(cwc_data['RO_liq_water_content_uncertainty'].T).astype(float)
            self.cloudsat['iwc'] = np.flipud(cwc_data['RO_ice_water_content'].T)/1000
            self.cloudsat['iwc_unc'] = np.flipud(cwc_data['RO_ice_water_content_uncertainty'].T).astype(float)
            self.cloudsat['lwc_unc'][self.cloudsat['lwc']<0] = np.nan
            self.cloudsat['iwc_unc'][self.cloudsat['iwc']<0] = np.nan
            self.cloudsat['lwc'][self.cloudsat['lwc']<0] = np.nan
            self.cloudsat['iwc'][self.cloudsat['iwc']<0] = np.nan
            
        else: 
            self.cloudsat['lwc'] = None
            self.cloudsat['lwc_unc'] = None
            self.cloudsat['iwc'] = None
            self.cloudsat['iwc_unc'] = None
            print("Cloud water content file doesn't match up with this radar pass.") 
            
        return self
    
    def trim_pass(self,area_of_interest): 
        """Trim the radar pass so that it only contains profiles in the 'area of interest' 
        (which is given in the default_crs).
        """
        
        def determine_trim_indices(radar_pass,area_of_interest): 
            """Determine which profile indices will be kept after the trimming."""
            points = [Point(lon,lat) for lon,lat in zip(radar_pass.cloudsat['lon'],
                                                        radar_pass.cloudsat['lat'])]
            traj_gdf = gpd.GeoDataFrame({'idx':list(range(len(radar_pass.cloudsat['lon']))),
                                         'geometry':points})
            traj_gdf.crs = {'init':'epsg:4326'}
            traj_gdf = traj_gdf.to_crs({'init':default_crs})
            traj_gdf = gpd.sjoin(traj_gdf,area_of_interest,how='inner',op='within')\
                          .drop('index_right',axis=1)
            inds_to_keep = np.asarray(traj_gdf['idx'])
            return inds_to_keep
        
        inds_to_keep = determine_trim_indices(self,area_of_interest)
        
        self.cloudsat['lon'] = (np.array(self.cloudsat['lon'])[inds_to_keep]).tolist()
        self.cloudsat['lat'] = (np.array(self.cloudsat['lat'])[inds_to_keep]).tolist()
        self.cloudsat['timestamp'] = (np.array(self.cloudsat['timestamp'])[inds_to_keep]).tolist()
        self.cloudsat['radar_refl'] = (np.array(self.cloudsat['radar_refl'])[:,inds_to_keep]).tolist()
        self.cloudsat['cpr_cloud_mask'] = (np.array(self.cloudsat['cpr_cloud_mask'])[:,inds_to_keep]).tolist()
        # Cloud mask data
        try: 
            self.cloudsat['cloud_class'] = (np.array(self.cloudsat['cloud_class'])[:,inds_to_keep]).tolist()
            self.cloudsat['cloud_type'] = (np.array(self.cloudsat['cloud_type'])[:,inds_to_keep]).tolist()
            self.cloudsat['precip_type'] = (np.array(self.cloudsat['precip_type'])[:,inds_to_keep]).tolist()
        except: 
            self.cloudsat['cloud_class'] = None
            self.cloudsat['cloud_type'] = None
            self.cloudsat['precip_type'] = None
        # Cloud water content data  
        try: 
            self.cloudsat['lwc'] = (np.array(self.cloudsat['lwc'])[:,inds_to_keep]).tolist()
            self.cloudsat['lwc_unc'] = (np.array(self.cloudsat['lwc_unc'])[:,inds_to_keep]).tolist()
            self.cloudsat['iwc'] = (np.array(self.cloudsat['iwc'])[:,inds_to_keep]).tolist()
            self.cloudsat['iwc_unc'] = (np.array(self.cloudsat['iwc_unc'])[:,inds_to_keep]).tolist()
        except: 
            self.cloudsat['lwc'] = None
            self.cloudsat['lwc_unc'] = None
            self.cloudsat['iwc'] = None
            self.cloudsat['iwc_unc'] = None
        self.cloudsat['trim_inds'] = inds_to_keep
        
        return self
    
    def add_era5_data(self):
        """Adds ERA5 fields to the radar pass. The date/hour for those fields is determined by 
        rounding the median time in the satellite's pass through the Arctic to the nearest hour. 
        
        PARAMETERS: 
        self 
        
        RETURNS:
        self, with one additional attribute called 'era5'. 'era5' is a dictionary that contains 
        all the era5 data, including the x coordinates, y coordinates, and atmospheric variables. 
        """
        
        def find_era_date_hour(radar_pass):
            """Identifies the date/hour in the ERA data that's closest to the middle of 
            the satellite's pass through the Arctic.
            """
            t = radar_pass.cloudsat['timestamp'][len(radar_pass.cloudsat['timestamp'])//2]
            if t.minute >= 30:
                rounded_time = t.replace(second=0, microsecond=0, minute=0, hour=t.hour+1)
            else:
                rounded_time = t.replace(second=0, microsecond=0, minute=0)
            return datetime.datetime.strftime(rounded_time,'%Y%m%d'), rounded_time.hour

        
        def open_atmospheric_vbls_file(radar_pass):
            """Opens ERA5 'atmospheric variables' file that corresponds to the given radar pass. 
            Returns the file handle and the hour of data in the file that corresponds 
            to the radar pass. 
            """
            dstr,hr = find_era_date_hour(radar_pass)
            file = 'ERA5_atmospheric_vbls_'+dstr+'.nc'
            atm_dataset = Dataset(os.path.join(data_path,'ERA5',file))

            return atm_dataset,hr

        
        def open_singlelevel_vbls_file(radar_pass): 
            """Opens ERA5 'single level variables' file that corresponds to the given radar pass. 
            Returns the file handle and the hour of data in the file that corresponds 
            to the radar pass. 
            """
            dstr,hr = find_era_date_hour(radar_pass)
            file = 'ERA5_singlelevel_vbls_'+dstr+'.nc'
            slev_dataset = Dataset(os.path.join(data_path,'ERA5',file))

            return slev_dataset,hr
        
        
        def reproject_from_epsg4326(lon,lat,field,dst_crs):
            """Reproject ERA5 field from lat/lon coordinates to specified destination crs."""
            # Parameters needed for the transformation. 
            width = field.shape[1]
            height = field.shape[0]
            left,bottom,right,top = lon[0],lat[-1],lon[-1],lat[0]
            src_crs = {'init':'epsg:4326'}

            # Calculate affine transformation matrix for the source field. 
            src_transform = from_bounds(left, bottom, right, top, width, height)

            # Calculate affine transformation matrix, width, and height for the destination field.
            dst_transform, dst_width, dst_height = calculate_default_transform(src_crs, dst_crs, width, 
                                                                               height, left = left, 
                                                                               bottom = bottom, 
                                                                               right = right, top = top)

            # Perform reprojection. 
            destination_array = np.zeros((dst_height,dst_width))
            reproject(source=field, destination=destination_array, src_transform=src_transform,
                      src_crs=src_crs, dst_transform=dst_transform, 
                      dst_crs=dst_crs, resampling=Resampling.nearest)

            # Replace fill values with NaN's. 
            destination_array[destination_array==1e20] = np.NaN

            # Get X and Y vectors for the transformed field. I could be doing this wrong, 
            # but my method definitely works for this specific Arctic case! 
            dst_x = np.linspace(dst_transform[2],dst_transform[5],dst_width)
            dst_y = np.linspace(dst_transform[2],dst_transform[5],dst_height)

            return dst_x,dst_y,destination_array

        
        def get_lat_lon_level_hours(atm_dataset): 
            """Extracts latitude, longitude, time, and level vectors from netCDF dataset."""
            lat = atm_dataset.variables['latitude'][:]
            lon = atm_dataset.variables['longitude'][:]
            hours = np.array([(pd.Timestamp('19000101')+(hrs_since*pd.Timedelta('1 hours'))).hour 
                              for hrs_since in atm_dataset.variables['time'][:]])
            if 'level' in atm_dataset.variables.keys(): 
                levels = atm_dataset.variables['level'][:]
            else: 
                levels = None
            return lat,lon,hours,levels


        def get_weather_variable_names(atm_dataset):
            """Retrieves names of netCDF variables that correspond to weather variables (and 
            therefore aren't the latitude, longitude, atmospheric level, or time).
            """
            return [var for var in atm_dataset.variables.keys() if var not in 
                   ['longitude', 'latitude', 'level', 'time']]

        
        def add_era5_atmospheric_vbls(atm_dataset,hr,era5={}):
            """Adds ERA5 variables that are only available at atmospheric levels 
            to the 'era5' dictionary.
            """

            lat,lon,hours,levels = get_lat_lon_level_hours(atm_dataset)
            wx_vbls = get_weather_variable_names(atm_dataset)

            for vbl in wx_vbls: 
                for level in levels:
                    i_time = np.where(hours==hr)[0][0]
                    i_level = np.where(levels==level)[0][0]
                    field = atm_dataset.variables[vbl][i_time,i_level,:,:]
                    reproj_x,reproj_y,reproj_field = reproject_from_epsg4326(lon,lat,field,
                                                                             {'init':default_crs})
                    if vbl in ['q','clwc','ciwc']: # Change units 
                        reproj_field = reproj_field*1000
                    elif vbl == 'z': 
                        reproj_field = reproj_field/9.81
                        
                    era5[vbl+'_'+str(level)] = np.flipud(reproj_field).tolist()
            
            # Convert vertical velocities to mm/s
            era5['w_700'] = (-np.array(era5['w_700'])*287*np.array(era5['t_700'])/(7e4*9.81)*1000).tolist()
            era5['w_850'] = (-np.array(era5['w_850'])*287*np.array(era5['t_850'])/(8.5e4*9.81)*1000).tolist()
            
            era5['x'],era5['y'] = reproj_x, reproj_y
            
            return era5

        
        def add_era5_singlelevel_vbls(slev_dataset,hr,era5={}):
            """Adds ERA5 variables that are only available at the surface to the 'era5' 
            dictionary.
            """

            lat,lon,hours,_ = get_lat_lon_level_hours(slev_dataset)
            wx_vbls = get_weather_variable_names(slev_dataset)

            for vbl in wx_vbls: 
                i_time = np.where(hours==hr)[0][0]
                field = slev_dataset.variables[vbl][i_time,:,:]
                reproj_x,reproj_y,reproj_field = reproject_from_epsg4326(lon,lat,field,
                                                                         {'init':default_crs})
                if vbl in ['msl','sp']: # Change units
                    reproj_field = reproj_field/100
                        
                era5[vbl] = np.flipud(reproj_field).tolist()

            if 'x' not in era5.keys() and 'y' not in era5.keys(): 
                era5['x'],era5['y'] = reproj_x.tolist(), reproj_y.tolist()

            return era5

        atm_dataset,hr = open_atmospheric_vbls_file(self)
        era5 = add_era5_atmospheric_vbls(atm_dataset,hr,era5={})
        slev_dataset,hr = open_singlelevel_vbls_file(self)
        era5 = add_era5_singlelevel_vbls(slev_dataset,hr,era5=era5)
        self.era5 = era5
        
        return self
    
    def get_json_serializable_obj(self):
        """Convert structure and all data fields to dictionaries and lists."""
        radar_dict = self.__dict__
        for key,val in radar_dict['cloudsat'].items(): 
            if isinstance(val,np.ndarray): 
                radar_dict['cloudsat'][key] = val.tolist()
        for key,val in radar_dict['era5'].items():
            if isinstance(val,np.ndarray): 
                radar_dict['era5'][key] = val.tolist()    
        return radar_dict
    
    def reduce_size_cloudsat(self, reduction_factor=3, 
                             trim_vbls=['cpr_cloud_mask','cloud_class','precip_type','lwc_unc','iwc_unc']): 
        """Reduce image sizes by taking every reduction_factor'th profile. Also drop some fields that
        are less useful.
        """
        
        reduce_1d = ['lon','lat','timestamp']
        for var in reduce_1d: 
            self.cloudsat[var] = (np.array(self.cloudsat[var])[::reduction_factor]).tolist()

        reduce_2d = ['radar_refl','cpr_cloud_mask','cloud_class','cloud_type','precip_type',
                     'lwc','lwc_unc','iwc','iwc_unc']
        for var in reduce_2d: 
            self.cloudsat[var] = (np.array(self.cloudsat[var])[:,::reduction_factor]).tolist()
            
        for vbl in trim_vbls: 
            del(self.cloudsat[vbl])
            
        return self
    
    def reduce_size_era5(self, reduction_factor=3, trim_vbls=['d2m','sp','z_850','z_700']): 
        """Reduce image sizes by taking every reduction_factor'th pixel in both the x and y directions. 
        Also drop some fields that are less useful.
        """
        reduce_1d = ['x','y']
        for var in reduce_1d: 
            self.era5[var] = (np.array(self.era5[var])[::reduction_factor]).tolist()
        
        reduce_2d = [var for var in radar_pass.era5.keys() if var not in reduce_1d]
        for var in reduce_2d: 
            self.era5[var] = (np.array(self.era5[var])[::reduction_factor,::reduction_factor]).tolist()
        
        for vbl in trim_vbls: 
            del(self.era5[vbl]) 
        
        return self

#===================================================================================


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