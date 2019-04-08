# -*- coding: utf-8 -*-
"""
Firt code created on Tue Aug 14 14:59:13 2018

@author: Rizzieri Pedruzzi - pedruzzi@email.unc.edu

CARS Functions to read the input files from input folder:
    activity_data
    emissions_factor
    fleet_mix
    temporal_profile
This functions will read the input data
"""
import os, sys, time, matplotlib, glob, fnmatch
import xarray as xr
import pandas as pd
import geopandas as gpd
import numpy as np
import datetime as dt
from shapely.geometry import Polygon
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import geopandas.tools
import shapely
from shapely.geometry import *

total_start_time = time.time()
case_name = 'Seoul_1x1km'
home_dir   = 'C:/Users/pedruzzi/OneDrive - University of North Carolina at Chapel Hill/0000_EUA_IE/001_mobile_source_inventory/CARS_source_code'
src_dir    = home_dir+'/src'
input_dir  = home_dir+'/input'
inter_dir  = home_dir+'/intermediate'
output_dir = home_dir+'/output_seoul'

if not os.path.exists(home_dir+'/intermediate'):
    print('*** intermediate directory does not exist: Creating a new one')
    os.makedirs(home_dir+'/intermediate/')

#fleet_mix_file          = 'age_distribution.csv'
fleet_mix_Age_file      = 'Seoul_age_distribution_2015.csv' #''age_distribution_Activity_Data_Seoul_2017.csv' #'
Emis_Factor_list        = ['gasoline.csv','diesel.csv','cng.csv','lpg.csv']
#EF_Gasoline_file        = 'gasoline.csv'
#EF_Diesel_file          = 'diesel.csv'
#EF_CNG_file             = 'cng.csv'
#EF_LPG_file             = 'lpg.csv'

avg_SPD_Dist_file      = 'avgSpeedDistribution_rev_00.csv'
ambient_temp = 15.0
plot_24 = 'no '#'yes'
link_shape              = '/shapes'+'/GIS Korea Road Link Activity Data'+ \
                          '/shape_seoul'+'/seoul_eup_links_Avg_VKT_UTM52N.shp'
                          
link_shape_att          = ['link_id'  , 'EMD_CD' , 'EMD_ENG_NM', 'EMD_KOR_NM',
                           'road_type', 'speed', 'length_2', 'Avg_VKT']
#
#link_shape              = '/shapes/GIS Korea Road Link Activity Data'+ \
#                          '/shape_seoul/seoul_eup_road_by_county_UTM52N.shp'
#link_shape_att          = ['LINK_ID'  , 'EMD_CD' , 'EMD_ENG_NM', 'ROAD_NAME',
#                           'ROAD_RANK', 'MAX_SPD', 'length_2', 'Avg_VKT']

county_shape              = '/shapes/GIS Korea Road Link Activity Data'+ \
                          '/shape_seoul/seoul_eup_UTM52N.shp'


#link_shape              = '/shapes/GIS Korea Road Link Activity Data'+ \
#                          '/shape_soul_gyeonggi/soul_gyeonggi_Links_UTM52N.shp'
#link_shape_att          = ['link_id'  , 'EMD_CD' , 'EMD_ENG_NM', 'EMD_ENG_NM',
#                           'link_type', 'speed', 'length',]
#
#county_shape              = '/shapes/GIS Korea Road Link Activity Data'+ \
#                          '/shape_soul_gyeonggi/soul_gyeonggi_eup_UTM52N.shp'

                          
temporal_profile_folder = input_dir+'/temporal_profile'
temporal_profile_file   = 'temporal_profile_SK.csv'
temporal_cross_ref_file = 'temporal_cross_ref_SK.csv'

temp = np.asarray([.0094, .0060, .0050, .0065, .0126, .0347, .0591, .0605, .0558, .0545,
        .0536, .0532, .0538, .0539, .0559, .0569, .0580, .0611, .0586, .0525, .0495,.0419, .0308, .0162])

grid_size = 1000

activity_file           =  'seoul_2017.csv' #'seoul_gyeonggi_AD.csv' #'seoul_all.csv' #  #'activity_data.csv'

class EmissionFactor_table:
    def __init__(self, dataframe, name ):
        self.dataframe = dataframe
        self.name      = name.split('.')[0]

class Activity_Data_table:
    def __init__(self, dataframe, vhc_name, years):
        self.data     = dataframe
        self.fullname = vhc_name
        self.years    = years

class Roads_Grid_table:
    def __init__(self, grid_dataframe, surrogate, roads_df):
        self.grid      = grid_dataframe
        self.surrogate = surrogate
        self.roads_df  = roads_df

class EF_Grid_table:
    def __init__(self, EF_dataframe, EF_years, VHC_fullname_EF, Fuel_EF, Polls_EF):
        self.data        = EF_dataframe
        self.EF_years    = EF_years
        self.EF_fullname = VHC_fullname_EF
        self.EF_fuels    = Fuel_EF
        self.EF_polls    = Polls_EF

class EF_Speed_Distribution:
    def __init__(self, SPD_dataframe, Speeds, Speed_Bins):
        self.data      = SPD_dataframe
        self.spd       = Speeds
        self.spd_bins  = Speed_Bins

class Emissions_table:
    def __init__(self, County_Emissions, Road_Emissions, County_Emissions_GeoRef, County, Years):
        self.county_emis = County_Emissions
        self.road_emis   = Road_Emissions
        self.county_geo  = County_Emissions_GeoRef
        self.county      = County
        self.years       = Years


# =============================================================================
# Function to read the Fleet Mix
# =============================================================================
def read_flet_mix(input_dir, fleet_mix_file, sep = ';'):
    #sep = ';'
    name = '{0}{1}{2}'.format(input_dir,'/age_distribuition/',fleet_mix_file)
    if os.path.exists(name) == True:
        print ('')
        print ('Reading Fleet Mix ...')
        print (name)
        fleet_mix = pd.read_csv(name, sep = sep).fillna(0)
#        fleet_mix['Name'] = fleet_mix.Vehicle.str.cat(vhc_type[['Types','Fuel']], sep=' ')
    else:
        print ('')
        print('*** ERROR ABORT ***:  Fleet Mix file "" ', fleet_mix_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read fleet mix file')

    return fleet_mix
# =============================================================================

FM = read_flet_mix(input_dir,fleet_mix_Age_file)
year_FM  = [int(i) for i in list(FM.columns[4:])]
vhc_type = FM.loc[:,['Vehicle','Types','Fuel']]
vhc_type['Name'] = vhc_type.Vehicle.str.cat(vhc_type[['Types','Fuel']], sep=' ')




# =============================================================================
# Function to read the temporal profile
# =============================================================================
def read_temporal_profile(input_dir, temporal_profile_file, sep = ';'):
    name = '{0}{1}{2}'.format(input_dir,'/temporal_profile/',temporal_profile_file)
    if os.path.exists(name) == True:
        print ('')
        print ('Reading Temporal profile ...')
        print (name)
        TP = pd.read_csv(name, sep = sep).fillna(np.nan)
    else:
        print ('')
        print('*** ERROR ABORT ***:  Temporal Profile ', temporal_profile_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read Temporal profile')

    return TP

TP = read_temporal_profile(input_dir,temporal_profile_file)

# =============================================================================
# Function to read the temporal profile cross reference
# =============================================================================
def read_temporal_Cross_ref(input_dir, temporal_cross_ref_file, sep = ';'):
    name = '{0}{1}{2}'.format(input_dir,'/temporal_profile/',temporal_cross_ref_file)
    if os.path.exists(name) == True:
        print ('')
        print ('Reading Temporal Cross reference ...')
        print (name)
        TPCR = pd.read_csv(name, sep = sep).fillna(np.nan)

    else:
        print ('')
        print('*** ERROR ABORT ***:  Temporal Profile cross reference ', temporal_cross_ref_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read Temporal Profile cross reference')

    return TPCR

TPCR = read_temporal_Cross_ref(input_dir, temporal_cross_ref_file)



# =============================================================================
# Function to read the Speed average distribution
# =============================================================================
def read_avgSpeedDistribution(input_dir, avg_Speed_Distribution_file, sep = ';'):
    start_time = time.time()
    input_dir = input_dir
    spd_file = 'avgSpeedDistribution_rev_00.csv' #avg_Speed_Distribution_file #ef_file #['gasoline.csv'] #
    sep = ';'
    final_df = pd.DataFrame()
    name = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',spd_file)
    if os.path.exists(name) == True:
        print ('')
        print ('Reading Average Speed Distribution table ...')
        print (name)

        ASD = pd.read_csv(name, sep = sep).fillna(np.nan)

    else:
        print ('')
        print('*** ERROR ABORT ***:  Temporal Profile cross reference ', temporal_cross_ref_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read Temporal Profile cross reference')

    out_spd_bins = pd.DataFrame({'spd_bins': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]})
    out_spd      = pd.DataFrame({'spd': list(ASD.Speed)})

    return EF_Speed_Distribution(ASD, out_spd, out_spd_bins)

avgSpeedDist = read_avgSpeedDistribution(input_dir, avg_SPD_Dist_file, sep = ';')


# =============================================================================
# Function to read the Activity Data
# =============================================================================
def read_activity_data_csv_SK(input_dir, ad_file, sep = ';'):
    start_time = time.time()
    ad_file = ad_file
#    ad_file = activity_file
#    sep = ','
    name = '{0}{1}{2}'.format(input_dir,'/activity_data/',ad_file)
    if os.path.exists(name) == True:
        print ('')
        print ('Reading Activity Data table ...')
        print (name)
        activity_data = (pd.read_csv(name, sep = sep, encoding = 'utf-8')).fillna(0)
        activity_data.loc[:,'Vehicle'] = activity_data.loc[:,'Vehicle'].str.lower()
        activity_data.loc[:,'Fuel']    = activity_data.loc[:,'Fuel'].str.lower()
        activity_data.loc[:,'Types']   = activity_data.loc[:,'Types'].str.lower()
        activity_data = activity_data.rename(columns={'Region_code' : 'region_cd'})
        activity_data.loc[:,'Manufacture_date']   = (activity_data.loc[:,'Manufacture_date'] / 10000).astype(int)
        activity_data['FullName'] = activity_data.Vehicle.str.cat(activity_data[['Types','Fuel']], sep='_')

##        gasoline_correction = (activity_data.loc[activity_data.Fuel == 'gasoline','Daily_VKT']) * 3 #mutiplying the gasoline VKT by 3 as BH asked
##        activity_data.loc[activity_data.Fuel == 'gasoline','Daily_VKT'] = gasoline_correction
        
        sedan_correction = (activity_data.loc[(activity_data.Vehicle == 'sedan'),'Daily_VKT']) * 3
        activity_data.loc[(activity_data.Vehicle == 'sedan'),'Daily_VKT'] = sedan_correction
        suv_correction = (activity_data.loc[(activity_data.Vehicle == 'suv'), 'Daily_VKT']) * 3
        activity_data.loc[(activity_data.Vehicle == 'suv'), 'Daily_VKT'] = suv_correction
        
        grouped = activity_data.groupby(['Manufacture_date','region_cd','FullName']).sum()
        grouped = grouped.unstack().fillna(0)
        grouped.columns = [x[1] for x in grouped.columns]
        grouped.reset_index(inplace=True)
        grouped.columns = grouped.columns.str.lower()
        grouped = grouped.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
        vhc_names = pd.DataFrame({'vhc_name'  : list(activity_data.FullName.unique())})
        vhc_years = pd.DataFrame({'vhc_years' : list(grouped.manufacture_date.unique())})
        out_table = Activity_Data_table(grouped, vhc_names,vhc_years)
    else:
        print('*** ERROR ABORT ***:  Emissions Factor file "" ', ad_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read Emissions Factor file')
    run_time = ((time.time() - start_time))
    print("--- %f seconds ---" % (run_time))
    print("--- %f minutes ---" % (run_time/60))
    print("--- %f Hours   ---" % (run_time/3600))
    return out_table
## =============================================================================

AD_SK = read_activity_data_csv_SK(input_dir, activity_file, sep = ',')


# =============================================================================
# Function to read link level shapefile
# =============================================================================
def roads_grid_surrogate_inf(input_dir, file_name, Link_ID_attr, 
                         Region_Code, Region_Name,
                         RD_name_attr, RD_type_attr, 
                         Speed_attr, Link_length, VKT_attr, Unit_meters = True):
    start_time = time.time()
    
    Link_ID_attr       = Link_ID_attr
    Region_CD          = Region_Code
    Region_NM          = Region_Name
    RD_name_attr       = RD_name_attr 
    RD_type_attr       = RD_type_attr
    Speed_attr         = Speed_attr
    Link_length        = Link_length
    VKT_attr           = VKT_attr
    file_name          = link_shape

    shp_file = '{0}{1}'.format(input_dir,file_name)
    if os.path.exists(shp_file) == True:
        print ('')
        print ('Reading Link Shapefile ...')
        print (shp_file)
        prj_file = shp_file.replace('.shp', '.prj')
        prj = [l.strip() for l in open(prj_file,'r')][0]
        lnk_shp = gpd.read_file(shp_file)
        out_roads = lnk_shp.loc[:,['geometry',Link_ID_attr, Region_CD, Region_NM, 
                                   RD_name_attr, RD_type_attr, Speed_attr, Link_length, VKT_attr]]
        
        number_links = np.arange(0,len(out_roads))
        # changing the name of columns to keep a standard
        out_roads = out_roads.rename(columns={Link_ID_attr       : 'link_id'})
        out_roads = out_roads.rename(columns={Region_CD          : 'region_cd'})
        out_roads = out_roads.rename(columns={Region_NM          : 'region_nm'})
        out_roads = out_roads.rename(columns={RD_name_attr       : 'road_name'})
        out_roads = out_roads.rename(columns={RD_type_attr       : 'road_type'})
        out_roads = out_roads.rename(columns={Speed_attr         : 'max_speed'})
        out_roads = out_roads.rename(columns={Link_length        : 'link_length'})
        out_roads = out_roads.rename(columns={VKT_attr           : 'vkt_avg'})
        out_roads['number_links']  = number_links
    
        out_roads['activity_data'] = (out_roads['link_length'] * 0.0).astype(float)
        out_roads['region_cd']     = out_roads['region_cd'].astype(int)
        out_roads['road_type']     = out_roads['road_type'].astype(int)
        out_roads['link_id']       = out_roads['link_id'].astype(np.int64)
        out_roads['max_speed']     = out_roads['max_speed'].astype(float)
        out_roads['link_length']   = out_roads['link_length'].astype(float)
        out_roads['number_links']  = out_roads['number_links'].astype(int)
        out_roads['geometry_BKP']  = out_roads.geometry
        out_roads['geometry']      = out_roads.buffer(0.1)
        out_roads['total_area']    = out_roads.area
        out_roads['link_split_total']    = (out_roads['link_length'] * 0.0).astype(float)
        out_roads['link_split_county']    = (out_roads['link_length'] * 0.0).astype(float)
        out_roads['vkt_split_county']    = (out_roads['link_length'] * 0.0).astype(float)
        reduc = 0.6                                  #60% reduction as BH asked
        rt = {101 : 80 *reduc, 102 : 60 *reduc, 103 : 60 *reduc, 104 : 50 *reduc,   #60% reduction as BH asked
              105 : 30 *reduc, 106 : 30 *reduc, 107 : 30 *reduc, 108 : 30 *reduc}
        
        for igeocd in out_roads.region_cd.unique():
            aux_split_county = out_roads.link_length.loc[out_roads.region_cd   == igeocd].values / \
                 (out_roads.link_length.loc[out_roads.region_cd  == igeocd]).sum()
            out_roads.loc[out_roads.region_cd == igeocd, ['link_split_county']] = aux_split_county
            
            vkt_split_county = out_roads.vkt_avg.loc[out_roads.region_cd   == igeocd].values / \
                 (out_roads.vkt_avg.loc[out_roads.region_cd  == igeocd]).sum()
            out_roads.loc[out_roads.region_cd == igeocd, ['vkt_split_county']] = vkt_split_county

            
        aux_split_total = out_roads.link_length.values / out_roads.link_length.sum()
        out_roads.loc[:,['link_split_total']] = aux_split_total
        for key, ispd in rt.items():
            out_roads.loc[(out_roads.loc[:,'road_type'] == key),'max_speed'] = ispd
        # Creating grid
        roads_bounds = out_roads.bounds
        xmin = roads_bounds.minx.min() - grid_size
        xmax = roads_bounds.maxx.max() + grid_size
        ymin = roads_bounds.miny.min() - grid_size
        ymax = roads_bounds.maxy.max() + grid_size
        cols =  int(abs(xmax - xmin) / grid_size)
        rows =  int(abs(ymax - ymin) / grid_size)
        print(cols,rows)
        lat_list = (np.arange(ymin,ymax,grid_size))
        lon_list = (np.arange(xmin,xmax,grid_size))
        polygons = []
        grid_row = []
        grid_col = []
        for j in range(0,rows):
            yini = lat_list[j]
            yfin = lat_list[j+1]
            for i in range(0,cols):
                grid_row.append(j+1)
                grid_col.append(i+1)
                xini = lon_list[i]
                xfin = lon_list[i+1]
                polygons.append(Polygon([(xini, yini), (xfin, yini), (xfin, yfin), (xini, yfin), (xini, yini)])) 
                
        crs = out_roads.crs #{'init': 'epsg:4326'}
        grid_ID = [x for x in range (1,len(polygons)+1)]
        grid = gpd.GeoDataFrame({'geometry':polygons, 'grid_id':grid_ID,
                                 'row':grid_row, 'col':grid_col}, crs=crs)
        grid.crs = crs
        # exporting grid as shapefile
        grid.to_file(filename = output_dir+'/grid_{0}.shp'.format(case_name), driver='ESRI Shapefile',crs_wkt=prj)
        #creating the surrogate
        surrogate = gpd.overlay(out_roads, grid, how='intersection').reset_index(drop=True)
        surrogate['split_area']    = surrogate.area
        surrogate['weight_factor'] = surrogate.area / surrogate.total_area #total_area attibute is comming from roads dataframe
        surrogate = surrogate.loc[:,['geometry', 'link_id', 'region_cd',
                                     'link_length', 'total_area', 'split_area',
                                     'grid_id','row', 'col', 'weight_factor']]
        surrogate.to_file(filename = output_dir+'/road_grid_surrogate_{0}.shp'.format(case_name), driver='ESRI Shapefile',crs_wkt=prj)

        out_roads = out_roads.drop(columns=['geometry'])
        out_roads = out_roads.rename(columns={'geometry_BKP': 'geometry'}).set_geometry('geometry')
    else:
        print('*** ERROR ABORT ***:  Shapefile "" ', shp_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read link Shapefile file')
    
    
    run_time = ((time.time() - start_time))
    print("--- %f seconds ---" % (run_time))
    print("--- %f minutes ---" % (run_time/60))
    print("--- %f Hours   ---" % (run_time/3600))

    return Roads_Grid_table(grid, surrogate, out_roads)
# =============================================================================
    
roads_RGS = roads_grid_surrogate_inf(input_dir,link_shape, link_shape_att[0],
                                    link_shape_att[1],link_shape_att[2],
                                    link_shape_att[3],link_shape_att[4],
                                    link_shape_att[5],link_shape_att[6], link_shape_att[7], Unit_meters = True ) 


# =============================================================================
# Function to read link level shapefile
# =============================================================================
def processing_County_shape(input_dir, file_name, Region_CD, Region_name_attr, 
                         Region_name_attr_SK):
    start_time = time.time()
    Region_Geocode       = Region_CD
    Region_name_attr     = Region_name_attr 
    Region_name_attr_SK  = Region_name_attr_SK
#    Link_ID_attr       = 'LINK_ID'
#    RD_name_attr       = 'ROAD_NAME'
#    RD_type_attr       = 'ROAD_RANK'
#    Activity_data_attr = 'SHAPE_STLe'
#    Speed_attr         = 'MAX_SPD'
#    Link_length        = 'SHAPE_STLe'
#    file_name          = link_shape
    cnty_file = '{0}{1}'.format(input_dir,file_name)
    if os.path.exists(cnty_file) == True:
        print ('')
        print ('Reading Link Shapefile to get bounds to create the GRID...')
        print (cnty_file)
        cnty_shp = gpd.read_file(cnty_file)
        cnty_shp = cnty_shp.rename(columns={Region_Geocode       : 'region_cd'})
        cnty_shp = cnty_shp.rename(columns={Region_name_attr     : 'region_nm'})
        cnty_shp = cnty_shp.rename(columns={Region_name_attr_SK  : 'region_nm_SK'})
    
        out_cnty = cnty_shp.loc[:,['geometry','region_cd','region_nm', 
                                   'region_name_SK']]
        
        out_cnty['region_cd'] = out_cnty['region_cd'].astype(int)
    else:
        print('*** ERROR ABORT ***:  Shapefile "" ', cnty_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read county Shapefile file')         
    
    run_time = ((time.time() - start_time))
    print("--- %f seconds ---" % (run_time))
    print("--- %f minutes ---" % (run_time/60))
    print("--- %f Hours   ---" % (run_time/3600))

    return out_cnty
# =============================================================================
    
county_df = processing_County_shape(input_dir, county_shape, 'EMD_CD', 'EMD_ENG_NM',
                             'EMD_KOR_NM') 




# =============================================================================
# Function to read the emissions factor from South Korea* (*particular case)
# =============================================================================
def read_emissions_factor_SK(input_dir, EmisFactor_list, sep = ';'):
    start_time = time.time()
    input_dir = input_dir
    ef_list = EmisFactor_list #ef_file #['gasoline.csv'] #
#    sep = ';'
    final_EF = pd.DataFrame()
    for ifile in range(0,len(ef_list)):
        print(ef_list[ifile])
        name = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',ef_list[ifile])
        if os.path.exists(name) == True:
            print ('')
            print ('Reading Emissions Factor table ...')
            print (name)
            emissions_factor = (pd.read_csv(name, sep = sep)).fillna(0)
            if (final_EF.shape[0] == 0) and (final_EF.shape[1] == 0):
                print('TEY')
                final_EF = pd.DataFrame(columns=emissions_factor.columns)
            iaux = 1
            drop_list = []
            aux = True
            while aux == True:
                pos = (emissions_factor.shape[0] - iaux)
                aux = 0 in list(emissions_factor.loc[pos,['Vehicle','Types','SCC','Pollutant']])
                if aux == True:
                    drop_list.append(pos)
                iaux +=1
            emissions_factor = emissions_factor.drop(index=drop_list)
            emissions_factor['Years'] =  emissions_factor['Years'] / 10000
            emissions_factor['Years'] =  emissions_factor.Years.astype(int)
            emissions_factor.loc[:,'Vehicle'] = emissions_factor.loc[:,'Vehicle'].str.lower()
            emissions_factor.loc[:,'Types'] = emissions_factor.loc[:,'Types'].str.lower()
            emissions_factor.loc[:,'Fuel'] = emissions_factor.loc[:,'Fuel'].str.lower()
            emissions_factor['FullName'] = emissions_factor.Vehicle.str.cat(emissions_factor[['Types','Fuel']], sep='_')

            vhc  = emissions_factor.FullName.unique()
            poll = list(emissions_factor.Pollutant.unique())
            yr   = list(np.sort(emissions_factor.Years.unique()))
    #        ivhc = 'sedan_supercompact_diesel'
    #        ipoll = 'CO'
    #        iyr = 1990
            for ipoll in poll:
                for ivhc in vhc:
                    for iyr in yr:
                        
                        dat_ef = emissions_factor[((emissions_factor.Vehicle   == ivhc.split('_')[0])  & 
                                                   (emissions_factor.Types     == ivhc.split('_')[1])  &
                                                   (emissions_factor.Fuel      == ivhc.split('_')[2])  &
                                                   (emissions_factor.Pollutant == ipoll) &
                                                   (emissions_factor.Years     == iyr))] 
                        dat_ef = dat_ef.reset_index(drop=True)
                        iscc = (emissions_factor.SCC[(emissions_factor.FullName == ivhc)].unique())
                        if (dat_ef.shape[0] == 0) and  (iscc.size != 0) :
                            print('*** WARNING *** There is no emissions factor of year {0} for vehicle {1}'.format(iyr,ivhc))
                        elif (2 <= dat_ef.Years.shape[0] <= 3) and (dat_ef.V[0] == 0) and (dat_ef.Temperatures[0] == 0):
                            df = dat_ef.iloc[[0]]
                            final_EF = final_EF.append(df,ignore_index=True)
                        elif (dat_ef.Years.shape[0] >= 4) and (dat_ef.V[0] != 0) and (dat_ef.Temperatures[0] == 0):
                            df = dat_ef.iloc[[0,1]]
                            final_EF = final_EF.append(df,ignore_index=True)
                        else:
                            df = dat_ef.loc[:]
                            final_EF = final_EF.append(df,ignore_index=True)
    
        else:
            print (' ')
            print('*** ERROR ABORT ***:  Emissions Factor file "" ', ef_list[ifile], ' "" does not exist!')
            sys.exit('CARS preProcessor can not read Emissions Factor file')
    
    final_EF = final_EF.reset_index(inplace=False, drop=True)
    EF_names = pd.DataFrame({'EF_fullname': list(final_EF.FullName.unique())})
    EF_years = pd.DataFrame({'EF_years'   : list(np.sort(final_EF.Years.unique()))})
    EF_fuels = pd.DataFrame({'EF_fuel'    : list(final_EF.Fuel.unique())})
    EF_polls = pd.DataFrame({'EF_poll'    : list(final_EF.Pollutant.unique())})
    run_time = ((time.time() - start_time))
    print("--- %f seconds ---" % (run_time))
    print("--- %f minutes ---" % (run_time/60))
    print("--- %f Hours   ---" % (run_time/3600))
    
    return EF_Grid_table(final_EF, EF_years, EF_names, EF_fuels, EF_polls)
# =============================================================================

EF_All_fuel_function = read_emissions_factor_SK(input_dir, Emis_Factor_list, sep = ';')


                              
def calculate_EmisFact(Input_dir, Emissions_Factor_Table ):
    EF_All_fuel = Emissions_Factor_Table.data.copy()
    input_dir = Input_dir
    start_time = time.time()
#    nspd = np.asarray([x for x in range(1,151)])
#    Spd Bins    :   1, 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,  14,  15, 16
#    Spd Bins mph: 2.5, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,  65,  70, 75
#    Spd Bins kmh:   4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 89, 97, 105, 113, 121
    nspd = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 89, 97, 105, 113, 121]
    EF_spd_limit = []
    for i in list(EF_All_fuel.V.unique()):
        if (type(i) == int) and (i != 0):
            print('*** WARNING *** Check your Emissions Factor input at EF temp depended ')
        elif type(i) != int:
            EF_spd_limit.append(i)
    EF_spd_limit = np.unique(EF_spd_limit)

    dat_EF = pd.DataFrame(columns=['FullName','Pollutant','V','Years','Temperatures','spd','emis_fact'])
    for ispd in nspd:
        aux_EF = EF_All_fuel.loc[:,['FullName','Pollutant','V','Years','Temperatures']]
        ef = (EF_All_fuel.a * (ispd**(EF_All_fuel.b))) + \
                                     (EF_All_fuel.c * ispd**(EF_All_fuel.d)) + \
                                     ((EF_All_fuel.f) * \
                                     (EF_All_fuel.k))
        
        aux_EF['spd']       = (aux_EF.Years * 0) + ispd
        aux_EF['emis_fact'] = ef
        
        dat_EF = dat_EF.append(aux_EF, ignore_index=True, sort=False)
    
    dat_EF['ambient_temp'] = ((dat_EF.spd * 0 ) + ambient_temp).astype(int)
    
    #Extracting the NOx Diesel EF from final EF table
    Temp_EF_Diesel = dat_EF[(dat_EF.Temperatures != 0) & 
                            (dat_EF.Pollutant == 'NOx')]
    Temp_EF_Diesel_drop_list = list(Temp_EF_Diesel.index)
    dat_EF = dat_EF.drop(index=Temp_EF_Diesel_drop_list).reset_index(drop=True)
    Temp_EF_Diesel = Temp_EF_Diesel.reset_index(drop=True)
    
    #Filtering the NOx Diesel EF based on the ambient temperature
    filter_temp_diesel = []
    for tlim in Temp_EF_Diesel.Temperatures.unique():
        if len(tlim) <= 4:
            if (tlim[:2] == 'GT') and (ambient_temp < int(tlim[2:])):
                idx = list(Temp_EF_Diesel[(Temp_EF_Diesel.Temperatures  == tlim)].index)
                filter_temp_diesel.extend(idx)
            elif (tlim[:2] == 'LE') and (ambient_temp > int(tlim[2:])):
                idx = list(Temp_EF_Diesel[(Temp_EF_Diesel.Temperatures  == tlim)].index)
                filter_temp_diesel.extend(idx)
        elif len(tlim) > 4:
            tlim2 = tlim.split('_')
            if (ambient_temp < int(tlim2[0][2:])) or (ambient_temp > int(tlim2[1][2:])):
                idx = list(Temp_EF_Diesel[(Temp_EF_Diesel.Temperatures  == tlim)].index)
                filter_temp_diesel.extend(idx)
                
    Temp_EF_Diesel = Temp_EF_Diesel.drop(index=filter_temp_diesel).reset_index(drop=True)
    dat_EF = dat_EF.append(Temp_EF_Diesel, ignore_index=True, sort=True)
    
    #The EF were calculated for all speeds for all equations - Filtering the EF
    # based of the equation and speed 
    drop_list = []
    for vlim in EF_spd_limit:
        if vlim[:2] == 'GT':
            idx = list(dat_EF[(dat_EF.V   == vlim) & 
                              (dat_EF.spd <= int(vlim[2:]))].index)
            drop_list.extend(idx)
        elif vlim[:2] == 'LE':
            idx = list(dat_EF[(dat_EF.V   == vlim) & 
                              (dat_EF.spd > int(vlim[2:]))].index)
            drop_list.extend(idx)
    dat_EF = dat_EF.drop(index=drop_list).reset_index(drop=True)
    
    #Filtering the duplicate EF of NOx Diesel - we are averaging them
    Temp_EF_Diesel = dat_EF[(dat_EF.Temperatures != 0) & 
                            (dat_EF.Pollutant == 'NOx')]
    Temp_EF_Diesel_drop_list = list(Temp_EF_Diesel.index)
    dat_EF = dat_EF.drop(index=Temp_EF_Diesel_drop_list).reset_index(drop=True)
    Temp_EF_Diesel = Temp_EF_Diesel.groupby(['FullName','Pollutant','Temperatures','V','Years','ambient_temp','spd']).mean()
    Temp_EF_Diesel = Temp_EF_Diesel.reset_index(drop=False)
    dat_EF = dat_EF.append(Temp_EF_Diesel, ignore_index=True, sort=True)
    dat_EF = dat_EF.sort_values(by=['FullName','spd'])
    dat_EF = dat_EF.reset_index(drop=True)
    dat_EF.loc[(dat_EF.emis_fact < 0),'emis_fact'] = 0
    EF_names = pd.DataFrame({'EF_fullname': list(np.sort(dat_EF.FullName.unique()))})
    EF_years = pd.DataFrame({'EF_years'   : list(np.sort(dat_EF.Years.unique()))})
    EF_fuels = pd.DataFrame({'EF_fuel'    : list(EF_All_fuel.Fuel.unique())})
    EF_polls = pd.DataFrame({'EF_poll'    : list(np.sort(dat_EF.Pollutant.unique()))})
    
    EmisFact_yr_spd = EF_Grid_table(dat_EF, EF_years, EF_names, EF_fuels, EF_polls)
    run_time = ((time.time() - start_time))
    print("--- %f seconds ---" % (run_time))
    print("--- %f minutes ---" % (run_time/60))
    print("--- %f Hours   ---" % (run_time/3600))
    
    return EmisFact_yr_spd

EmisFactor_yr_spd = calculate_EmisFact(input_dir, EF_All_fuel_function)



def calc_County_Emissions(Input_dir, Emissions_Factor_DataFrame, 
                          Activity_Data_DataFrame, Roads_DataFrame, County_DataFrame ):
    print('*****************************************')
    print('*** Starting calculation of emissions ***')
    print('***          Please wait ...          ***')
    print('*****************************************')
    start_time = time.time()
    input_dir = Input_dir          #Input_dir
    EF_yr_spd  = EmisFactor_yr_spd #Emissions_Factor_DataFrame
    AD_yr_vhc  = AD_SK             #Activity_Data_DataFrame
    roads_DF   = roads_RGS         #Roads_DataFrame
#    county_df  = county_df         #County_DataFrame
    
  
    temp_EF_df = EF_yr_spd.data.copy()
    nspd = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 89, 97, 105, 113, 121]
    EFyears = list(np.sort(temp_EF_df.Years.unique()))
    ADyears = list(np.sort(AD_yr_vhc.data.manufacture_date.unique()))
    EF_yr_min = np.min(EFyears)
    AD_yr_max = np.max(ADyears)
    nyears = list(np.arange(EF_yr_min,AD_yr_max+1))
    final_EF_df = pd.DataFrame(columns=['Years','Pollutant','spd']+list(temp_EF_df.FullName.unique()))
    yrlist = []
    spdlist = []
    for jspd in nspd:
        for jyr in nyears:
            spdlist.append(jspd)
            yrlist.append(jyr)
    for ipol in list(temp_EF_df.Pollutant.unique()):
        inter_EF = pd.DataFrame({'Years' : yrlist,
                              'spd'   : spdlist,
                              'Pollutant' : [ipol for i in range(0,len(yrlist))]})
        for ivhc in list(temp_EF_df.FullName.unique()):
            aux_df = temp_EF_df.loc[(temp_EF_df.FullName == ivhc)  & 
                                    (temp_EF_df.Pollutant == ipol), ['Years','emis_fact','spd']].reset_index(drop=True)
            aux_df = aux_df.rename(columns={'emis_fact':ivhc})
            aux_df.Years = aux_df.Years.astype(int)
            aux_df.spd = aux_df.spd.astype(int)
            inter_EF = pd.merge(inter_EF,aux_df, left_on = ['Years','spd'], right_on = ['Years','spd'], how='left') #.fillna(method='bfill')
            #inter_EF = inter_EF.fillna(method='ffill')
        final_EF_df = final_EF_df.append(inter_EF,ignore_index=True, sort=False)
        for ispd in nspd:
            final_EF_df.loc[(final_EF_df.Pollutant == ipol) & (final_EF_df.spd == ispd)] = \
            final_EF_df.loc[(final_EF_df.Pollutant == ipol) & (final_EF_df.spd == ispd)].fillna(method='bfill').fillna(method='ffill')
    
    final_EF_df = final_EF_df.fillna(0)
    
    
    aux_emis_AD = AD_yr_vhc.data.copy()
    EF_yr_min = np.min(nyears)
    AD_years = list(np.sort(AD_yr_vhc.data.manufacture_date.unique()))
    AD_yrlist = []
    region_list = []
    for iregion in list(np.sort(AD_yr_vhc.data.region_cd.unique())):
        for zyr in AD_years:
            AD_yrlist.append(zyr)
            region_list.append(iregion)
    
    final_AD_df = pd.DataFrame({'manufacture_date' : AD_yrlist,
                                'region_cd'        : region_list})
    aux_emis_AD.loc[:,['manufacture_date','region_cd']] = aux_emis_AD.loc[:,['manufacture_date','region_cd']].astype(int)
    final_AD_df = pd.merge(final_AD_df,aux_emis_AD, left_on = ['manufacture_date','region_cd'], 
                        right_on = ['manufacture_date','region_cd'], how='left').fillna(0.0)
    
    
    inter_AD = final_AD_df.loc[final_AD_df.manufacture_date <= EF_yr_min]
    drop_list_AD = list(inter_AD.index)
    inter_AD = inter_AD.groupby(['region_cd']).sum()
    inter_AD = inter_AD.reset_index(drop=False)
    inter_AD.loc[:,'manufacture_date'] = EF_yr_min
    
    final_AD_df = final_AD_df.drop(index=drop_list_AD).reset_index(drop=True)
    final_AD_df = final_AD_df.append(inter_AD,ignore_index=True, sort=False)
    final_AD_df = final_AD_df.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
    
    Emissions_by_yr_county = pd.DataFrame(columns=['pollutant'] + list(final_AD_df.columns))
    #jpol = 'CO'
    for jpol in list(EF_yr_spd.data.Pollutant.unique()):
        aux_emis_AD = final_AD_df.copy()
        aux_emis_AD['pollutant'] = [jpol for i in aux_emis_AD.index]
        Emissions_by_yr_county = Emissions_by_yr_county.append(aux_emis_AD,ignore_index=True, sort=False)
    Emissions_by_yr_county = Emissions_by_yr_county.sort_values(by=['region_cd','pollutant','manufacture_date']).reset_index(drop=True)
    Emissions_by_yr_county.loc[:,list(AD_yr_vhc.fullname.vhc_name)] = 0
    
    aux_emis_AD = final_AD_df.copy()
    
    for iroad in np.sort(roads_DF.roads_df.road_type.unique()): #[101,103,104,107]: #  
        aux_emis_EF = final_EF_df.copy()
        aux_RD = roads_DF.roads_df.loc[(roads_DF.roads_df.road_type == iroad),:].reset_index(drop=True)
        road_weight = aux_RD.groupby(['region_cd']).sum().reset_index(drop=False)
        road_weight = road_weight.loc[:,['region_cd','link_split_county']]
        aux_avgSpd = avgSpeedDist.data.loc[:,['Speed',str(iroad)]]
        for ispd in avgSpeedDist.data.Speed:
            aux_emis_EF.loc[(aux_emis_EF.spd.isin([ispd])),list(temp_EF_df.FullName.unique())] = \
            aux_emis_EF.loc[(aux_emis_EF.spd.isin([ispd])),list(temp_EF_df.FullName.unique())] * \
            np.float(aux_avgSpd.loc[(aux_avgSpd.Speed.isin([ispd])),str(iroad)])
        
        weight_EF = pd.DataFrame(columns=aux_emis_EF.columns)
        for iyr in aux_emis_EF.Years.unique():
            aux_EF = aux_emis_EF.loc[aux_emis_EF.Years ==  iyr]
            aux_EF = aux_EF.groupby(['Pollutant','Years']).sum()
            aux_EF = aux_EF.reset_index(drop=False)
            weight_EF = weight_EF.append(aux_EF,ignore_index=True, sort=False)
        weight_EF = weight_EF.sort_values(by=['Pollutant','Years']).reset_index(drop=True)
        weight_EF = weight_EF.drop(columns=['spd'])
        weight_EF.loc[:,'Years'] = weight_EF.loc[:,'Years'].astype(int)
        
        for ipol in list(EF_yr_spd.data.Pollutant.unique()):
            nregion = list(road_weight.region_cd.unique())
            aux_AD_county = aux_emis_AD.loc[(aux_emis_AD.region_cd.isin(nregion))].reset_index(drop=True)
            aux_AD_county.loc[:,'manufacture_date'] = aux_AD_county.loc[:,'manufacture_date'].astype(int)
            apply_RD_weight = pd.DataFrame({'region_cd' : aux_AD_county.region_cd})
            apply_RD_weight = pd.merge(apply_RD_weight,road_weight, left_on = 'region_cd', 
                              right_on = ['region_cd'], how='left')
            apply_EF = pd.DataFrame({'Years'     : list(aux_AD_county.manufacture_date),
                                     'Pollutant' : [ipol for i in aux_AD_county.index]})
            apply_EF = pd.merge(apply_EF,weight_EF, left_on = ['Years','Pollutant'], 
                        right_on = ['Years','Pollutant'], how='left')
            nvhc = list(np.sort(AD_yr_vhc.fullname.vhc_name))
            year_list_EF = list(aux_AD_county.manufacture_date.unique())
            aux_AD_county.loc[:,nvhc] = aux_AD_county.loc[:,nvhc].apply(lambda x: np.asarray(x) * apply_RD_weight.link_split_county.values)
            aux_AD_county.loc[:,nvhc] = aux_AD_county.loc[:,nvhc] * apply_EF.loc[:,nvhc]
    
            Emissions_by_yr_county.loc[(Emissions_by_yr_county.region_cd.isin(nregion)) &
                (Emissions_by_yr_county.manufacture_date.isin(year_list_EF)) &
                (Emissions_by_yr_county.pollutant.isin([ipol])), nvhc] = \
            (Emissions_by_yr_county.loc[(Emissions_by_yr_county.region_cd.isin(nregion)) &
                (Emissions_by_yr_county.manufacture_date.isin(year_list_EF)) &
                (Emissions_by_yr_county.pollutant.isin([ipol])), nvhc].reset_index(drop=True) + \
                aux_AD_county.loc[:,nvhc].reset_index(drop=True)).values
    
    nvhc = list(np.sort(AD_yr_vhc.fullname.vhc_name))
    total_county = Emissions_by_yr_county.groupby(['region_cd','pollutant']).sum()
    total_county = total_county.reset_index(drop=False)
    Emissions_by_link = pd.merge(roads_DF.roads_df,total_county, left_on = ['region_cd'], 
                        right_on = ['region_cd'], how='left')
    
    Emissions_by_link.loc[:,nvhc] = Emissions_by_link.loc[:,nvhc].apply(lambda x: np.asarray(x) * \
                         Emissions_by_link.vkt_split_county.values)    

    #Calculation the total emissions for county and link level
    Emissions_by_yr_county['total_emis'] = Emissions_by_yr_county.loc[:,nvhc].sum(axis=1)
    Emissions_by_link['total_emis']      = Emissions_by_link.loc[:,nvhc].sum(axis=1)
    total_county['total_emis']           = total_county.loc[:,nvhc].sum(axis=1)
    county_table = pd.DataFrame({'region_cd' : Emissions_by_yr_county.region_cd.unique()})
    years_table  = pd.DataFrame({'years' : Emissions_by_yr_county.manufacture_date.unique()})

#    cnt_info = county_df.loc[:,['geometry','region_cd']]
    total_county = pd.merge(total_county,county_df, left_on = ['region_cd'], 
                        right_on = ['region_cd'], how='left')
    total_county = gpd.GeoDataFrame(total_county, crs=county_df.crs, geometry='geometry')

    print('*****************************************')
    print('*** Emissions calculation is done     ***')
    print('*****************************************')        
    run_time = ((time.time() - start_time))
    print("--- %f seconds ---" % (run_time))
    print("--- %f minutes ---" % (run_time/60))
    print("--- %f Hours   ---" % (run_time/3600))
        
    return Emissions_table(Emissions_by_yr_county, Emissions_by_link, 
                           total_county, county_table, years_table)

County_Emissions = calc_County_Emissions(input_dir, EmisFactor_yr_spd,
                                       AD_SK, roads_RGS, county_df)


#cnt_info = county_df.loc[:,['geometry','region_cd']]
#total_county = pd.merge(total_county,cnt_info, left_on = ['region_cd'], 
#                    right_on = ['region_cd'], how='left')
#
#
#
#
#aaa = County_Emissions.road_emis.groupby(['region_cd','pollutant']).sum()
#aaa = aaa.reset_index(drop=False)
#bbb = aaa.loc[(aaa.region_cd == 11260105) & (aaa.pollutant == 'CO')]
#
#
#
ccc = County_Emissions.road_emis.loc[(County_Emissions.road_emis.pollutant == 'CO')]
ccc.plot(column='total_emis', cmap='jet', linewidth=0.8)

ddd = County_Emissions.county_geo.loc[(County_Emissions.county_geo.pollutant == 'CO')]
ddd.plot(column='total_emis', cmap='jet')

