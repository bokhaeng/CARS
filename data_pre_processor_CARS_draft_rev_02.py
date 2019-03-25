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
#link_shape              = '/shapes'+'/GIS Korea Road Link Activity Data'+ \
#                          '/ITS Standard Road Links and Nodes'+'/SK_Link_WGS84_CG.shp'

link_shape              = '/shapes/GIS Korea Road Link Activity Data'+ \
                          '/shape_seoul/seoul_eup_road_by_county_UTM52N.shp'
link_shape_att          = ['LINK_ID'  , 'EMD_CD' , 'EMD_ENG_NM', 'ROAD_NAME',
                           'ROAD_RANK', 'MAX_SPD', 'length_2',]

county_shape              = '/shapes/GIS Korea Road Link Activity Data'+ \
                          '/shape_seoul/seoul_eup_UTM52N.shp'
                          
temporal_profile_folder = input_dir+'/temporal_profile'
temporal_profile_file   = 'temporal_profile_SK.csv'
temporal_cross_ref_file = 'temporal_cross_ref_SK.csv'

temp = np.asarray([.0094, .0060, .0050, .0065, .0126, .0347, .0591, .0605, .0558, .0545,
        .0536, .0532, .0538, .0539, .0559, .0569, .0580, .0611, .0586, .0525, .0495,.0419, .0308, .0162])

grid_size = 1000

activity_file           = 'seoul_2017.csv' #'seoul_all.csv' #  #'activity_data.csv'

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
        activity_data = (pd.read_csv(name, sep = sep)).fillna(0)
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
                         Speed_attr, Link_length, Unit_meters = True):
    start_time = time.time()
    
    Link_ID_attr       = Link_ID_attr
    Region_CD          = Region_Code
    Region_NM          = Region_Name
    RD_name_attr       = RD_name_attr 
    RD_type_attr       = RD_type_attr
    Speed_attr         = Speed_attr
    Link_length        = Link_length
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
                                   RD_name_attr, RD_type_attr, Speed_attr, Link_length]]
        
        number_links = np.arange(0,len(out_roads))
        # changing the name of columns to keep a standard
        out_roads = out_roads.rename(columns={Link_ID_attr       : 'link_id'})
        out_roads = out_roads.rename(columns={Region_CD          : 'region_cd'})
        out_roads = out_roads.rename(columns={Region_NM          : 'region_nm'})
        out_roads = out_roads.rename(columns={RD_name_attr       : 'road_name'})
        out_roads = out_roads.rename(columns={RD_type_attr       : 'road_type'})
        out_roads = out_roads.rename(columns={Speed_attr         : 'max_speed'})
        out_roads = out_roads.rename(columns={Link_length        : 'link_length'})
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
        reduc = 0.6                                  #60% reduction as BH asked
        rt = {101 : 80 *reduc, 102 : 60 *reduc, 103 : 60 *reduc, 104 : 50 *reduc,   #60% reduction as BH asked
              105 : 30 *reduc, 106 : 30 *reduc, 107 : 30 *reduc, 108 : 30 *reduc}
        
        for igeocd in out_roads.region_cd.unique():
            aux_split_county = out_roads.link_length.loc[out_roads.region_cd   == igeocd].values / \
                 (out_roads.link_length.loc[out_roads.region_cd  == igeocd]).sum()
            out_roads.loc[out_roads.region_cd == igeocd, ['link_split_county']] = aux_split_county
            
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
                                    link_shape_att[5],link_shape_att[6], Unit_meters = True ) 


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


#pollutants = EF_All_fuel.Pollutant.unique() #['NOx', 'CO', 'VOC', 'PM10', 'PM2.5']

#iyr  =2011

act_data_by_region = AD_SK.data.drop(columns='manufacture_date')
act_data_by_region = act_data_by_region.groupby(['region_cd']).sum()
act_data_by_region.reset_index(inplace=True)
act_data_by_link =  pd.concat([roads_RGS.roads_df, 
                               pd.DataFrame(np.zeros((roads_RGS.roads_df.shape[0],
                                                      AD_SK.fullname.shape[0])),columns=AD_SK.fullname.vhc_name.unique())], axis=1)


for igeocd in AD_SK.data.region_cd.unique():
    if igeocd in list(act_data_by_link.region_cd.unique()):
        for ivhc in AD_SK.fullname.vhc_name:
            if ivhc in list(act_data_by_link.columns):
                aux = act_data_by_link.total_area.loc[act_data_by_link.region_cd   == igeocd] / \
                     (act_data_by_link.total_area.loc[act_data_by_link.region_cd  == igeocd]).sum() * \
                      np.asarray(act_data_by_region[ivhc].loc[act_data_by_region.region_cd == igeocd])
                act_data_by_link.loc[act_data_by_link.region_cd == igeocd,[ivhc]] = aux
            else:
                print('*** WARNING ***:  There is no activity data for vehicle {0} for County {1}'.format(ivhc,igeocd))
            
    else:
        print('*** WARNING *** Region code {0} is not is not in link shapefile.'.format(igeocd))
        print('*** WARNING *** The activity data of region code {0} will be not count into emissions'.format(igeocd))
        print('*** WARNING *** Please, correct your link shapefile to avoid wrong calculations!!')

act_data_age_dist = pd.DataFrame(columns= AD_SK.data.columns)
for ifips in AD_SK.data.region_cd.unique():
    aux_year   = AD_SK.data.manufacture_date[(AD_SK.data.region_cd == ifips)]
    aux_region = AD_SK.data.region_cd[(AD_SK.data.region_cd == ifips)]
    aux_data   = (AD_SK.data[(AD_SK.data.region_cd == ifips)]) / AD_SK.data[(AD_SK.data.region_cd == ifips)].sum()
    aux_data.manufacture_date = aux_year
    aux_data.region_cd      = aux_region
    aux_data = aux_data.fillna(0.0)
    act_data_age_dist = act_data_age_dist.append(aux_data,ignore_index=True)


def calculate_EmisFact_OLD(Input_dir, Emissions_Factor_Table):
    EF_All_fuel = Emissions_Factor_Table
    input_dir = Input_dir
    start_time = time.time()
    array_col =['year','pollutant','spd']+list(EF_All_fuel.EF_fullname.EF_fullname)
    nspd = np.asarray([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                       65, 70, 80, 90, 100, 110, 120, 130, 140, 150])
    yr_loop = np.asarray([iyr for iyr in range(EF_All_fuel.EF_years.EF_years.min(),EF_All_fuel.EF_years.EF_years.max()+1)])
    pol_loop = np.asarray( EF_All_fuel.EF_polls.EF_poll)
    yr_list = []
    spd_list = []
    pol_list = []
    for z in pol_loop:
        for i in yr_loop:
            for j in nspd:
                pol_list.append(z)
                yr_list.append(i)
                spd_list.append(j)
    aux_array = np.zeros(((nspd.shape[0] * \
                          EF_All_fuel.EF_polls.EF_poll.shape[0] * \
                          EF_All_fuel.EF_years.EF_years.shape[0]), len(array_col)))
    
    EmisFact_yr_spd = pd.DataFrame(aux_array, columns = array_col )
    EmisFact_yr_spd.year      = yr_list
    EmisFact_yr_spd.spd       = spd_list
    EmisFact_yr_spd.pollutant = pol_list
    for ipoll in EF_All_fuel.EF_polls.EF_poll:
        print('*** Processing emissions factor for {0}'.format(ipoll))
        for ivhc in EF_All_fuel.EF_fullname.EF_fullname:
            spd = {'spd' : nspd,
                   ivhc  : (np.zeros(nspd.shape[0]))}
            spd_EF = pd.DataFrame(spd)
            yr_vhc_list = pd.DataFrame({'year' : EF_All_fuel.data.Years.loc[(EF_All_fuel.data.FullName == ivhc) &
                                                                            (EF_All_fuel.data.Pollutant == ipoll)].unique()})
            yr_vhc_min = yr_vhc_list.year.min()
            yr_vhc_max = yr_vhc_list.year.max()
    
            for iyr in range(EF_All_fuel.EF_years.EF_years.min(),EF_All_fuel.EF_years.EF_years.max()+1):
                EF_real_yr = []
                if iyr in list(yr_vhc_list.year):
                    iyr_EF = iyr
                elif iyr >= yr_vhc_max:
                    iyr_EF = yr_vhc_list.year[(yr_vhc_list.year-iyr).abs().argsort()[:2]].max()
                else: 
                    iyr_EF = yr_vhc_list.year[(yr_vhc_list.year-iyr).abs().argsort()[:2]].min()
                dat_ef = EF_All_fuel.data[((EF_All_fuel.data.Vehicle  == ivhc.split('_')[0]) & 
                                          (EF_All_fuel.data.Types     == ivhc.split('_')[1]) &
                                          (EF_All_fuel.data.Fuel      == ivhc.split('_')[2]) &
                                          (EF_All_fuel.data.Pollutant == ipoll)              &
                                          (EF_All_fuel.data.Years     == iyr_EF))]
                dat_ef = dat_ef.reset_index(drop=True)
                ef_list = []
                if dat_ef.V.shape[0] > 1 and (dat_ef.V[0] != 0):
                    if   'L' in list(dat_ef.V[0].split('E')):
                        spd_LEGT = float((dat_ef.V[0].split('E'))[1])
                        aux_LE = dat_ef.V[0]
                        aux_GT = dat_ef.V[1]
                        idx0 = (dat_ef.V[dat_ef.V == aux_LE].index)[0]
                        idx1 = (dat_ef.V[dat_ef.V == aux_GT].index)[0]
                        
                    elif 'G' in list(dat_ef.V[0].split('T')):
                        spd_LEGT = float((dat_ef.V[0].split('T'))[1])
                        aux_LE = dat_ef.V[1]
                        aux_GT = dat_ef.V[0]
                        idx0 = (dat_ef.V[dat_ef.V == aux_LE].index)[0]
                        idx1 = (dat_ef.V[dat_ef.V == aux_GT].index)[0]
                    
                    for ispd, idx in zip(spd_EF.spd, spd_EF.index):
                        if ispd <= spd_LEGT:
                            ef = ((np.float(dat_ef.a[idx0]) * ispd**(np.float(dat_ef.b[idx0]))) + \
                                  (np.float(dat_ef.c[idx0]) * ispd**(np.float(dat_ef.d[idx0]))) + \
                                  (np.float(dat_ef.f[idx0]))) * \
                                  (np.float(dat_ef.k[idx0]))
                            ef_list.append(ef)
                            EF_real_yr.append(iyr)
    
                        elif ispd > spd_LEGT:
                            ef = ((np.float(dat_ef.a[idx1]) * ispd**(np.float(dat_ef.b[idx1]))) + \
                                  (np.float(dat_ef.c[idx1]) * ispd**(np.float(dat_ef.d[idx1]))) + \
                                  (np.float(dat_ef.f[idx1]))) * \
                                  (np.float(dat_ef.k[idx1]))
                            ef_list.append(ef)
                            EF_real_yr.append(iyr)
                elif dat_ef.V.shape[0] > 1 and (dat_ef.V[0] == 0):
                    idx0 = dat_ef.index[0]
                    for ispd, idx in zip(spd_EF.spd, spd_EF.index):
                        ef = ((np.float(dat_ef.a[idx0]) * ispd**(np.float(dat_ef.b[idx0]))) + \
                              (np.float(dat_ef.c[idx0]) * ispd**(np.float(dat_ef.d[idx0]))) + \
                              (np.float(dat_ef.f[idx0]))) * \
                              (np.float(dat_ef.k[idx0]))
                        ef_list.append(ef)
                        EF_real_yr.append(iyr)
                elif dat_ef.shape[0] == 0 :
                    print('*** WARNING *** The emission factor for vehicle {0} of year {1} is zero!'.format(ivhc,iyr))
                    print('*** WARNING *** Check the emissions factor input table')
                    for ispd, idx in zip(spd_EF.spd, spd_EF.index):
                        ef = 0.0
                        ef_list.append(ef)
                        EF_real_yr.append(iyr)
                else: 
                    idx0 = dat_ef.index[0]
                    for ispd, idx in zip(spd_EF.spd, spd_EF.index):
                        ef = ((np.float(dat_ef.a[idx0]) * ispd**(np.float(dat_ef.b[idx0]))) + \
                              (np.float(dat_ef.c[idx0]) * ispd**(np.float(dat_ef.d[idx0]))) + \
                              (np.float(dat_ef.f[idx0]))) * \
                              (np.float(dat_ef.k[idx0]))
                        ef_list.append(ef)
                        EF_real_yr.append(iyr)
    
                spd_EF[ivhc]      = np.asarray(ef_list)
                spd_EF['year']    = EF_real_yr
    
                EmisFact_yr_spd[ivhc][((EmisFact_yr_spd.pollutant == ipoll) & 
                                       (EmisFact_yr_spd.year      == iyr))] = spd_EF[ivhc].values
    EmisFact_yr_spd.to_csv(input_dir+'/intermediate/EmisFactor_yr_spd.csv',sep=';', index=False)            
    run_time = ((time.time() - start_time))
    print("--- %f seconds ---" % (run_time))
    print("--- %f minutes ---" % (run_time/60))
    print("--- %f Hours   ---" % (run_time/3600))
    
    return EmisFact_yr_spd
#EmisFactor_yr_spd = calculate_EmisFact_OLD(input_dir, EF_All_fuel_function)
#EmisFactor_yr_spd = pd.read_csv(inter_dir+'/EmisFactor_yr_spd.csv', sep = ';')



                              
def calculate_EmisFact(Input_dir, Emissions_Factor_Table, ):
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



#teste = EmisFactor_yr_spd.data.copy()
#teste = teste.groupby(['FullName','Pollutant','Years']).mean()
#teste = teste.reset_index(drop=False)
npol = list(EmisFactor_yr_spd.data.Pollutant.unique())

start_time = time.time()

Emissions_by_link = pd.concat([roads_RGS.roads_df, 
                               pd.DataFrame(np.zeros((roads_RGS.roads_df.shape[0],
                                                      AD_SK.fullname.shape[0])),columns=AD_SK.fullname.vhc_name.unique())], axis=1)

EF_by_link = pd.DataFrame(columns=list(Emissions_by_link.columns))

#geocd = [11260105]
#igeocd = 11260105
#ispd = 4.0
#iroad = 101
#ivhc = 'bus_highway_diesel'#'sedan_compact_gasoline'
#ipol = 'CO'  #'CO'

speeds = np.sort(Emissions_by_link.max_speed.unique())
#aux_emis_AD = AD_SK.data.copy()
#aux_emis_AD = aux_emis_AD.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)

Emissions_by_yr_county = pd.DataFrame(columns=list(aux_emis_AD.columns) + ['pollutant'])

for ipol in list(EmisFactor_yr_spd.data.Pollutant.unique()):
    aux_emis_AD = AD_SK.data.copy()
    aux_emis_AD['Pollutant'] = [ipol for i in aux_emis_AD.index]
    Emissions_by_yr_county = Emissions_by_yr_county.append(aux_emis_AD,ignore_index=True, sort=False)
Emissions_by_yr_county = Emissions_by_yr_county.sort_values(by=['region_cd','pollutant','manufacture_date']).reset_index(drop=True)
Emissions_by_yr_county.loc[:,list(AD_SK.fullname.vhc_name)] = 0

aux_emis_AD = AD_SK.data.copy()

for iroad in np.sort(Emissions_by_link.road_type.unique()): #[101,103,104,107]: #
    start_time_loop = time.time()
    aux_emis_EF = EmisFactor_yr_spd.data.copy()
    aux_RD = Emissions_by_link.loc[(Emissions_by_link.road_type == iroad),:].reset_index(drop=True)
    
    
    road_weight = aux_RD.groupby(['region_cd']).sum().reset_index(drop=False)
    
    
    
    aux_avgSpd = avgSpeedDist.data.loc[:,['Speed',str(iroad)]]
    for ispd in avgSpeedDist.data.Speed:
        aux_emis_EF.loc[(aux_emis_EF.spd.isin([ispd])),'emis_fact'] = \
        aux_emis_EF.loc[(aux_emis_EF.spd.isin([ispd])),'emis_fact'].values * \
        aux_avgSpd.loc[(aux_avgSpd.Speed.isin([ispd])),str(iroad)].values
    
    aux_emis_EF = aux_emis_EF.drop(columns=['Temperatures','V','ambient_temp','spd'])
    aux_emis_EF = aux_emis_EF.groupby(['FullName', 'Pollutant','Years']).sum()
    aux_emis_EF = aux_emis_EF.reset_index(drop=False)
    

    for ipol in list(EmisFactor_yr_spd.data.Pollutant.unique()):
        aux_AD_county = aux_emis_AD.loc[(aux_emis_AD.region_cd.isin(list(road_weight.region_cd.unique())))].reset_index(drop=True)

    #    icnty = 11260105
        for icnty in list(road_weight.region_cd.unique()): #geocd: #
            aux_df = aux_AD_county.loc[(aux_AD_county.region_cd.isin([icnty])),:].reset_index(drop=True)
            
            aux_df.loc[(aux_df.region_cd.isin([icnty])),list(AD_SK.fullname.vhc_name)] = \
            aux_AD_county.loc[(aux_AD_county.region_cd.isin([icnty])),list(AD_SK.fullname.vhc_name)].values * \
            road_weight.loc[(road_weight.region_cd.isin([icnty])),'link_split_county'].values
            
            for ivhc in list(AD_SK.fullname.vhc_name):
                aux_df_yr_EF = aux_emis_EF.loc[(aux_emis_EF.FullName.isin([ivhc])) &
                                               (aux_emis_EF.Pollutant.isin([ipol])), ['Years','emis_fact']].reset_index(drop=True)
                aux_df_yr_EF[ivhc] = aux_df_yr_EF.emis_fact.values
                aux_df_yr_EF = aux_df_yr_EF.drop(columns = 'emis_fact').reset_index(drop=True) 
                aux_emisfact = pd.DataFrame({'manufacture_date' : aux_df.manufacture_date})
                aux_emisfact = pd.merge(aux_emisfact,aux_df_yr_EF, left_on = 'manufacture_date', right_on = 'Years', how='left').fillna(method='bfill')
                aux_emisfact = aux_emisfact.fillna(method='ffill')
                aux_df.loc[:,[ivhc]] = aux_df.loc[:,ivhc].values * aux_emisfact[ivhc].values
            
            
            year_list_EF = list(aux_df.manufacture_date)
            Emissions_by_yr_county.loc[(Emissions_by_yr_county.region_cd.isin([icnty])) &
            (Emissions_by_yr_county.manufacture_date.isin(year_list_EF)) &
            (Emissions_by_yr_county.pollutant.isin([ipol])), list(AD_SK.fullname.vhc_name)] += \
                                     aux_df.loc[:,list(AD_SK.fullname.vhc_name)]
        
    run_time_loop = ((time.time() - start_time_loop))
    print("--- %f seconds LOOP roads ---" % (run_time_loop))

run_time = ((time.time() - start_time))
print("--- %f seconds ---" % (run_time))
print("--- %f minutes ---" % (run_time/60))
print("--- %f Hours   ---" % (run_time/3600))
    





    yr_list  = []
    pol_list = []

    for z in aux_emis_EF.Pollutant.unique():
        for i in np.sort(aux_emis_AD.manufacture_date.unique()):
            pol_list.append(z)
            yr_list.append(i)

    new_EF_table = pd.DataFrame({'Years'     : yr_list,
                                 'Pollutant' : pol_list})
    for ivhc in list(AD_SK.fullname.vhc_name):
                aux_emisfact = aux_emis_EF.loc[(aux_emis_EF.FullName.isin([ivhc])) &
                                               (aux_emis_EF.Pollutant.isin([ipol])), ['Years','emis_fact']].reset_index(drop=True)



    
    aux_RD.link_split_county[(aux_RD.max_speed == ispd)].sum()





grouped = EmisFactor_yr_spd.data.copy()
grouped.index = grouped.FullName
grouped = grouped.unstack().fillna(0)
grouped.columns = [x[1] for x in grouped.columns]
grouped.reset_index(inplace=True)







    
    
#    aux_emis_AD.loc[:,AD_SK.fullname.vhc_name] *=  spd_weight
    aux_emisfact = EmisFactor_yr_spd.data.loc[(EmisFactor_yr_spd.data.Pollutant == ipol)]
    aux_emisfact = aux_emisfact.groupby(['FullName','Pollutant','Temperatures','V','Years','ambient_temp']).mean()
    aux_emisfact = aux_emisfact.reset_index(drop=False)
    
    aux_emis_AD.loc[:,AD_SK.fullname.vhc_name] *=  spd_weight
    aux_emisfact = EmisFactor_yr_spd.data.loc[(EmisFactor_yr_spd.data.Pollutant == ipol) &
                                              (EmisFactor_yr_spd.data.spd == ispd)]
    aux_emisfact = aux_emisfact.drop(columns = ['Pollutant', 'Temperatures', 'V', 'spd'])
    inter_emis_yr_county = aux_emis_AD.copy()
    inter_emis_yr_county.loc[:,AD_SK.fullname.vhc_name] = 0
    inter_EF = pd.DataFrame(columns = ['manufacture_date'])
    inter_EF.manufacture_date = aux_emis_AD.manufacture_date.values
    for ivhc in AD_SK.fullname.vhc_name:
        aux2 = aux_emisfact.loc[(aux_emisfact.FullName == ivhc),['Years','emis_fact']]
        aux2[ivhc] = aux2.emis_fact.values
        aux2 = aux2.drop(columns = 'emis_fact').reset_index(drop=True)
        aux2.Years = aux2.Years.astype(int)
        inter_EF = pd.merge(inter_EF,aux2, left_on = 'manufacture_date', right_on = 'Years', how='left').fillna(method='bfill')
        inter_EF = inter_EF.fillna(method='ffill')
        inter_EF = inter_EF.drop(columns = 'Years')
        aux_emis_AD.loc[:,ivhc] = aux_emis_AD[ivhc].values * inter_EF[ivhc].values
        Emissions_by_link.loc[(Emissions_by_link.region_cd == igeocd) &
                              (Emissions_by_link.max_speed == ispd), [ivhc]] = \
                        aux_emis_AD[ivhc].sum() * \
                        (Emissions_by_link.loc[(Emissions_by_link.region_cd == igeocd) &
                                               (Emissions_by_link.max_speed == ispd), ['link_length']] / \
                         Emissions_by_link.loc[(Emissions_by_link.region_cd == igeocd) &
                                               (Emissions_by_link.max_speed == ispd), ['link_length']].sum().values).values
    
    inter_emis_yr_county.loc[:,AD_SK.fullname.vhc_name] += aux_emis_AD.loc[:,AD_SK.fullname.vhc_name]



Emissions_by_yr_county = Emissions_by_yr_county.append(inter_emis_yr_county, ignore_index=True, sort=True)


run_time = ((time.time() - start_time))
print("--- %f seconds ---" % (run_time))
print("--- %f minutes ---" % (run_time/60))
print("--- %f Hours   ---" % (run_time/3600))






















for igeocd in geocd: #
    speeds = np.sort(Emissions_by_link.max_speed.unique())
    for ispd in speeds:
        aux_emis_AD = AD_SK.data[(AD_SK.data.region_cd == igeocd)]
        aux_RD = Emissions_by_link[(Emissions_by_link.region_cd == igeocd)]
        spd_weight = aux_RD.link_split_county[(aux_RD.max_speed == ispd)].sum()
        aux_emis_AD.loc[:,AD_SK.fullname.vhc_name] *=  spd_weight
        aux_emisfact = EmisFactor_yr_spd.data.loc[(EmisFactor_yr_spd.data.Pollutant == ipol) &
                                                  (EmisFactor_yr_spd.data.spd == ispd)]
        aux_emisfact = aux_emisfact.drop(columns = ['Pollutant', 'Temperatures', 'V', 'spd'])
        inter_emis_yr_county = aux_emis_AD.copy()
        inter_emis_yr_county.loc[:,AD_SK.fullname.vhc_name] = 0
        inter_EF = pd.DataFrame(columns = ['manufacture_date'])
        inter_EF.manufacture_date = aux_emis_AD.manufacture_date.values
        for ivhc in AD_SK.fullname.vhc_name:
            aux2 = aux_emisfact.loc[(aux_emisfact.FullName == ivhc),['Years','emis_fact']]
            aux2[ivhc] = aux2.emis_fact.values
            aux2 = aux2.drop(columns = 'emis_fact').reset_index(drop=True)
            aux2.Years = aux2.Years.astype(int)
            inter_EF = pd.merge(inter_EF,aux2, left_on = 'manufacture_date', right_on = 'Years', how='left').fillna(method='bfill')
            inter_EF = inter_EF.fillna(method='ffill')
            inter_EF = inter_EF.drop(columns = 'Years')
            aux_emis_AD.loc[:,ivhc] = aux_emis_AD[ivhc].values * inter_EF[ivhc].values
            Emissions_by_link.loc[(Emissions_by_link.region_cd == igeocd) &
                                  (Emissions_by_link.max_speed == ispd), [ivhc]] = \
                            aux_emis_AD[ivhc].sum() * \
                            (Emissions_by_link.loc[(Emissions_by_link.region_cd == igeocd) &
                                                   (Emissions_by_link.max_speed == ispd), ['link_length']] / \
                             Emissions_by_link.loc[(Emissions_by_link.region_cd == igeocd) &
                                                   (Emissions_by_link.max_speed == ispd), ['link_length']].sum().values).values
        
        inter_emis_yr_county.loc[:,AD_SK.fullname.vhc_name] += aux_emis_AD.loc[:,AD_SK.fullname.vhc_name]



    Emissions_by_yr_county = Emissions_by_yr_county.append(inter_emis_yr_county, ignore_index=True, sort=True)


run_time = ((time.time() - start_time))
print("--- %f seconds ---" % (run_time))
print("--- %f minutes ---" % (run_time/60))
print("--- %f Hours   ---" % (run_time/3600))


ivhc = 'sedan_midsize_lpg' #'sedan_compact_diesel'
ipol = 'CO'
iyr = 2010
dat_EF = EmisFactor_yr_spd.data
aaa = dat_EF[(dat_EF.Pollutant == ipol) &
                             (dat_EF.FullName == ivhc)] &
                             (dat_EF.Years == iyr)]
#bbb = EmisFactor_yr_spd.loc[(EmisFactor_yr_spd.pollutant == ipol) &
#                        (EmisFactor_yr_spd.year == iyr), ['spd',ivhc]]
#ccc = Temp_EF_Diesel[(Temp_EF_Diesel.Pollutant == ipol) &
#                             (Temp_EF_Diesel.FullName == ivhc)]
#
#
#plt.plot(aaa.spd, aaa.emis_fact)
#plt.plot(bbb.spd,bbb[ivhc])
#
#vhc_list = ['sedan_supercompact_diesel',
#             'sedan_compact_diesel',
#             'sedan_midsize_diesel',
#             'sedan_fullsize_diesel']
##             'van_compact_diesel',
##             'van_midsize_diesel',
##             'van_fullsize_diesel',
##             'truck_compact_diesel',
##             'truck_midsize_diesel',
##             'truck_fullsize_diesel',
##             'suv_compact_diesel',
##             'suv_midsize_diesel']
#
#
#for ipol in ['NOx']:
#    for ivhc in vhc_list:
#        print(ivhc)
#        for iyr in dat_EF.Years.unique():
#            aaa = dat_EF[(dat_EF.Pollutant == ipol) &
#                         (dat_EF.FullName == ivhc) &
#                         (dat_EF.Years == iyr)]
#            bbb = EmisFactor_yr_spd.loc[(EmisFactor_yr_spd.pollutant == ipol) &
#                                    (EmisFactor_yr_spd.year == iyr), ['spd',ivhc]]
#            
#            if aaa.shape[0] != 0:
#                fig, ax = plt.subplots(1, figsize = (13,8))
#                plt.plot(aaa.spd, aaa.emis_fact)
#                plt.plot(bbb.spd,bbb[ivhc])
#                ax.set_title('{0}_{1}_{2}'.format(ivhc, ipol, iyr), fontsize = 16)
#        
#                name = '{0}_{1}_{2}.png'.format(ivhc, ipol, iyr)
#                fig.savefig(output_dir+'/'+name, bbox_inches = 'tight',dpi=150)
#                plt.close()
#    