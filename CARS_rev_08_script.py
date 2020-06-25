# -*- coding: utf-8 -*-
"""
@author: Rizzieri Pedruzzi - pedruzzi@email.unc.edu or rizzieri.pedruzzi@gmail.com
"""
import os, sys, time, matplotlib, glob
import pandas as pd
import geopandas as gpd
import numpy as np
import datetime as dt
from shapely.geometry import Polygon
from shapely import wkt
import pyproj, csv
import calendar
import matplotlib.pyplot as plt
from netCDF4 import Dataset
#matplotlib.use('Agg')
total_start_time = time.time()

'''
*****      Case name and path to the folders - Please, set for all directories      *****
'''
case_name = 'EPA_Spd_RR' #'_Seoul_Old_AD'#
home_dir   = r'C:/Users/rizzi/OneDrive - University of North Carolina at Chapel Hill/0000_EUA_IE/001_mobile_source_inventory/CARS_source_code/CARS_test_Country_KNU'
src_dir    = home_dir+'/src'
input_dir  = home_dir+'/input_country'
inter_dir  = home_dir+'/intermediate'#+case_name
output_dir = home_dir+'/output_'+case_name
met_dir    = input_dir+'/metdata'

'''
*****      Dates and Runlenght      *****
'''
STDATE = '2017-01-01'  # start date
ENDATE = '2017-01-02'  # end date
STTIME =  00            # start time 
RUNLEN =  24            # run length


'''
Do you want to calculate the county emissions or it was already calculated?
If you want to calculate them, set calculate_county_emissions = 'yes'
if no, set calculate_county_emissions = 'no'
If it is set to 'no', CARS will look at the output tables of the previous 
emissions calculation.
'''

calculate_county_emissions = 'yes'



'''
*****      INPUT FILES, SHAPEFILES and attibutes      *****
Set in this part the:
Activity data file;
List of Emissions factors tables (eg. in python, the list is set by brackets [])
Average Speed Distribution to calculate the 16 Speed bins emissions factors
Link and county Shape file with their attibute values
**** WARNING ****: the attibute information are from the attribute of the shapefile
make sure you are setting the correct name.
'''
temp_max               = 17.8
temp_mean              = 12.9
temp_min               =  7.8

'''
Processing the Activity data
Please, set the name of the activity data csv file and if you want to process it.
If the user set to not process it (e.g. 'no'), CARS will assume that it was already
been processed and it will look for the intermediate file of acticity data.  
'''
process_activity_data  = 'yes' #'no' 
activity_file          =  'cars_input_v3.2.2.csv' #'cars_input_v3.1.csv' #'seoul_2017_OLD.csv' #



'''
Processing the Emissions Factor
Please, insert the names of the emissions factor csv files in the list Emis_Factor_list.
If the user set to not process it (e.g. 'no'), CARS will assume that it was already
been processed and it will look for the intermediate files of emissions factor.  
'''
process_emissions_factor  = 'yes' #'no' 
Emis_Factor_list          = ['cng_v3.csv','diesel_v3.csv','gasoline_v3.csv','lpg_v3.csv',
                          'H-CNG_v3.csv','H-Diesel_v3.csv','H-Gasoline_v3.csv' ,'H-LPG_v3.csv'] #['gasoline.csv','diesel.csv','cng.csv','lpg.csv']




'''
Processing the Roads and County shape file
Please, insert the names of the Roads and county shape files in 
link_shape and county_shape and the respective attribute table names.
If the user set to not process it (e.g. 'no'), CARS will assume that it was already
been processed and it will look for the intermediate files of emissions factor.  
'''
process_shapes  = 'yes' #'no' 
link_shape             = '/shapes'+'/Road_by_County_SK_AADT_UTM52N.shp'
link_shape_att         = ['LINK_ID'  , 'EMD_CD' , 'EMD_ENG_NM', 'EMD_ENG_NM',
                           'ROAD_RANK', 'spd', 'SHAPE_STLe', 'VKT']
                       
county_shape           = '/shapes/TL_SCCO_EMD.shp'
county_shape_att       = ['EMD_CD', 'EMD_ENG_NM', 'EMD_KOR_NM']



'''
Cold start vehicles and Average Speed profile are mandatory for running CARS
Please, insert the names csv files in the Cold_Start_list and avg_SPD_Dist_file
These two files will be used if   calculate_county_emissions is set to 'yes'
'''
Cold_Start_list        = ['cold_start_vehicles.csv']
avg_SPD_Dist_file      = 'avgSpeedDistribution_rev_01_EPA.csv'  #'avgSpeedDistribution_Average_Country.csv' #  'avgSpeedDistribution_rev_flat.csv' #



'''
Processing Roads restriction for vehicles which are not allow to run in certain roads
Please, insert the names csv files in the road_restriction and set 
process_road_restriction to 'yes'
This step is optional, if you set it to 'no', CARS will assume that all vehicles
run in all roads
These step in processed if calculate_county_emissions is set to 'yes'
'''
process_road_restriction = 'yes' #'no' 
road_restriction       = 'road_restriction_by_vehicles_country.csv'



'''
Do you want to apply Deterioration factor into emissions factors?
If YES, set Deterioration_factor     = 'yes' and insert the deterioration
list file. Remember that vehciles names should match with emissions factor
and Activity data
'''
Deterioration_factor     = 'yes'
Deterioration_list     = ['degradation_rate_Diesel.csv','degradation_rate_Gasoline.csv','degradation_rate_LPG.csv']


'''
Do you want to apply Control factors into emissions?
If YES, set control_strategy     = 'yes' and insert the control
list file. Remember that vehciles names should match with emissions factor
and Activity data
'''
control_strategy     = 'no' #'yes' or 'no'
control_list     = ['control_factors_emergency_sma.csv'] #['control_factors.csv']


'''
'**********      Gridding options      **********'
If the user DO NOT want to generate the grid based on the GRIDDESC,
set the grid_size in meters (e.g 1000) and the grid4AQM = 'no'
'''
grid_size = 27000

'''
CARS has the option to create gridded, time variated and chemical speciated
emissions for air quality modeling.
To run it, set grid4AQM to 'yes'
and point the path for TEMPORAL and CHEMICAL profiles
If you want to output the NetCDF I/O API format, set IOAPI_out to 'yes'
Set the gridfile_name pointing to GRIDDESC file.
Set the Radius_Earth (default=6370000.0)
Set the Datum
'''
grid4AQM = 'no'     #'yes' or 'no'
IOAPI_out = 'no'
gridfile_name = met_dir+'/GRIDDESC_NIER_09_01'
Radius_Earth = 6370000.0
Datum = 'WGS84'

'''
'*****      Temporal profiles names      *****'
'''
temporal_profile_folder = input_dir+'/temporal_profile'
temporal_monthly_file  = 'monthly_profile.csv'
temporal_week_file     = 'week_profile.csv'
temporal_weekday_file  = 'weekday_profile.csv'
temporal_weekend_file  = 'weekend_profile.csv'
temporal_CrossRef      = 'temporal_profile_CrossRef_country.csv' #'temporal_profile_CrossRef_country.csv'

'''
'*****      Chemical speciation file and cross reference      *****'
'''
chemical_profile_folder = input_dir+'/chemical_speciation'
chemical_profile  = 'gspro_cmaq_cb6_2014fa_nata_cb6cmaq_14j_criteria.txt'
speciation_CrossRef   = 'chem_profile_CrossRef.csv'


'''
Do you want to CARS outpout the GRID shapefile?
If yes, set the outGridShape = 'yes'
This option was added because sometimes generate the shapefile
can take couple minutes.
'''
outGridShape = 'no'



''' IT IS NOT IMPLEMENTED
Do you want to estimate emissions for the future ?
If YES, set future_case     = 'yes' and make sure there are activity data and
emissions factor for the vehicles you want to estimate.
Remember that vehciles names should match with emissions factor and Activity data
'''
future_case      = 'no' #'yes' or 'no'


'''
***** PLOTTING *****
Do you want to generate the figures?
One option is the Adjustment of the plot scale. This was added to allow the user
change the scale of all plots, because sometimes, there are one point with high emissions
than the rest of the domain, so there is the need to renormalize the scale.
The  adj_scale will multiply the maximum value of the plot area, so it will reduce
or increase the maximum value of the scale.
To use it, set the adj_scale to the value you want (e.g. 0.5) The default value is 0.4

*** WARNING ***: The plotting takes a bit longer to finish especially the 24 hours plot
'''
plot_figures = 'no'   #'yes' or 'no'
plot_24      = 'no'        #'yes' or 'no' #24 hours animated plot
adj_scale    = 0.4 
cmap         = 'jet'





'''
----------------------------------------------------------------------------------------
End of users input. Please, DO NOT CHANGE the code below this point.
Any bugs or issues report at 
https://github.com/rpedruzzi/CARS

Calling the main code of CARS...
'''
if __name__== "__main__":
    src_dir    = home_dir+'/src'
    sys.path.append(src_dir)
    from CARS_v1_08 import *
    def main(home_dir, input_dir, output_dir, ):
        print('')
        print('******************************************************************')
        print('** Comprehensive Automobile Emissions Research Simulator - CARS **')
        print('******************************************************************')
        print('')
        print('*****    Loading classes, functions and primary input data   *****')
        print('*****                      Please wait                       *****')
        print('')
        
        if not os.path.exists(home_dir+'/intermediate'):
            print('*** intermediate directory does not exist: Creating a new one')
            os.makedirs(home_dir+'/intermediate/')
        
        if not os.path.exists(output_dir):
            print('*** Outdir directory does not exist: Creating a new one')
            os.makedirs(output_dir)
    
        if not os.path.exists(input_dir+'/emissions_factor/'):
            print('emissions_factor directory does not exist: Creating a new one')
            print('Remember to put the emissions factor files, average speed distribuition file,')
            print('cold start vehicle files, degradation factor file and control factors file')
            print('inside the emissions_factor directory')
            os.makedirs(input_dir+'/emissions_factor/')
    
        if not os.path.exists(output_dir+'/plots'):
            print('*** intermediate directory does not exist: Creating a new one')
            os.makedirs(output_dir+'/plots/')
            
        if not os.path.exists(output_dir+'/LOGS'):
                print('*** LOG directory does not exist: Creating a new one')
                os.makedirs(output_dir+'/LOGS')
        
        
        GrayJet = createGrayJet()
        run_period = run_period_times(STDATE, ENDATE, STTIME, RUNLEN)
        
        #Reading chemical profile csv file. Calling read_chemical_info function
        avgSpeedDist = read_avgSpeedDistribution(input_dir, avg_SPD_Dist_file)
    
    
        #Reading Activity data e csv file. Calling read_activity_data_csv_SK function
        if (process_activity_data.lower() == 'yes') or (process_activity_data.lower() == 'y'):
            AD_SK = read_activity_data_csv_SK(input_dir, inter_dir, case_name, activity_file, ENDATE)
            
        elif (process_activity_data.lower() == 'no') or (process_activity_data.lower() == 'n'):
            print('')
            print('******************************************************************')
            print('*****          process_activity_data is set to NO            *****')
            print('*****             reading the intemediate files              *****')
            print('******************************************************************')
            print('')
    
            if (os.path.exists(inter_dir+'/outdf_AD_{0}.csv'.format(case_name)))  and \
                (os.path.exists(inter_dir+'/vhc_names_AD_{0}.csv'.format(case_name))) and \
                (os.path.exists(inter_dir+'/vhc_years_AD_{0}.csv'.format(case_name))) and \
                (os.path.exists(inter_dir+'/vhc_fuels_AD_{0}.csv'.format(case_name))) and \
                (os.path.exists(inter_dir+'/vhc_cars_AD_{0}.csv'.format(case_name)))  and \
                (os.path.exists(inter_dir+'/vhc_numbers_AD_{0}.csv'.format(case_name))):
    
                print('*****          Reading file. Please wait...                 *****')
                outdf     = pd.read_csv(inter_dir+'/outdf_AD_{0}.csv'.format(case_name), sep = None, engine='python')
                vhc_names = pd.read_csv(inter_dir+'/vhc_names_AD_{0}.csv'.format(case_name))
                vhc_years = pd.read_csv(inter_dir+'/vhc_years_AD_{0}.csv'.format(case_name))
                vhc_fuels = pd.read_csv(inter_dir+'/vhc_fuels_AD_{0}.csv'.format(case_name))
                vhc_cars  = pd.read_csv(inter_dir+'/vhc_cars_AD_{0}.csv'.format(case_name))
                count_outdf = pd.read_csv(inter_dir+'/vhc_numbers_AD_{0}.csv'.format(case_name), sep = None, engine='python')       
                AD_SK = Activity_Data_table(outdf, vhc_names, vhc_years, vhc_fuels, vhc_cars, count_outdf)
    
            else:    
                print('')
                print('******************************************************************')
                print('*** ERROR ABORT ***: missing intermediate files of activity data *****')
                print('***** Please, set process_activity_data to yes and run it    *****')
                print('******************************************************************')
                print('')
                sys.exit()
        else:
            AD_SK = read_activity_data_csv_SK(input_dir, inter_dir, case_name, activity_file, ENDATE)
        
        #Reading and calculating Emissions factor. Calling read_emissions_factor_SK and calculate_EmisFact function
        if (process_emissions_factor.lower() == 'yes') or (process_emissions_factor.lower() == 'y'):
            EF_All_fuel_function = read_emissions_factor_SK(input_dir, Emis_Factor_list, inter_dir, case_name, output_dir)
            EmisFactor_yr_spd    = calculate_EmisFact(input_dir, EF_All_fuel_function, temp_min, temp_mean, temp_max,
                                                      case_name, inter_dir, output_dir)
            
            if (Deterioration_factor.lower() == 'yes') or (Deterioration_factor.lower() == 'y'):
                EmisFactor_yr_spd = apply_deterioration_emissions_factor_SK(input_dir, inter_dir, output_dir, 
                                                                            case_name, STDATE, ENDATE, 
                                                                        Deterioration_list, EmisFactor_yr_spd)
        elif (process_emissions_factor.lower() == 'no') or (process_emissions_factor.lower() == 'n'):
            print('')
            print('******************************************************************')
            print('*****        process_emissions_factor is set to NO           *****')
            print('*****             reading the intemediate files              *****')
            print('******************************************************************')
            print('')
            
            if  (os.path.exists(inter_dir+'/'+'dat_EF_{0}.csv'.format(case_name)  )) and \
                (os.path.exists(inter_dir+'/'+'EF_years_{0}.csv'.format(case_name))) and \
                (os.path.exists(inter_dir+'/'+'EF_names_{0}.csv'.format(case_name))) and \
                (os.path.exists(inter_dir+'/'+'EF_fuels_{0}.csv'.format(case_name))) and \
                (os.path.exists(inter_dir+'/'+'EF_polls_{0}.csv'.format(case_name))):
    
                print('*****          Reading files. Please wait...                 *****')
                dat_EF   = pd.read_csv(inter_dir+'/'+'dat_EF_{0}.csv'.format(case_name)  , sep = None, engine='python')
                EF_years = pd.read_csv(inter_dir+'/'+'EF_years_{0}.csv'.format(case_name))
                EF_names = pd.read_csv(inter_dir+'/'+'EF_names_{0}.csv'.format(case_name))
                EF_fuels = pd.read_csv(inter_dir+'/'+'EF_fuels_{0}.csv'.format(case_name))
                EF_polls = pd.read_csv(inter_dir+'/'+'EF_polls_{0}.csv'.format(case_name))
          
                EmisFactor_yr_spd = EF_Grid_table(dat_EF, EF_years, EF_names, EF_fuels, EF_polls)
    
            else:    
                print('')
                print('*************************************************************************')
                print('*** ERROR ABORT ***: missing intermediate files of emissions factor *****')
                print('*****    Please, set process_emissions_factor to yes and run it     *****')
                print('*************************************************************************')
                print('')
                sys.exit()
        else:
            EF_All_fuel_function = read_emissions_factor_SK(input_dir, Emis_Factor_list)
            EmisFactor_yr_spd    = calculate_EmisFact(input_dir, EF_All_fuel_function, 
                                                      temp_min, temp_mean, temp_max,
                                                      case_name, inter_dir, output_dir)
            
            if (Deterioration_factor.lower() == 'yes') or (Deterioration_factor.lower() == 'y'):
                EmisFactor_yr_spd = apply_deterioration_emissions_factor_SK(input_dir, inter_dir, output_dir, 
                                                                            case_name, STDATE, ENDATE, 
                                                                            Deterioration_list, EmisFactor_yr_spd)
    
        #Comparing vehicles between Activity data Emissions factor. Calling check_VHC_AD_EF
        check_VHC_AD_EF(EmisFactor_yr_spd, AD_SK, output_dir)
    
    
        #Reading and processinf Roads and county shape file. Calling roads_RGS and county_SHP
        if (process_shapes.lower() == 'yes') or (process_shapes.lower() == 'y'):
            roads_RGS = roads_grid_surrogate_inf(input_dir, inter_dir, output_dir, link_shape, link_shape_att, case_name, grid_size, grid4AQM,
                                         gridfile_name, Radius_Earth, Unit_meters = True ) 
            county_SHP = processing_County_shape(input_dir, inter_dir, output_dir, county_shape, 
                                                 county_shape_att, case_name, gridfile_name, Radius_Earth,grid4AQM) 
        elif (process_shapes.lower() == 'no') or (process_shapes.lower() == 'n'):
            print('')
            print('******************************************************************')
            print('*****             process_shapes is set to NO                *****')
            print('*****             reading the intemediate files              *****')
            print('******************************************************************')
            print('')
            
    
            if  (os.path.exists(inter_dir+'/'+'grid_{0}.csv'.format(case_name)  )) and \
                (os.path.exists(inter_dir+'/'+'surrogate_{0}.csv'.format(case_name))) and \
                (os.path.exists(inter_dir+'/'+'out_roads_{0}.csv'.format(case_name))) and \
                (os.path.exists(inter_dir+'/'+'out_district_{0}.csv'.format(case_name))):
                
                print('*****          Reading files. Please wait...                 *****')
                with open(inter_dir+'/'+'out_roads_CRS_{0}.txt'.format(case_name), 'r') as crsfile:
                    Lines = crsfile.readlines() 
                crs = Lines[0]
                import csv
                maxInt = sys.maxsize
                while True:
                    # decrease the maxInt value by factor 10 
                    # as long as the OverflowError occurs.
                    try:
                        csv.field_size_limit(maxInt)
                        break
                    except OverflowError:
                        maxInt = int(maxInt/10)
                grid      = pd.read_csv(inter_dir+'/'+'grid_{0}.csv'.format(case_name)      , sep = None, engine='python')
                grid['geometry'] = grid['geometry'].apply(wkt.loads)
                grid      = gpd.GeoDataFrame(grid, crs=crs, geometry='geometry')
                surrogate = pd.read_csv(inter_dir+'/'+'surrogate_{0}.csv'.format(case_name) , sep = None, engine='python')
                surrogate['geometry'] = surrogate['geometry'].apply(wkt.loads)
                surrogate = gpd.GeoDataFrame(surrogate, crs=crs, geometry='geometry')
                out_roads = pd.read_csv(inter_dir+'/'+'out_roads_{0}.csv'.format(case_name) , sep = None, engine='python')
                out_roads['geometry'] = out_roads['geometry'].apply(wkt.loads)
                out_roads = gpd.GeoDataFrame(out_roads, crs=crs, geometry='geometry')
                # county_SHP = pd.read_csv(inter_dir+'/'+'out_district_{0}.csv'.format(case_name), sep = None, engine='python')
                county_SHP = processing_County_shape(input_dir, inter_dir, output_dir, county_shape, county_shape_att, 
                                                     case_name, gridfile_name, Radius_Earth,) 
                roads_RGS  = Roads_Grid_table(grid, surrogate, out_roads, crs)
            else:    
                print('')
                print('*************************************************************************')
                print('*** ERROR ABORT ***: missing intermediate files of Roads and county *****')
                print('*****         Please, set process_shapes to yes and run it          *****')
                print('*************************************************************************')
                print('')
                sys.exit()
        else:
            roads_RGS = roads_grid_surrogate_inf(input_dir, inter_dir, output_dir, link_shape, link_shape_att, case_name, grid_size, grid4AQM,
                                 gridfile_name, Radius_Earth, Unit_meters = True ) 
            county_SHP = processing_County_shape(input_dir, inter_dir, output_dir, county_shape, 
                                                 county_shape_att, case_name, gridfile_name, Radius_Earth,) 
    
        # Calculating the county emissions
        if (calculate_county_emissions.lower() == 'yes') or (calculate_county_emissions.lower() == 'y'):
            County_Emissions = calc_County_Emissions(input_dir, inter_dir, output_dir, case_name, STDATE, ENDATE, 
                                                     temp_min, temp_mean, temp_max,
                                                     EmisFactor_yr_spd, AD_SK, 
                                                     roads_RGS, county_SHP, Cold_Start_list, 
                                                     avgSpeedDist, road_restriction, process_road_restriction)
        elif (calculate_county_emissions.lower() == 'no') or (calculate_county_emissions.lower() == 'n'):
            print('')
            print('******************************************************************')
            print('*****     calculate_county_emissions is set to NO            *****')
            print('*****             reading the intemediate files              *****')
            print('******************************************************************')
            print('')
    
            if  (os.path.exists(output_dir+'/'+'District_Total_Normalized_Emissions_by_Year_{0}.csv'.format(case_name)))  and \
                (os.path.exists(output_dir+'/'+'District_Total_Emissions_Tons_per_Year_{0}.csv'.format(case_name))) and \
                (os.path.exists(output_dir+'/'+'District_Total_Normalized_Emissions_Tons_per_Year_{0}.csv'.format(case_name))) and \
                (os.path.exists(inter_dir+'/'+'district_table_{0}.csv'.format(case_name))) and \
                (os.path.exists(inter_dir+'/'+'years_table_{0}.csv'.format(case_name)))  and \
                (os.path.exists(inter_dir+'/'+'vhc_table_{0}.csv'.format(case_name)))  and \
                (os.path.exists(inter_dir+'/'+'road_by_district_{0}'.format(case_name))):
    
                print('*****          Reading file. Please wait...                 *****')
                Emissions_by_yr_county_WGT = pd.read_csv(output_dir+'/'\
                                                         +'District_Total_Normalized_Emissions_by_Year_{0}.csv'.format(case_name), sep=',')
                total_county     = pd.read_csv(output_dir+'/'+'District_Total_Emissions_Tons_per_Year_{0}.csv'.format(case_name), sep=',')
                total_county_WGT = pd.read_csv(output_dir+'/'+'District_Total_Normalized_Emissions_Tons_per_Year_{0}.csv'.format(case_name), sep=',')
                county_table     = pd.read_csv(inter_dir+'/'+'district_table_{0}.csv'.format(case_name), sep = ',')
                years_table      = pd.read_csv(inter_dir+'/'+'years_table_{0}.csv'.format(case_name), sep = ',')
                vhc_table        = pd.read_csv(inter_dir+'/'+'vhc_table_{0}.csv'.format(case_name), sep = ',')
                road_by_county   = pd.read_csv(inter_dir+'/'+'road_by_district_{0}'.format(case_name), sep = ',')
                County_Emissions = Emissions_table(Emissions_by_yr_county_WGT, total_county, 
                                                   total_county_WGT, county_table, years_table, 
                                                   vhc_table, road_by_county)
            else:    
                print('')
                print('*************************************************************************')
                print('*** ERROR ABORT ***: missing intermediate files of county emissions *****')
                print('*****   Please, set calculate_county_emissions to yes and run it    *****')
                print('*************************************************************************')
                print('')
                sys.exit()
        else:
            County_Emissions = calc_County_Emissions(input_dir, inter_dir, output_dir, case_name, STDATE, ENDATE, 
                                                     temp_min, temp_mean, temp_max,
                                                     EmisFactor_yr_spd, AD_SK, 
                                                     roads_RGS, county_SHP, Cold_Start_list, 
                                                     avgSpeedDist, road_restriction, process_road_restriction)
    
        # Applying control strategy
        if (control_strategy.lower() == 'yes') or (control_strategy.lower() == 'y'):
            aplly_control(input_dir, inter_dir, output_dir, case_name, STDATE, ENDATE, 
                          control_list, County_Emissions, county_SHP)
        else:
            print('')
            print('******************************************************************')
            print('*****           control_strategy is set to NO                *****')
            print('******************************************************************')
            print('')
    
        # Processing chemical speciation, temporal allocation and GRID info
        if (grid4AQM.lower() == 'yes') or (grid4AQM.lower() == 'y'):
            print('')
            print('******************************************************************')
            print('*****                grid4AQM is set to YES                  *****')
            print('*****  reading TEMPORAL and CHEMICAL profiles and GRID info  *****')
            print('******************************************************************')
            print('')
    
            # Reading grid info. Calling read_griddesc function
            GRID_info = read_griddesc(input_dir, gridfile_name, grid4AQM = grid4AQM)
            
            #Reading temporal profile csv file. Calling read_temporal_info function
            TempPro = read_temporal_info(temporal_profile_folder, temporal_monthly_file, temporal_week_file,
                                   temporal_weekday_file, temporal_weekend_file, temporal_CrossRef, run_period)
        
            #Reading chemical profile csv file. Calling read_chemical_info function
            Chemical_Spec_Table = read_chemical_info(chemical_profile_folder, chemical_profile, speciation_CrossRef)
    
            County_Emissions_ChemSpec = chemical_speciation(input_dir, County_Emissions,
                                                    AD_SK, Chemical_Spec_Table)
    
            gridded_emissions = gridded_emis_NC(input_dir, inter_dir, output_dir, case_name, STDATE, ENDATE,
                                                AD_SK, TempPro, County_Emissions_ChemSpec, roads_RGS, 
                                                run_period, GRID_info, IOAPI_out = 'Y', grid4AQM = 'Y')
            
        else:
            print('')
            print('******************************************************************')
            print('*****               grid4AQM is set to NO                    *****')
            print('******************************************************************')
            print('')
            # Reading grid info. Calling read_griddesc function
            GRID_info = read_griddesc(input_dir, gridfile_name, grid4AQM = grid4AQM)
            
            #Reading temporal profile csv file. Calling read_temporal_info function
            TempPro = read_temporal_info(temporal_profile_folder, temporal_monthly_file, temporal_week_file,
                                   temporal_weekday_file, temporal_weekend_file, temporal_CrossRef, run_period)
        
            #Reading chemical profile csv file. Calling read_chemical_info function
            Chemical_Spec_Table = read_chemical_info(chemical_profile_folder, chemical_profile, speciation_CrossRef)
    
            County_Emissions_ChemSpec = chemical_speciation(input_dir, County_Emissions,
                                                    AD_SK, Chemical_Spec_Table)
    
            gridded_emissions = gridded_emis_NC(input_dir, inter_dir, output_dir, case_name, STDATE, ENDATE,
                                                AD_SK, TempPro, County_Emissions_ChemSpec, roads_RGS, 
                                                run_period, GRID_info, IOAPI_out = 'N', grid4AQM = 'N')
    
        if (plot_figures.lower() == 'yes') or (plot_figures.lower() == 'y'):
            plot_figures_chart(input_dir, inter_dir, output_dir, case_name, STDATE, ENDATE,
                               County_Emissions, County_Emissions_ChemSpec, AD_SK, gridded_emissions, 
                               roads_RGS, avgSpeedDist, run_period,
                               plot_figures = plot_figures, plot_24 = plot_24, 
                               pol_list = [], species_list = ['NO2', 'CO'])
        else:
            print('')
            print('******************************************************************')
            print('*****            plot_figures is set to NO                   *****')
            print('******************************************************************')
            print('')
    
    
        total_run = ((time.time() - total_start_time))
        print('')  
        print('---  CARS Total Elapsed time in seconds = {0}     ---'.format(total_run))
        print('---  CARS Total Elapsed time in minutes = {0}     ---'.format(total_run/60))
        print('---  CARS Total Elapsed time in   hours = {0}     ---'.format(total_run/3600))    
        print('')
        print('******************************************************************')
        print('*****                CARS finished to run                    *****')
        print('******************************************************************')
        print('')  
    main(home_dir, input_dir, output_dir, )

