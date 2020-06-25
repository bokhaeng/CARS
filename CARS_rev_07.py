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
home_dir   = r'C:/Users/pedruzzi/OneDrive - University of North Carolina at Chapel Hill/0000_EUA_IE/001_mobile_source_inventory/CARS_source_code/CARS_test_Country_KNU'
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




#    now = dt.datetime.now()
'''
----------------------------------------------------------------------------------------
End of users input. Please, DO NOT CHANGE the code below this point.
Any bugs or issues report at 
https://github.com/rpedruzzi/CARS
'''

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



class EmissionFactor_table:
    def __init__(self, dataframe, name ):
        self.dataframe = dataframe
        self.name      = name.split('.')[0]

class Activity_Data_table:
    def __init__(self, dataframe, vhc_name, years, fuels, vehicle, vhc_count):
        self.data      = dataframe
        self.fullname  = vhc_name
        self.years     = years
        self.fuels     = fuels
        self.vhc       = vehicle
        self.vhc_count = vhc_count

class Roads_Grid_table:
    def __init__(self, grid_dataframe, surrogate, roads_df, crs):
        self.grid      = grid_dataframe
        self.surrogate = surrogate
        self.roads_df  = roads_df
        self.crs       = crs

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
    def __init__(self, County_Emissions, County_Emissions_GeoRef, 
                 County_Emissions_GeoRef_WGT,County, Years, vhc_name, road_by_county):
        self.county_by_yr     = County_Emissions
        self.county_total     = County_Emissions_GeoRef
        self.county_total_WGT = County_Emissions_GeoRef_WGT
        self.counties         = County
        self.years            = Years
        self.fullname         = vhc_name
        self.road_by_county   = road_by_county

class Emissions_ChemSpec:
    def __init__(self, ChemSpec_emissions, Grams_pollutants, Moles_pollutants):
        self.chemspec_emissions  = ChemSpec_emissions
        self.grams_pol           = Grams_pollutants
        self.moles_pol           = Moles_pollutants



class Temporal_table:
    def __init__(self, temporal_monthly, temporal_week, temporal_weekday, temporal_weekend, temporal_CrossRef, Diurnal_profile):
        self.monthly  = temporal_monthly
        self.week     = temporal_week
        self.weekday  = temporal_weekday
        self.weekend  = temporal_weekend
        self.crossref = temporal_CrossRef
        self.diurnalPro = Diurnal_profile

class Chemical_Speciation_table:
    def __init__(self, chemical_profile, speciation_CrossRef):
        self.chempro  = chemical_profile
        self.crossref = speciation_CrossRef

class GRID_info_table:
    def __init__(self, NTHIK, NCOLS, NROWS, NLAYS, GDTYP, P_ALP, P_BET, P_GAM,
                 XCENT, YCENT, XORIG, YORIG, XCELL, YCELL, GDNAM):
        self.NTHIK = NTHIK
        self.NCOLS = NCOLS
        self.NROWS = NROWS
        self.NLAYS = NLAYS
        self.GDTYP = GDTYP
        self.P_ALP = P_ALP
        self.P_BET = P_BET
        self.P_GAM = P_GAM
        self.XCENT = XCENT
        self.YCENT = YCENT
        self.XORIG = XORIG
        self.YORIG = YORIG
        self.XCELL = XCELL
        self.YCELL = YCELL
        self.GDNAM = GDNAM
    


# =============================================================================
# Function to create the GrayJet colormap
# =============================================================================
def createGrayJet():
    from matplotlib import cm
    from matplotlib.colors import ListedColormap
    import numpy as np
    jet = cm.get_cmap('jet', 256)
    newcolors = jet(np.linspace(0, 1, 256))
    newcolors[0, :] = np.array([256/256, 256/256, 256/256, 1])    #white
    newcolors[1:4, :] = np.array([200/256, 200/256, 200/256, 1])  #lightgray
    GrayJet = ListedColormap(newcolors)
    return GrayJet
# =============================================================================


# =============================================================================
# Function to define the run period
# =============================================================================
def run_period_times(STDATE, ENDATE, STTIME, RUNLEN):
    STDATE = STDATE+' {:>02.2s}'.format(str(STTIME))
    RUNLEN = RUNLEN+1
    df_times = pd.date_range(start=STDATE, freq='H', periods=RUNLEN )
    run_period = pd.DataFrame({'DateTime' : df_times})
    run_period['dayofweek'] = run_period.DateTime.dt.strftime('%a').str.lower()  #%a or %A
    run_period['month'] = run_period.DateTime.dt.strftime('%b').str.lower()  #%a or %A
    run_period[['Day','Hour']] = run_period['DateTime'].dt.strftime('%Y-%m-%d_%H:%M:%S').str.split('_',expand=True)
    run_period['jul_day']  = run_period.DateTime.dt.strftime('%Y%j').astype(int)
    run_period['jul_hour'] = run_period.DateTime.dt.strftime('%H0000').astype(int)
    run_period['TFLAG']    = list(zip(run_period.jul_day, run_period.jul_hour))

    return run_period

        
# =============================================================================
# Function to read Temporal information
# =============================================================================
def read_griddesc(input_dir, gridfile_name, grid4AQM = 'no'):
    start_time = time.time()
    if ((grid4AQM.lower() == 'yes') or (grid4AQM.lower() == 'y')) and (gridfile_name != ''):


        input_dir = input_dir
        gridfile_name = gridfile_name
        print('******************************************************************')
        print('*****                Processing GRIDDESC                     *****')
        print('*****                    Please wait ...                     *****')
        print('******************************************************************')
        print('')
    
    
        if os.path.exists(gridfile_name) == True:
            print ('')
            print ('Reading GridDesc file to generate gridded emissions ...')
            print (gridfile_name)
            griddesc = pd.read_csv(gridfile_name, header=None,engine='python')
            
            GDTYP, P_ALP, P_BET, P_GAM, XCENT, YCENT = griddesc.loc[2,0].split()
            P_ALP, P_BET, P_GAM, XCENT, YCENT = map(float, [P_ALP, P_BET, P_GAM, XCENT, YCENT])
            GDTYP = int(GDTYP)
            GDNAM = griddesc.loc[4,0].split()[0].replace("'",'')
            
            COORD_NAME, XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK = griddesc.loc[5,0].split()
            XORIG, YORIG, XCELL, YCELL = map(float, [XORIG, YORIG, XCELL, YCELL])
            NCOLS, NROWS, NTHIK = map(int, [NCOLS, NROWS, NTHIK])
            NLAYS = 1 #for now NLAYS is equal 1 because CARS generates emissions only for the first level
        
        else:
            print ('')
            print('*** ERROR ABORT ***:  Griddesc file {0} does not exist!'.format(gridfile_name))
            sys.exit('CARS can not read Griddesc file')
 

        run_time = ((time.time() - start_time))
        print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
        print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
        print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
        print('')
        print('******************************************************************')
        print('*****             Processing GRIDDESC is done               *****')
        print('******************************************************************')
        print('')    
       
        return GRID_info_table(NTHIK, NCOLS, NROWS, NLAYS, GDTYP, P_ALP, P_BET, P_GAM, 
                               XCENT, YCENT, XORIG, YORIG, XCELL, YCELL, GDNAM)

    else:
        run_time = ((time.time() - start_time))
        print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
        print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
        print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
        print('')
        print('******************************************************************')
        print('*****             Processing GRIDDESC is done               *****')
        print('******************************************************************')
        print('')    

        return 





# =============================================================================
# Function to read Temporal information
# =============================================================================
def read_temporal_info(Temporal_Profile_Folder, Temporal_Monthly, Temporal_Week,
                       Temporal_WeekDay, Temporal_WeekEnd, Temporal_CrossRef, Run_Period):
    start_time = time.time()
    print('******************************************************************')
    print('*****           Processing temporal profile                  *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    temp_dir = Temporal_Profile_Folder   #temporal_profile_folder #
#    sep = ';'
#    Temporal_Monthly   = 'monthly_profile.csv'
#    Temporal_Week      = 'week_profile.csv'
#    Temporal_WeekDay   = 'weekday_profile.csv'
#    Temporal_WeekEnd   = 'weekend_profile.csv'
#    Temporal_CrossRef  = 'temporal_profile_CrossRef_country.csv' #'temporal_profile_CrossRef.csv'
    month    = '{0}{1}'.format(temp_dir+'/',Temporal_Monthly)
    week     = '{0}{1}'.format(temp_dir+'/',Temporal_Week)
    weekday  = '{0}{1}'.format(temp_dir+'/',Temporal_WeekDay)
    weekend  = '{0}{1}'.format(temp_dir+'/',Temporal_WeekEnd)
    crossref = '{0}{1}'.format(temp_dir+'/',Temporal_CrossRef)
    Times = Run_Period  #run_period #
    files = [month, week, weekday, weekend, crossref]
    for ifl in files:
        if os.path.exists(ifl) == False:
            print ('')
            print('*** ERROR ABORT ***:  Temporal file ', ifl, ' "" does not exist!')
            sys.exit('CARS preProcessor can not read Temporal profile')
        else:
            print ('')
            print ('*** Reading Temporal information : ***')
            print (ifl)

    with open(month, 'r') as csvfile:
        sep = csv.Sniffer().sniff(csvfile.read(40960)).delimiter

    month_out    = pd.read_csv(month, sep = None, engine = 'python').fillna(np.nan)
    week_out     = pd.read_csv(week, sep = None, engine = 'python').fillna(np.nan)
    weekday_out  = pd.read_csv(weekday, sep = None, engine = 'python').fillna(np.nan)
    weekend_out  = pd.read_csv(weekend, sep = None, engine = 'python').fillna(np.nan)
    crossref_out = pd.read_csv(crossref, sep = None, engine = 'python').fillna(np.nan)

    month_out.columns    = month_out.columns.str.lower()
    week_out.columns     = week_out.columns.str.lower()
    weekday_out.columns  = weekday_out.columns.str.lower()
    weekend_out.columns  = weekend_out.columns.str.lower()
    crossref_out.columns = crossref_out.columns.str.lower()
    crossref_out['fullname'] = crossref_out.vehicle.str.cat(crossref_out[['engine','fuel']], sep='_')

    CR_month = crossref_out.loc[:,['fullname','road_type','month']]
    CR_week  = crossref_out.loc[:,['fullname','road_type','week']]
    CR_wday  = crossref_out.loc[:,['fullname','road_type','weekday']]
    CR_wend  = crossref_out.loc[:,['fullname','road_type','weekend']]
    
    CR_month = pd.merge(CR_month, month_out  , left_on='month'   , right_on='profile', how='left')
    CR_week  = pd.merge(CR_week , week_out   , left_on='week'    , right_on='profile', how='left')
    CR_wday  = pd.merge(CR_wday , weekday_out, left_on='weekday' , right_on='profile', how='left')
    CR_wend  = pd.merge(CR_wend , weekend_out, left_on='weekend' , right_on='profile', how='left')
    
    diurnal_temp = pd.DataFrame(columns=['Day']+list(CR_wday.columns)) #['DateTime']+list(CR_wday.fullname)) #
    for i in Times.Day.unique():
        x,y,z = i.split('-')
        imonth = (dt.date(int(x), int(y), int(z))).strftime('%b').lower()
        DofW   = (dt.date(int(x), int(y), int(z))).strftime('%a').lower()
        mon_week_F = pd.merge(CR_month.loc[:,['fullname','road_type',imonth]] ,
                               CR_week.loc[:,['fullname','road_type',DofW]] ,
                               left_on=['fullname','road_type'] , 
                               right_on=['fullname','road_type'], how='left')
        ndays = pd.Period(i).days_in_month #getting the number of days in a month to apply into weekl profile
        mon_week_F['factor'] = mon_week_F[imonth] * ( mon_week_F[DofW] * (7/ndays))
        colsD = [str(y) for y in range(0,24)]
        if (DofW == 'sun') or (DofW == 0) or (DofW == 'sat') or (DofW == 6):
            aux_diurnal = CR_wend.copy()
            aux_diurnal['Day'] = [i for x in aux_diurnal.index]
            aux_diurnal.loc[:,colsD] = aux_diurnal.loc[:,colsD].apply(lambda x: np.asarray(x) * mon_week_F['factor'].values)
        else:
            aux_diurnal = CR_wday.copy()
            aux_diurnal['Day'] = [i for x in aux_diurnal.index]
            aux_diurnal.loc[:,colsD] = aux_diurnal.loc[:,colsD].apply(lambda x: np.asarray(x) * mon_week_F['factor'].values)
        
        diurnal_temp = diurnal_temp.append(aux_diurnal,ignore_index=True, sort=False)
    
    diurnal_temp = diurnal_temp.loc[:,['Day', 'fullname', 'road_type', 'profile']+list((np.arange(0,24)).astype(str))]
        
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    print('')
    print('******************************************************************')
    print('*****             Temporal calculation is done               *****')
    print('******************************************************************')
    print('')    
    return Temporal_table(month_out, week_out, weekday_out, weekend_out, crossref_out, diurnal_temp)





# =============================================================================
# Function to read Chemical speciation
# =============================================================================
def read_chemical_info(Chemical_Speciation_Folder, chemical_profile, Speciation_CrossRef):
    start_time = time.time()
    print('******************************************************************')
    print('*****             Reading chemical profile                   *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    
    
    spec_dir = Chemical_Speciation_Folder  #chemical_profile_folder#
    chempro  = '{0}{1}'.format(spec_dir+'/',chemical_profile)
    crossref = '{0}{1}'.format(spec_dir+'/',Speciation_CrossRef) #speciation_CrossRef)#
    

    files = [chempro, crossref]
    for ifl in files:
        if os.path.exists(ifl) == False:
            print ('')
            print('*** ERROR ABORT ***:  Chemical Speciation file ', ifl, ' "" does not exist!')
            sys.exit('CARS preProcessor can not read Chemical Speciation file')
        else:
            print ('')
            print ('*** Reading Chemical Speciation information ... ***')
            print (ifl)

    with open(chempro, 'r') as csvfile:
        sep = csv.Sniffer().sniff(csvfile.read(100000)).delimiter  
        
    chempro_out  = pd.read_csv(chempro, sep = sep, comment='#', 
                               header=None, usecols=[0,1,2,3,4],
                               names=['profile', 'pollutant', 'species', 'fraction', 'mw'],
                               engine ='python').fillna(0)
    chempro_out['pollutant'] = chempro_out['pollutant'].str.replace('_','.')
    chempro_out['pollutant'] = chempro_out['pollutant'].str.upper()
    crossref_out = pd.read_csv(crossref, sep = None, engine='python', dtype=str).fillna(0)

    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    print('')
    print('******************************************************************')
    print('*****           Reading chemical profile is done             *****')
    print('******************************************************************')
    print('')    

    return Chemical_Speciation_table(chempro_out, crossref_out)




# =============================================================================
# Function to read the Speed average distribution
# =============================================================================
def read_avgSpeedDistribution(input_dir, avg_Speed_Distribution_file):
    start_time = time.time()
    print('******************************************************************')
    print('*****             Reading Speed profile                      *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    input_dir = input_dir
    spd_file = avg_Speed_Distribution_file #avg_SPD_Dist_file #
    name = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',spd_file)
    if os.path.exists(name) == True:
        print ('')
        print ('Reading Average Speed Distribution table ...')
        print (name)

        # with open(name, 'r') as csvfile:
        #     sep = csv.Sniffer().sniff(csvfile.read(4096)).delimiter

        ASD = pd.read_csv(name, sep = None, engine='python').fillna(np.nan)
        ASD.columns    = ASD.columns.str.lower()
        # for ird  in [101, 102, 103, 104, 105, 106, 107, 108]:
        #     ASD.loc[:,str(ird)] = ASD.loc[:,str(ird)] * ASD.loc[:,'speed']
        #     ASD.loc[:,str(ird)] = ASD.loc[:,str(ird)] / ASD.loc[:,str(ird)].sum()
    else:
        print ('')
        print('*** ERROR ABORT ***:  Average speed Distribution ', spd_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read Average speed Distribution')

    out_spd_bins = pd.DataFrame({'spd_bins': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]})
    out_spd      = pd.DataFrame({'spd': list(ASD.speed)})

    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    print('')
    print('******************************************************************')
    print('*****            Speed profile profile is done               *****')
    print('******************************************************************')
    print('')    

    return EF_Speed_Distribution(ASD, out_spd, out_spd_bins)



# =============================================================================
# Function to read the Activity Data
# =============================================================================
def read_activity_data_csv_SK(input_dir, ad_file, End_date, future_case = 'NO'):
    start_time = time.time()
    print('******************************************************************')
    print('*****           Processing Activity data                     *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    ad_file = ad_file #activity_file #
    end_date = End_date #ENDATE #
    future_case = future_case
#    now = dt.datetime.now()
#    ad_file = activity_file
#    end_date = ENDATE
    
    if len(end_date.split('-')) < 3:
        print ('')
        print ('***** ERROR *****: Please, Start date and End Data should be written as YYYY-MM-DD')
        print ('')
        sys.exit()
    else:
        current_yr = int(end_date.split('-')[0])

    name = '{0}{1}{2}'.format(input_dir,'/activity_data/',ad_file)
    with open(name, 'r') as csvfile:
        sep = csv.Sniffer().sniff(csvfile.read(4096)).delimiter
        
    outdf = pd.DataFrame(columns=['fuel','vehicle','engine','daily_vkt','region_cd',
                                  'manufacture_date'])
    count_outdf = pd.DataFrame()
    if os.path.exists(name) == True:
        print ('')
        print ('Reading Activity Data table ...')
        print (name)
        vhc_fuels = pd.DataFrame()
        vhc_cars  = pd.DataFrame()
        for activity_data in pd.read_csv(name, chunksize=1500000, sep = sep, usecols = [0,1,2,3,4,5]):
            print(activity_data.shape)
            activity_data = activity_data.fillna(0)
#        activity_data = pd.read_csv(name, chunksize=1500000, sep = sep, usecols = [0,1,2,3,4,5])
            '''
            header
               0,      1,    2,        3,          4,               5,
            Fuel,vehicle,engine,Daily_VKT,Region_code,Manufacture_date,
            '''
            activity_data.columns    = activity_data.columns.str.lower()
            activity_data.columns = ['fuel','vehicle','engine','daily_vkt','region_cd',
                                     'manufacture_date']

            activity_data.loc[:,'vehicle'] = activity_data.loc[:,'vehicle'].str.lower()
            activity_data.loc[:,'fuel']    = activity_data.loc[:,'fuel'].str.lower()
            activity_data.loc[:,'engine']   = activity_data.loc[:,'engine'].str.lower()
            activity_data.loc[:,'manufacture_date']   = (activity_data.loc[:,'manufacture_date'] / 10000).astype(int)
            activity_data['fullname'] = activity_data.vehicle.str.cat(activity_data[['engine','fuel']], sep='_')
            activity_data['count'] = 1
            
            vhc_fuels = vhc_fuels.append(pd.DataFrame({'fuels'  : list(activity_data['fuel'].unique())}),ignore_index=True, sort=False)
            vhc_cars  = vhc_cars.append(pd.DataFrame({'vhc'  : list(activity_data['vehicle'].unique())}),ignore_index=True, sort=False)
            if (activity_data.loc[:,'manufacture_date'].max() > current_yr) and \
            (future_case == 'NO') or (future_case == 'no') or (future_case == 'N') or (future_case == 'n'):
                AD_max_yr = activity_data.loc[:,'manufacture_date'].max()
                idx_drop = list(activity_data.loc[activity_data['manufacture_date'] > current_yr].index)
                print('***** WARNING *****')
                print('*** There are Activity data for future years . Check your input data ***')
                print('Current year {0}   :   Activity data max year {1}'.format(current_yr,AD_max_yr))
                print('*** Deleting these years ... ***')
                activity_data = activity_data.drop(index=idx_drop).reset_index(drop=True)
            
            if activity_data.loc[:,'daily_vkt'].min() < 0:
                print('***** WARNING *****')
                print('*** There are zero/null Activity data. Check your input data ***')
                print('*** Deleting these data ... ***')
                idxneg_drop = list(activity_data.loc[activity_data['daily_vkt'] < 0].index)
                activity_data = activity_data.drop(index=idxneg_drop).reset_index(drop=True)

            count_vhc = activity_data.groupby(['region_cd','fullname','manufacture_date']).sum().reset_index(drop=False)
            count_vhc = count_vhc.drop(columns=['daily_vkt'])
            count_vhc = count_vhc.sort_values(by=['region_cd','fullname','manufacture_date'])
            count_vhc = count_vhc.pivot_table(values='count',
                                              index=['region_cd','manufacture_date'],
                                              columns='fullname', 
                                              aggfunc=np.sum, fill_value=0).reset_index(drop=False)
            
            if count_outdf.shape[0] == 0:
                count_outdf = count_vhc
            else:
                count_outdf = count_outdf.append(count_vhc, sort=True)
                
            activity_data = activity_data.drop(columns=['count'])
            grouped = activity_data.groupby(['manufacture_date','region_cd','fullname']).sum().reset_index(drop=False)
            grouped = grouped.pivot_table(values='daily_vkt', 
                                          index=['manufacture_date','region_cd'], 
                                          columns='fullname', 
                                          aggfunc=np.sum).reset_index(drop=False).fillna(0.0)
            grouped = grouped.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
            outdf = outdf.append(grouped, ignore_index=True, sort=False)
#            vhc_names = pd.DataFrame({'vhc_name'  : list(activity_data.FullName.unique())})
#            vhc_years = pd.DataFrame({'vhc_years' : list(grouped.manufacture_date.unique())})
    #        out_table = Activity_Data_table(grouped, vhc_names,vhc_years)
    else:
        print('*** ERROR ABORT ***:  Emissions Factor file "" ', ad_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read Emissions Factor file')

    count_outdf = count_outdf.groupby(['region_cd','manufacture_date']).sum().reset_index(drop=False)
    count_outdf = count_outdf.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
    
    outdf = outdf.groupby(['manufacture_date','region_cd']).sum().reset_index(drop=False)
    outdf = outdf.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
    
    nvhc = list(np.sort(outdf.columns))
    nvhc.remove('manufacture_date')
    nvhc.remove('region_cd')
    vhc_names = pd.DataFrame({'vhc_name'  : nvhc})
    vhc_years = pd.DataFrame({'vhc_years' : list(outdf.manufacture_date.unique())})
    
    vhc_fuels = pd.DataFrame({'fuels' : list(vhc_fuels.fuels.unique())})
    vhc_cars  = pd.DataFrame({'vhc'   : list(vhc_cars.vhc.unique())})
    # There are "generic" counties condes endind with 0000. So we aggregate it
    # and split into the counties wich starts with the same four digits based on
    # the vehicle with the greatest VKT

    outdf.loc[:,'region_cd'] = outdf.loc[:,'region_cd'].astype(str)
    gendf_list_drop = outdf.loc[(outdf.region_cd.str.endswith(('0000','000')))].index.to_list()
    print('')
    print('***** WARNING *****: There counties with generic entry, meaning it ends with "0000" ')
    print('***** WARNING *****: These VKT values will be assigned into counties with similar region code')
    print('')
    gendf = outdf.loc[outdf.region_cd.str.endswith(('0000','000'))].reset_index(drop=True)
    gendf['region_cd_4dig'] = gendf.loc[:,'region_cd'].str.slice(stop=4)
    gendf_years = list(gendf.manufacture_date.unique())
    outdf   = outdf.drop(index=gendf_list_drop).reset_index(drop=True)
    
    sumdict = outdf.loc[:,nvhc].sum().to_dict()
    vhc_max_vkt = max(sumdict, key=sumdict.get) #getting the name og the vehivle with greatest VKT
    # Creating the split factor based on the VKT of vehicle vhc_max_vkt 
    splitdf  = outdf.loc[outdf.manufacture_date.isin(gendf_years),['manufacture_date','region_cd',vhc_max_vkt]].reset_index(drop=True)
    splitdf  = splitdf.rename(columns={vhc_max_vkt : 'split_factor'})
    splitdf['region_cd_4dig'] = splitdf.loc[:,'region_cd'].str.slice(stop=4)
    auxsplit = splitdf.groupby(['region_cd_4dig']).sum().reset_index(drop=False)
    auxsplit = auxsplit.rename(columns={'split_factor' : 'total_by_cnt'})
    auxsplit = auxsplit.drop(columns=['manufacture_date'])
    splitdf  = pd.merge(splitdf, auxsplit, right_on='region_cd_4dig', left_on='region_cd_4dig', how='left')
    splitdf.loc[:,'split_factor'] = (splitdf.loc[:,'split_factor'] / splitdf.loc[:,'total_by_cnt']).fillna(0.0)
    
    # Getting the total VKT by the "generic" county code
    gendf = gendf.groupby(['region_cd_4dig']).sum().reset_index(drop=False)
    gendf = gendf.drop(columns=['manufacture_date'])
    gendf = pd.merge(gendf, splitdf.loc[:,['manufacture_date','region_cd','split_factor','region_cd_4dig']],
                   left_on  ='region_cd_4dig', right_on ='region_cd_4dig', how='right').fillna(0.0)
    # Applying the split factor after merge the df
    gendf.loc[:,nvhc] = gendf.loc[:,nvhc].apply(lambda x: np.asarray(x) * gendf.split_factor.values).fillna(0.0)
    gendf = gendf.drop(columns=['region_cd_4dig', 'split_factor'])
    
    gendf = pd.merge(outdf.loc[:,['manufacture_date', 'region_cd']], gendf, 
                   left_on  = ['manufacture_date', 'region_cd'],
                   right_on = ['manufacture_date', 'region_cd'], how='left').fillna(0.0)
    
    # Summing the splitted "Generic VKT" into the output dataframe
    outdf.loc[:,nvhc] = outdf.loc[:,nvhc].add(gendf.loc[:,nvhc], axis=1)
    outdf = outdf.groupby(['manufacture_date', 'region_cd']).sum().reset_index(drop=False)

    # ------------------------------------------------------------------------#
    # Doing the same split to the vehicle count
    count_outdf.loc[:,'region_cd'] = count_outdf.loc[:,'region_cd'].astype(str)
    count_list_drop = count_outdf.loc[count_outdf.region_cd.str.endswith('0000')].index.to_list()
    print('')
    print('***** WARNING *****: There counties with generic entry, meaning it ends with "0000" ')
    print('***** WARNING *****: These VKT values will be assigned into counties with similar region code')
    print('')
    count_gendf = count_outdf.loc[count_outdf.region_cd.str.endswith('0000')].reset_index(drop=True)
    count_gendf['region_cd_4dig'] = count_gendf.loc[:,'region_cd'].str.slice(stop=4)
    count_outdf   = count_outdf.drop(index=count_list_drop).reset_index(drop=True)
    
    sumdict = count_outdf.loc[:,nvhc].sum().to_dict()
    vhc_max_vkt = max(sumdict, key=sumdict.get) #getting the name og the vehivle with greatest VKT
    # Creating the split factor based on the VKT of vehicle vhc_max_vkt 
    count_splitdf  = count_outdf.loc[:,['manufacture_date','region_cd',vhc_max_vkt]]
    count_splitdf  = count_splitdf.rename(columns={vhc_max_vkt : 'split_factor'})
    count_splitdf['region_cd_4dig'] = count_splitdf.loc[:,'region_cd'].str.slice(stop=4)
    count_auxsplit = count_splitdf.groupby(['region_cd_4dig']).sum().reset_index(drop=False)
    count_auxsplit = count_auxsplit.rename(columns={'split_factor' : 'total_by_cnt'})
    count_auxsplit = count_auxsplit.drop(columns=['manufacture_date'])
    count_splitdf  = pd.merge(count_splitdf, count_auxsplit, right_on='region_cd_4dig', left_on='region_cd_4dig', how='left')
    count_splitdf.loc[:,'split_factor'] = (count_splitdf.loc[:,'split_factor'] / count_splitdf.loc[:,'total_by_cnt']).fillna(0.0)
    
    # Getting the total VKT by the "generic" county code
    count_gendf = count_gendf.groupby(['region_cd_4dig']).sum().reset_index(drop=False)
    count_gendf = count_gendf.drop(columns=['manufacture_date'])
    count_gendf = pd.merge(count_gendf, count_splitdf.loc[:,['manufacture_date','region_cd','split_factor','region_cd_4dig']],
                   left_on  ='region_cd_4dig', right_on ='region_cd_4dig', how='right').fillna(0.0)
    # Applying the split factor after merge the df
    count_gendf.loc[:,nvhc] = count_gendf.loc[:,nvhc].apply(lambda x: np.asarray(x) * count_gendf.split_factor.values).fillna(0.0)
    count_gendf = count_gendf.drop(columns=['region_cd_4dig', 'split_factor'])
    
    # Summing the splitted "Generic VKT" into the output dataframe
    count_outdf.loc[:,nvhc] = count_outdf.loc[:,nvhc].add(count_gendf.loc[:,nvhc], axis=1)
    count_outdf = count_outdf.groupby(['manufacture_date', 'region_cd']).sum().reset_index(drop=False)

    outdf.to_csv(inter_dir+'/'+'outdf_AD_{0}.csv'.format(case_name)        , sep=',', index=False)
    vhc_names.to_csv(inter_dir+'/'+'vhc_names_AD_{0}.csv'.format(case_name), sep=',', index=False)
    vhc_years.to_csv(inter_dir+'/'+'vhc_years_AD_{0}.csv'.format(case_name), sep=',', index=False)
    vhc_fuels.to_csv(inter_dir+'/'+'vhc_fuels_AD_{0}.csv'.format(case_name), sep=',', index=False)
    vhc_cars.to_csv(inter_dir+'/'+'vhc_cars_AD_{0}.csv'.format(case_name)  , sep=',', index=False)
    count_outdf.to_csv(inter_dir+'/'+'vhc_numbers_AD_{0}.csv'.format(case_name), sep=',', index=False)
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    print('')
    print('******************************************************************')
    print('*****          Processing Activity data  is done             *****')
    print('******************************************************************')
    print('')    
    return Activity_Data_table(outdf, vhc_names, vhc_years, vhc_fuels, vhc_cars, count_outdf)
## =============================================================================



# =============================================================================
# Function to read link level shapefile
# =============================================================================
def roads_grid_surrogate_inf(input_dir, File_Name, Link_Shape_att, GridDesc_file='', Radius_Earth='', outGridShape = 'yes', Unit_meters = True):
    start_time = time.time()
    print('******************************************************************')
    print('*****         Processing Link/Road shapefile                 *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    Link_Shape_att = link_shape_att
    Link_ID_attr       = Link_Shape_att[0]   #= Link_ID_attr
    Region_CD          = Link_Shape_att[1]   #= Region_Code
    Region_NM          = Link_Shape_att[2]   #= Region_Name
    RD_name_attr       = Link_Shape_att[3]   #= RD_name_attr 
    RD_type_attr       = Link_Shape_att[4]   #= RD_type_attr
    Speed_attr         = Link_Shape_att[5]   #= Speed_attr
    Link_length        = Link_Shape_att[6]   #= Link_length
    VKT_attr           = Link_Shape_att[7]   #= VKT_attr
    file_name          = File_Name #link_shape # 
    gridfile_name      = GridDesc_file #gridfile_name #
    outGridShape       = outGridShape
    if Radius_Earth != '':
        radius = Radius_Earth
    else:
        radius = 6370000.0
    if ((grid4AQM.lower() == 'yes') or (grid4AQM.lower() == 'y')) and (gridfile_name == ''):
        print('')
        print('*** ERROR ***: Error on processing grid to CMAQ')
        print('*** ERROR ***: Please set the gridfile_name and grid4AQM the correctly')
        print('')
        sys.exit()

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
        out_roads['geometry']      = out_roads.buffer(0.2) #changing link to area to do the overlay
        out_roads['total_area']    = out_roads.area
        out_roads['link_split_total']    = (out_roads['link_length'] * 0.0).astype(float)
        out_roads['link_split_county']    = (out_roads['link_length'] * 0.0).astype(float)
        out_roads['vkt_split_county']    = (out_roads['link_length'] * 0.0).astype(float)
        reduc = 0.6                                  #60% reduction as BH asked
        # rt = {101 : 80 *reduc, 102 : 60 *reduc, 103 : 60 *reduc, 104 : 50 *reduc,   #60% reduction as BH asked
        #       105 : 30 *reduc, 106 : 30 *reduc, 107 : 30 *reduc, 108 : 30 *reduc}

        # for key, ispd in rt.items():
        #     out_roads.loc[(out_roads.loc[:,'road_type'] == key),'max_speed'] = ispd
        # Creating grid
        if ((grid4AQM.lower() == 'yes') or (grid4AQM.lower() == 'y')) and (gridfile_name != ''):
            if os.path.exists(gridfile_name) == True:
                print ('')
                print ('Reading GridDesc file to generate gridded emissions ...')
                print ('')
                print (gridfile_name)
                griddesc = pd.read_csv(gridfile_name, header=None,engine='python')
                
                COORDTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT = griddesc.loc[2,0].split()
                P_ALP, P_BET, P_GAM, XCENT, YCENT = map(float, [P_ALP, P_BET, P_GAM, XCENT, YCENT])
                GRID_NAME = griddesc.loc[4,0].split()[0].replace("'",'')
                
                COORD_NAME, XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK = griddesc.loc[5,0].split()
                XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK = \
                map(float, [XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK])
                
                if   (COORDTYPE == '1'):  #lat-lon
                    proj_grid = pyproj.Proj(init='epsg:4326')
                elif (COORDTYPE == '2'):  #lambert
                    proj_grid = pyproj.Proj('+proj=lcc +ellps=sphere +lat_1={0} +lat_2={1} +lat_0={3} +lon_0={2} +a={4} +R={4} +no_defs'.format(P_ALP, P_BET, XCENT, YCENT,radius)) #lcc
                elif (COORDTYPE == '5') and (P_ALP <0):
                     proj_grid = pyproj.Proj('+proj=utm +zone={0} +south +datum=WGS84 +to_meter=1 +no_defs'.format(P_ALP, radius))
                elif (COORDTYPE == '5') and (P_ALP >= 0):
                     proj_grid = pyproj.Proj('+proj=utm +zone={0} +datum=WGS84 +to_meter=1 +no_defs'.format(P_ALP,radius))
                elif (COORDTYPE == '6'): #polar
                    proj_grid = pyproj.Proj('+proj=stere +lat_0={0} +lat_ts={1} +lon_0={2} +lat_0={3} +R={4} +ellps=WGS84 +datum=WGS84 +to_meter=1 +no_defs'.format(P_ALP, P_BET, P_GAM, XCENT, YCENT, radius))

                cols =  int(NCOLS)
                rows =  int(NROWS)
                print("")
                print('***** The domain has {0} columns and {1} rows '.format(cols,rows))
                print("")
                proj_grid(XORIG, YORIG, inverse=True)
                proj_grid(XCENT, YCENT, inverse=True)
                
                xmin, ymin = [(float(XORIG)), (float(YORIG))]
                xmax, ymax = [xmin + ((cols+1) * XCELL), ymin + ((rows+1) * YCELL)]
                proj_grid(xmax,ymax, inverse=True)
                lat_list = (np.arange(ymin,ymax,YCELL))
                lon_list = (np.arange(xmin,xmax,XCELL))
                polygons = []
                grid_row = []
                grid_col = []
                print('')
                print('Creating grid cells ...')
                print('')
                for j in range(0,rows):
                    yini = lat_list[j]
                    yfin = lat_list[j+1] 
                    for i in range(0,cols):
                        grid_row.append(j+1)
                        grid_col.append(i+1)
                        xini = lon_list[i]
                        xfin = lon_list[i+1]
                        polygons.append(Polygon([(xini, yini), (xfin, yini),
                                                 (xfin, yfin), (xini, yfin), (xini, yini)])) 

                grid_ID = [x for x in range (1,len(polygons)+1)]
                grid = gpd.GeoDataFrame({'geometry': polygons, 'grid_id': grid_ID,
                                         'row'     : grid_row, 'col'    : grid_col})
                grid.crs = proj_grid.srs
                
                # exporting grid as shapefile
                if (outGridShape == 'yes') or (outGridShape == 'YES') or \
                (outGridShape == 'y')   or (outGridShape == 'Y'):
                     grid.to_file(filename = output_dir+'/grid_{0}.shp'.format(case_name),
                                  driver='ESRI Shapefile', crs_wkt=grid.crs)
                
                out_roads = out_roads.to_crs(proj_grid.srs)
            else:
                print ('')
                print ('There is no GridCro File in {0}'.format(gridfile_name))
                sys.exit()
        else:
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
            print('')
            print('Creating grid cells ...')
            print('')
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
            if (outGridShape.lower() == 'yes') or (outGridShape.lower() == 'y'):
                grid.to_file(filename = output_dir+'/grid_{0}.shp'.format(case_name), driver='ESRI Shapefile',crs_wkt=prj)

        print('')
        print('Calculating the surrograte based on the link shapefile ...')
        print('')
        #creating the surrogate
        surrogate = gpd.overlay(out_roads, grid, how='intersection').reset_index(drop=True)
        surrogate['split_area'] = surrogate.area
        surrogate['vkt_norm'] = ((surrogate.loc[:,'vkt_avg'] * \
                                  surrogate.loc[:,'split_area']) / \
                                  surrogate.loc[:,'total_area']).astype(float)
        
        surrogate['weight_factor'] = surrogate['vkt_norm'] * 0.0

        for igeocd in out_roads.region_cd.unique():
            srgt_split_vkt = surrogate.vkt_norm.loc[surrogate.region_cd == igeocd].values / \
                            (surrogate.vkt_norm.loc[surrogate.region_cd == igeocd]).sum()
            surrogate.loc[surrogate.region_cd == igeocd, ['weight_factor']] = srgt_split_vkt
            
            vkt_split_county = out_roads.vkt_avg.loc[out_roads.region_cd == igeocd].values / \
                              (out_roads.vkt_avg.loc[out_roads.region_cd == igeocd]).sum()
            out_roads.loc[out_roads.region_cd == igeocd, ['vkt_split_county']] = vkt_split_county

            lnk_split_county = out_roads.link_length.loc[out_roads.region_cd == igeocd].values / \
                              (out_roads.link_length.loc[out_roads.region_cd == igeocd]).sum()
            out_roads.loc[out_roads.region_cd == igeocd, ['link_split_county']] = lnk_split_county


        surrogate = surrogate.groupby(['region_cd','grid_id']).sum()
        surrogate = surrogate.reset_index(drop=False)
        surrogate = surrogate.loc[:,['region_cd', 'grid_id', 'link_length',
                                     'vkt_avg', 'total_area', 'split_area',
                                     'link_split_total','link_split_county',
                                     'vkt_split_county', 'weight_factor']]
        
        surrogate = pd.merge(grid, surrogate, how='left', on=['grid_id'])
        surrogate.loc[:,~surrogate.columns.str.contains('geometry')] = \
            surrogate.loc[:,~surrogate.columns.str.contains('geometry')].fillna(0)
        out_roads = out_roads.drop(columns=['geometry'])
        out_roads = out_roads.rename(columns={'geometry_BKP': 'geometry'}).set_geometry('geometry')
        out_roads['region_cd'] = out_roads['region_cd'].astype(int)
        surrogate['region_cd'] = surrogate['region_cd'].astype(int)
    else:
        print('*** ERROR ABORT ***:  Shapefile "" ', shp_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read link Shapefile file')
    crs = out_roads.crs
    if isinstance(crs, dict):
        crs = pyproj.Proj(out_roads.crs).srs
    crs_nameOut = inter_dir+'/'+'out_roads_CRS_{0}.txt'.format(case_name)
    with open(crs_nameOut, 'w') as crsfile:
        crsfile.write(crs)

    grid.to_csv(inter_dir+'/'+'grid_{0}.csv'.format(case_name)        , sep=',', index=False)
    surrogate.to_csv(inter_dir+'/'+'surrogate_{0}.csv'.format(case_name)        , sep=',', index=False)
    out_roads.to_csv(inter_dir+'/'+'out_roads_{0}.csv'.format(case_name)        , sep=',', index=False)
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    print('')
    print('******************************************************************')
    print('*****      Processing Link/Road shapefile is done            *****')
    print('******************************************************************')
    print('')    

    return Roads_Grid_table(grid, surrogate, out_roads, crs)
# =============================================================================
    





# =============================================================================
# Function to read link level shapefile
# =============================================================================
def processing_County_shape(input_dir, file_name, County_Shape_att, GridDesc_file='', Radius_Earth=''):
    start_time = time.time()
    print('******************************************************************')
    print('*****           Processing county shapefile                  *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    Region_Geocode       = County_Shape_att[0] # county_shape_att[0] #
    Region_name_attr     = County_Shape_att[1] # county_shape_att[1] # 
    Region_name_attr_SK  = County_Shape_att[2] # county_shape_att[2] #
#    Link_ID_attr       = 'LINK_ID'
#    RD_name_attr       = 'ROAD_NAME'
#    RD_type_attr       = 'ROAD_RANK'
#    Activity_data_attr = 'SHAPE_STLe'
#    Speed_attr         = 'MAX_SPD'
#    Link_length        = 'SHAPE_STLe'
    file_name          = county_shape
    if Radius_Earth != '':
        radius = Radius_Earth
    else:
        radius = 6370000.0

    if ((grid4AQM.lower() == 'yes') or (grid4AQM.lower() == 'y')) and (gridfile_name != ''):
        if os.path.exists(gridfile_name) == True:
            print ('')
            print ('Reading GridDesc file to generate gridded emissions ...')
            print (gridfile_name)
            griddesc = pd.read_csv(gridfile_name, header=None,engine='python')
            
            COORDTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT = griddesc.loc[2,0].split()
            P_ALP, P_BET, P_GAM, XCENT, YCENT = map(float, [P_ALP, P_BET, P_GAM, XCENT, YCENT])
            GRID_NAME = griddesc.loc[4,0].split()[0].replace("'",'')
            
            COORD_NAME, XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK = griddesc.loc[5,0].split()
            XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK = \
            map(float, [XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK])
            
            if   (COORDTYPE == '1'):  #lat-lon
                proj_grid = pyproj.Proj(init='epsg:4326')
            elif (COORDTYPE == '2'):  #lambert
                proj_grid = pyproj.Proj('+proj=lcc +ellps=sphere +lat_1={0} +lat_2={1} +lat_0={3} +lon_0={2} +a={4} +R={4} +no_defs'.format(P_ALP, P_BET, XCENT, YCENT,radius)) #lcc
            elif (COORDTYPE == '5') and (P_ALP <0):
                 proj_grid = pyproj.Proj('+proj=utm +zone={0} +south +datum=WGS84 +to_meter=1 +no_defs'.format(P_ALP, radius))
            elif (COORDTYPE == '5') and (P_ALP >= 0):
                 proj_grid = pyproj.Proj('+proj=utm +zone={0} +datum=WGS84 +to_meter=1 +no_defs'.format(P_ALP,radius))
            elif (COORDTYPE == '6'): #polar
                proj_grid = pyproj.Proj('+proj=stere +lat_0={0} +lat_ts={1} +lon_0={2} +lat_0={3} +R={4} +ellps=WGS84 +datum=WGS84 +to_meter=1 +no_defs'.format(P_ALP, P_BET, P_GAM, XCENT, YCENT, radius))

            cnty_file = '{0}{1}'.format(input_dir,file_name)
            if os.path.exists(cnty_file) == True:
                print ('')
                print ('Reading County Shapefile to get bounds to create the GRID...')
                print (cnty_file)
                cnty_shp = gpd.read_file(cnty_file)
                cnty_shp = cnty_shp.rename(columns={Region_Geocode       : 'region_cd'})
                cnty_shp = cnty_shp.rename(columns={Region_name_attr     : 'region_nm'})
                cnty_shp = cnty_shp.rename(columns={Region_name_attr_SK  : 'region_nm_SK'})
            
                out_cnty = cnty_shp.loc[:,['geometry','region_cd','region_nm', 
                                           'region_nm_SK']]
                
                out_cnty['region_cd'] = out_cnty['region_cd'].astype(int)
                out_cnty = out_cnty.to_crs(proj_grid.srs)

    cnty_file = '{0}{1}'.format(input_dir,file_name)
    if os.path.exists(cnty_file) == True:
        print('******************************************************************')
        print('*****               Reading County Shapefile ...             *****')
        print('*****                    Please wait ...                     *****')
        print('******************************************************************')
        print('')
        print (cnty_file)

        cnty_shp = gpd.read_file(cnty_file)
        cnty_shp = cnty_shp.rename(columns={Region_Geocode       : 'region_cd'})
        cnty_shp = cnty_shp.rename(columns={Region_name_attr     : 'region_nm'})
        cnty_shp = cnty_shp.rename(columns={Region_name_attr_SK  : 'region_nm_SK'})
    
        out_cnty = cnty_shp.loc[:,['geometry','region_cd','region_nm', 
                                   'region_nm_SK']]
        
        out_cnty['region_cd'] = out_cnty['region_cd'].astype(int)

    else:
        print('*** ERROR ABORT ***:  Shapefile "" ', cnty_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read county Shapefile file')         
    
    
    out_cnty.to_csv(inter_dir+'/'+'out_district_{0}.csv'.format(case_name), sep=',', index=False)
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    print('')
    print('******************************************************************')
    print('*****        Processing county shapefile is done             *****')
    print('******************************************************************')
    print('')    

    return out_cnty
# =============================================================================
    



# =============================================================================
# Function to read the emissions factor from South Korea* (*particular case)
# =============================================================================
def read_emissions_factor_SK(input_dir, EmisFactor_list):
    start_time = time.time()
    print('******************************************************************')
    print('*****         Reading Emissions Factors tables ...           *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    input_dir = input_dir
    ef_list = EmisFactor_list #ef_file #['gasoline.csv'] #
#    ef_list = ['diesel_v3.csv']
#    sep = ';'
    final_EF = pd.DataFrame()
    for ifile in range(0,len(ef_list)):
#        print(ef_list[ifile])
        name = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',ef_list[ifile])
        if os.path.exists(name) == True:
            print('')
            print (name)
            with open(name, 'r') as csvfile:
                sep = csv.Sniffer().sniff(csvfile.read(4096)).delimiter

            emissions_factor = (pd.read_csv(name, sep = sep)).fillna(0)
            emissions_factor.columns = emissions_factor.columns.str.lower()
            if (final_EF.shape[0] == 0) and (final_EF.shape[1] == 0):
                final_EF = pd.DataFrame(columns=emissions_factor.columns)
            iaux = 1
            drop_list = []
            aux = True
            while aux == True:
                pos = (emissions_factor.shape[0] - iaux)
                aux = 0 in list(emissions_factor.loc[pos,['vehicle','engine','scc','pollutant']])
                if aux == True:
                    drop_list.append(pos)
                iaux +=1
            emissions_factor = emissions_factor.drop(index=drop_list)
            emissions_factor['years']           =  emissions_factor['years'] / 10000
            emissions_factor['years']           =  emissions_factor.years.astype(int)
            emissions_factor.loc[:,'pollutant'] = emissions_factor.loc[:,'pollutant'].str.upper()
            emissions_factor.loc[:,'vehicle']   = emissions_factor.loc[:,'vehicle'].str.lower()
            emissions_factor.loc[:,'engine']     = emissions_factor.loc[:,'engine'].str.lower()
            emissions_factor.loc[:,'fuel']      = emissions_factor.loc[:,'fuel'].str.lower()
            emissions_factor['fullname']        = emissions_factor.vehicle.str.cat(emissions_factor[['engine','fuel']], sep='_')
            emissions_factor.loc[:,['a', 'b', 'c', 'd', 'f', 'k']] =  emissions_factor.loc[:,['a', 'b', 'c', 'd', 'f', 'k']].astype(float)
            
            vhc  = emissions_factor.fullname.unique()
            poll = list(emissions_factor.pollutant.unique())
            yr   = list(np.sort(emissions_factor.years.unique()))
    #        ivhc = 'sedan_supercompact_diesel'
    #        ipoll = 'CO'
    #        iyr = 1990
            log_name = 'LOG_read_Emissions_Factor_functions.txt'
            with open(output_dir+'/LOGS/'+log_name, 'w') as EF_log:
                for ivhc in vhc:
                    for ipoll in poll:
                        for iyr in yr:
                            dat_ef = emissions_factor[((emissions_factor['fullname']   == ivhc)  & 
                                                       (emissions_factor['pollutant'] == ipoll) &
                                                       (emissions_factor['years']     == iyr))] 
                            dat_ef = dat_ef.reset_index(drop=True)
                            iscc = (emissions_factor.scc[(emissions_factor.fullname == ivhc)].unique())
                            if (dat_ef.shape[0] == 0) and  (iscc.size != 0) :
                                EF_log.write('*** WARNING *** There is no emissions factor of year {0} for vehicle {1}\n'.format(iyr,ivhc))
                            elif (dat_ef.years.shape[0] == 2) and (dat_ef.v[0] == 0) and (dat_ef.temperatures[0] == 0):
                                df = dat_ef.iloc[[0]]
                                final_EF = final_EF.append(df,ignore_index=True, sort=False)

                            elif (dat_ef.years.shape[0] == 3) and (True in list(dat_ef.v.isin([0]))) and (dat_ef.temperatures[0] == 0):
                                df = dat_ef.loc[~dat_ef.v.isin([0])]
                                final_EF = final_EF.append(df,ignore_index=True, sort=False)

                            elif (dat_ef.years.shape[0] >= 4) and (dat_ef.v[0] != 0) and (dat_ef.temperatures[0] == 0):
                                df = dat_ef.iloc[[0,1]]
                                final_EF = final_EF.append(df,ignore_index=True, sort=False)
                            else:
                                df = dat_ef.loc[:]
                                final_EF = final_EF.append(df,ignore_index=True, sort=False)
    
        else:
            print (' ')
            print('*** ERROR ABORT ***:  Emissions Factor file "" ', ef_list[ifile], ' "" does not exist!')
            sys.exit('CARS preProcessor can not read Emissions Factor file')
    
    final_EF = final_EF.reset_index(inplace=False, drop=True)
    EF_names = pd.DataFrame({'EF_fullname': list(final_EF.fullname.unique())})
    EF_years = pd.DataFrame({'EF_years'   : list(np.sort(final_EF.years.unique()))})
    EF_fuels = pd.DataFrame({'EF_fuel'    : list(final_EF.fuel.unique())})
    EF_polls = pd.DataFrame({'EF_poll'    : list(final_EF.pollutant.unique())})
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    print('')
    print('******************************************************************')
    print('*****          Reading Emissions Factors is done             *****')
    print('******************************************************************')
    print('')    
    
    return EF_Grid_table(final_EF, EF_years, EF_names, EF_fuels, EF_polls)
# =============================================================================






# =============================================================================
# Function to calculate the emissions factor from South Korea
# =============================================================================                              
def calculate_EmisFact(Input_dir, Emissions_Factor_Table, Out_dir ):
    start_time = time.time()
    print('******************************************************************')
    print('*****          Calculating Emissions Factors ...             *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    input_dir = Input_dir
    output_dir = Out_dir
    EF_All_fuel = Emissions_Factor_Table.data.copy() #EF_All_fuel_function.data.copy() #
    print ('')
    print ('Calculating Emissions Factor from EF functions ...')
    print ('')
#    nspd = np.asarray([x for x in range(1,151)])
#    Spd Bins    :   1, 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,  14,  15, 16
#    Spd Bins mph: 2.5, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,  65,  70, 75
#    Spd Bins kmh:   4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 89, 97, 105, 113, 121
    nspd = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 89, 97, 105, 113, 121]
    EF_spd_limit = []
    for i in list(EF_All_fuel.v.unique()):
        if (type(i) == int) and (i != 0):
            print('')
            print('*** WARNING *** Check your Emissions Factor input at EF temp depended ')
            print('')
        elif type(i) != int:
            EF_spd_limit.append(i)
    EF_spd_limit = np.unique(EF_spd_limit)

    dat_EF = pd.DataFrame(columns=['fullname','pollutant','v','years','temperatures','spd','emis_fact'])
    for ispd in nspd:
        aux_EF = EF_All_fuel.loc[:,['fullname','pollutant','v','years','temperatures']]
        aux_EF['emis_fact'] = EF_All_fuel.k * ((EF_All_fuel.a * (ispd ** EF_All_fuel.b)) + \
                                               (EF_All_fuel.c * (ispd ** EF_All_fuel.d)) + \
                                               (EF_All_fuel.f))

        aux_EF['spd']       = ispd
        dat_EF = dat_EF.append(aux_EF, ignore_index=True, sort=False)
    
    dat_EF['ambient_temp'] = ((dat_EF.spd * 0 ) + temp_mean).astype(int)
    dat_EF['spd']          = dat_EF['spd'].astype(int)
    dat_EF['years']        = dat_EF['years'].astype(int)
    dat_EF['pollutant']    = dat_EF['pollutant'].str.upper()
    #Extracting the NOX Diesel EF from final EF table
    Temp_EF_Diesel = dat_EF[(dat_EF.temperatures != 0) & 
                            (dat_EF.pollutant == 'NOX')]
    Temp_EF_Diesel_drop_list = list(Temp_EF_Diesel.index)
    dat_EF = dat_EF.drop(index=Temp_EF_Diesel_drop_list).reset_index(drop=True)
    Temp_EF_Diesel = Temp_EF_Diesel.reset_index(drop=True)
    
    #Filtering the NOX Diesel EF based on the ambient temperature
    filter_temp_diesel = []
    for tlim in Temp_EF_Diesel.temperatures.unique():
        if len(tlim) <= 4:
            if (tlim[:2] == 'GT') and (temp_mean < int(tlim[2:])):
                idx = list(Temp_EF_Diesel[(Temp_EF_Diesel.temperatures  == tlim)].index)
                filter_temp_diesel.extend(idx)
            elif (tlim[:2] == 'LE') and (temp_mean > int(tlim[2:])):
                idx = list(Temp_EF_Diesel[(Temp_EF_Diesel.temperatures  == tlim)].index)
                filter_temp_diesel.extend(idx)
        elif len(tlim) > 4:
            tlim2 = tlim.split('_')
            if (temp_mean < int(tlim2[0][2:])) or (temp_mean > int(tlim2[1][2:])):
                idx = list(Temp_EF_Diesel[(Temp_EF_Diesel.temperatures  == tlim)].index)
                filter_temp_diesel.extend(idx)
                
    Temp_EF_Diesel = Temp_EF_Diesel.drop(index=filter_temp_diesel).reset_index(drop=True)
    dat_EF = dat_EF.append(Temp_EF_Diesel, ignore_index=True, sort=True)
    
    #The EF were calculated for all speeds for all equations - Filtering the EF
    # based of the equation and speed 
    drop_list = []
    for vlim in EF_spd_limit:
        if (len(vlim.split('_')) == 1):
            if vlim[:2] == 'GT':
                idx = list(dat_EF[(dat_EF.v   == vlim) & 
                                  (dat_EF.spd <= int(vlim[2:]))].index)
                drop_list.extend(idx)
            elif vlim[:2] == 'LE':
                idx = list(dat_EF[(dat_EF.v   == vlim) & 
                                  (dat_EF.spd > int(vlim[2:]))].index)
                drop_list.extend(idx)
        
        elif (len(vlim.split('_')) == 2):
            if (vlim.split('_')[0][:2] == 'GT') and (vlim.split('_')[1][:2] == 'LE'):
                idGT = list(dat_EF.loc[(dat_EF.v   == vlim) & 
                                      (dat_EF.spd <= int(vlim.split('_')[0][2:]))].index)
                drop_list.extend(idGT)
                
                idLE = list(dat_EF.loc[(dat_EF.v   == vlim) & 
                                      (dat_EF.spd > int(vlim.split('_')[1][2:]))].index)
                drop_list.extend(idLE)
            
    dat_EF = dat_EF.drop(index=drop_list).reset_index(drop=True)
    
    #Filtering the duplicate EF of NOX Diesel - we are averaging them
    Temp_EF_Diesel = dat_EF[(dat_EF.temperatures != 0) & 
                            (dat_EF.pollutant == 'NOX')]
    Temp_EF_Diesel_drop_list = list(Temp_EF_Diesel.index)
    dat_EF = dat_EF.drop(index=Temp_EF_Diesel_drop_list).reset_index(drop=True)
    Temp_EF_Diesel = Temp_EF_Diesel.groupby(['fullname','pollutant','temperatures','v','years','ambient_temp','spd']).mean()
    Temp_EF_Diesel = Temp_EF_Diesel.reset_index(drop=False)
    dat_EF = dat_EF.append(Temp_EF_Diesel, ignore_index=True, sort=True)
    dat_EF = dat_EF.sort_values(by=['fullname','spd'])
    dat_EF = dat_EF.reset_index(drop=True)
    dat_EF.loc[(dat_EF.emis_fact < 0),'emis_fact'] = 0
    dat_EF.loc[:,['ambient_temp','emis_fact']] = dat_EF.loc[:,['ambient_temp','emis_fact']].astype(float)
    dat_EF.loc[:,['spd','years']] = dat_EF.loc[:,['spd','years']].astype(int)

    EF_names = pd.DataFrame({'EF_fullname': list(np.sort(dat_EF.fullname.unique()))})
    EF_years = pd.DataFrame({'EF_years'   : list(np.sort(dat_EF.years.unique()))})
    EF_fuels = pd.DataFrame({'EF_fuel'    : list(EF_All_fuel.fuel.unique())})
    EF_polls = pd.DataFrame({'EF_poll'    : list(np.sort(dat_EF.pollutant.unique()))})
    
    dat_EF.to_csv(output_dir+'/'+'EmisFact_by_YR_SPD.csv',sep=',', index=False)
    
    dat_EF.to_csv(inter_dir+'/'+'dat_EF_{0}.csv'.format(case_name)      , sep=',', index=False)
    EF_years.to_csv(inter_dir+'/'+'EF_years_{0}.csv'.format(case_name) , sep=',', index=False)
    EF_names.to_csv(inter_dir+'/'+'EF_names_{0}.csv'.format(case_name) , sep=',', index=False)
    EF_fuels.to_csv(inter_dir+'/'+'EF_fuels_{0}.csv'.format(case_name) , sep=',', index=False)
    EF_polls.to_csv(inter_dir+'/'+'EF_polls_{0}.csv'.format(case_name) , sep=',', index=False)
    
    EmisFact_yr_spd = EF_Grid_table(dat_EF, EF_years, EF_names, EF_fuels, EF_polls)
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    print('')
    print('******************************************************************')
    print('*****       Calculating Emissions Factors is done            *****')
    print('******************************************************************')
    print('')    
    
    return EmisFact_yr_spd









# =============================================================================
# Function to check any missing emissions factor
# =============================================================================
def check_VHC_AD_EF(Emissions_Factor_DataFrame, Activity_Data_DataFrame, Out_dir):
    start_time = time.time()
    print('******************************************************************')
    print('*****        Checking Vehicles and Activity data             *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    EFdf = Emissions_Factor_DataFrame
    ADdf = Activity_Data_DataFrame
    OutD = Out_dir
    try:
        ADdf
        a = True
    except: a = False
    
    try:
        EFdf
        b = True
    except: b = False
    
    if (a == True) and (type(ADdf) is (Activity_Data_table)) and \
       (b == True) and (type(EFdf) is (EF_Grid_table)):
        EFvhc = EFdf.EF_fullname.EF_fullname.tolist()
        ADvhc = ADdf.fullname.vhc_name.tolist()
        print('')
        print('******************************************************************')
        print('*****      Checking if all Vehicles have Emission Factors      ***')
        print('*****            Please wait and check the LOG...            *****')
        print('******************************************************************')
        print('')
        log_name = 'LOG_Missing_Vehicle_EmissionFactor.txt'
        with open(OutD+'/LOGS/'+log_name, 'w') as log:
            log.write('# ***** Vehicles without Emissions Factors *****\n')
            log.write('# \n')
            log.write('vehicle_engine_fuel\n')
            for ivhc in ADvhc:
                if ivhc not in EFvhc:
                    log.write('{0}\n'.format(ivhc))
    else:
        print('')
        print('*** ERROR ABORT ***: Check if the Activity data and Emissions Factor were processed')
        print('')
        sys.exit()
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    print('')
    print('******************************************************************')
    print('*****      Checking Vehicles and Activity data is done       *****')
    print('******************************************************************')
    print('')    

    return       



# =============================================================================
# Function to apply Deterioration rate into the emissions factor from South Korea
# =============================================================================                              

def apply_deterioration_emissions_factor_SK(input_dir, Deterioration_list, EmisFactor_yr_spd, output_dir):
    start_time = time.time()
    print('') 
    print('******************************************************************')
    print('*****         Deterioration_factor set to YES ...              *****')
    print('*****        Applying Deterioration factor to EF             *****')
    print('******************************************************************')
    print('')
    print('***** WARNING: The pure emissions factor will be discarded and only \
          the deteriorated EF will be used, unless there is no deterioration for an \
          specific vehicle or/and pollutant *****')    

    input_dir = input_dir
    output   = output_dir
    def_list = Deterioration_list #ef_file #['gasoline.csv'] #
    EF_yr_spd  = EmisFactor_yr_spd
#        def_list = ['degradation_rate_LPG.csv'] #['degradation_rate_Diesel.csv']
    final_DEF = pd.DataFrame()
    for ifile in range(0,len(def_list)):
#        print(ef_list[ifile])
        name = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',def_list[ifile])
        if os.path.exists(name) == True:
            print('')
            print('******************************************************************')
            print('*****      Reading Deterioration factor tables ...           *****')
            print('*****                    Please wait ...                     *****')
            print('******************************************************************')
            print('')
            print (name)
            with open(name, 'r') as csvfile:
                sep = csv.Sniffer().sniff(csvfile.read(4096)).delimiter

            det_factor = (pd.read_csv(name, sep = sep)).fillna(0)
            det_factor.columns = det_factor.columns.str.lower()
            det_factor = det_factor.rename(columns={'manufacture date' : 'manufacture_date'})
            if (final_DEF.shape[0] == 0) and (final_DEF.shape[1] == 0):
                final_DEF = pd.DataFrame(columns=det_factor.columns)
            iaux = 1
            drop_list = []
            aux = True
            while aux == True:
                pos = (det_factor.shape[0] - iaux)
                aux = 0 in list(det_factor.loc[pos,['vehicle','engine','scc','pollutant']])
                if aux == True:
                    drop_list.append(pos)
                iaux +=1
            det_factor = det_factor.drop(index=drop_list)
            det_factor['manufacture_date'] =  det_factor['manufacture_date'] / 10000
            det_factor['manufacture_date'] =  det_factor['manufacture_date'].astype(int)
            det_factor.loc[:,'vehicle'] = det_factor.loc[:,'vehicle'].str.lower()
            det_factor.loc[:,'engine'] = det_factor.loc[:,'engine'].str.lower()
            det_factor.loc[:,'fuel'] = det_factor.loc[:,'fuel'].str.lower()
            det_factor['fullname'] = det_factor.vehicle.str.cat(det_factor[['engine','fuel']], sep='_')
            final_DEF = final_DEF.append(det_factor, ignore_index=True, sort=False)
    
    final_DEF = final_DEF.sort_values(by=['manufacture_date','fullname']).reset_index(drop=True)
    final_DEF = final_DEF.set_index(['fullname','pollutant','manufacture_date'])
    final_DEF = final_DEF.drop(columns=['vehicle','engine','fuel','scc'])
    final_DEF = final_DEF.stack().reset_index(drop=False)
    final_DEF = final_DEF.rename(columns={'level_3' : 'age', 0: 'det_fac'})
    final_DEF.loc[:,['age', 'det_fac']] = final_DEF.loc[:,['age', 'det_fac']].astype(float)
    final_DEF['years'] = final_DEF['manufacture_date'] - final_DEF['age']
    
    # Droping the duplicantes by applying the average values
    final_DEF = final_DEF.groupby(['fullname', 'pollutant','years']).mean().reset_index(drop=False)
    
    def_vhc = final_DEF.fullname.unique().tolist()
    dat_EF = EF_yr_spd.data.loc[EF_yr_spd.data.fullname.isin(def_vhc)]
    drop_DEF_list = dat_EF.index.tolist()
    dat_EF = dat_EF.reset_index(drop=True)
    dat_EF.loc[:,['ambient_temp','emis_fact']] = dat_EF.loc[:,['ambient_temp','emis_fact']].astype(float)
    dat_EF.loc[:,['spd','years']] = dat_EF.loc[:,['spd','years']].astype(int)
    
    aux_EF_DEF = pd.merge(dat_EF, final_DEF, left_on=['fullname', 'pollutant','years'],
                   right_on = ['fullname', 'pollutant','years'], how='left')

    aux_EF_DEF = aux_EF_DEF.groupby(['fullname', 'pollutant','spd']).apply(lambda x: x.ffill().bfill())
    aux_EF_DEF = aux_EF_DEF.fillna(1.0)
    
    merge_EF_Final = EF_yr_spd.data.copy()
    merge_EF_Final = merge_EF_Final.drop(index=drop_DEF_list).reset_index(drop=True)
    merge_EF_Final = merge_EF_Final.append(aux_EF_DEF, ignore_index=True, sort=False)
    merge_EF_Final = merge_EF_Final.sort_values(by=['years','fullname']).reset_index(drop=True)
    merge_EF_Final = merge_EF_Final.fillna(1.0)
    merge_EF_Final['emis_fact'] = merge_EF_Final['emis_fact'] * merge_EF_Final['det_fac']
    
    out_cols = EF_yr_spd.data.columns.tolist()
    merge_EF_Final = merge_EF_Final.loc[:,out_cols + ['det_fac']]
    
    EF_names = pd.DataFrame({'EF_fullname': list(np.sort(merge_EF_Final.fullname.unique()))})
    EF_years = pd.DataFrame({'EF_years'   : list(np.sort(merge_EF_Final.years.unique()))})
    EF_fuels = pd.DataFrame({'EF_fuel'    : list(EF_yr_spd.EF_fuels.EF_fuel)})
    EF_polls = pd.DataFrame({'EF_poll'    : list(np.sort(merge_EF_Final.pollutant.unique()))})
    
    merge_EF_Final.to_csv(output_dir+'/'+'EmisFact_by_YR_SPD_Deteriorated.csv',sep=',', index=False)
    
    merge_EF_Final.to_csv(inter_dir+'/'+'dat_EF_{0}.csv'.format(case_name)      , sep=',', index=False)
    EF_years.to_csv(inter_dir+'/'+'EF_years_{0}.csv'.format(case_name) , sep=',', index=False)
    EF_names.to_csv(inter_dir+'/'+'EF_names_{0}.csv'.format(case_name) , sep=',', index=False)
    EF_fuels.to_csv(inter_dir+'/'+'EF_fuels_{0}.csv'.format(case_name) , sep=',', index=False)
    EF_polls.to_csv(inter_dir+'/'+'EF_polls_{0}.csv'.format(case_name) , sep=',', index=False)

    
    EmisFactor_yr_spd = EF_Grid_table(merge_EF_Final, EF_years, EF_names, EF_fuels, EF_polls)
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    print('')
    print('******************************************************************')
    print('*****              Deterioration_factor is done                *****')
    print('******************************************************************')
    print('')    
    
    return EmisFactor_yr_spd





# =============================================================================
# Function to calculate the county emissions 
# =============================================================================
def calc_County_Emissions(Input_dir, Emissions_Factor_DataFrame, 
                          Activity_Data_DataFrame, Roads_DataFrame, County_DataFrame,
                          Cold_Start_input, avgSpeedDist, Out_dir ):   
    start_time = time.time()
    print('')
    print('******************************************************************')
    print('*****           Starting calculation of emissions            *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')
    input_dir  = Input_dir
    output_dir = Out_dir                     #
    EF_yr_spd  = Emissions_Factor_DataFrame  #
    AD_yr_vhc  = Activity_Data_DataFrame     #
    roads_DF   = Roads_DataFrame             #
    county_df  = County_DataFrame            # 
    cold_start = Cold_Start_input

    # input_dir  = input_dir         
    # output_dir = output_dir
    # avgSpeedDist = avgSpeedDist
    # EF_yr_spd  = EmisFactor_yr_spd 
    # AD_yr_vhc  = AD_SK             
    # roads_DF   = roads_RGS         
    # county_df  = county_SHP        
    # cold_start = Cold_Start_list  
    
    if calendar.isleap(pd.to_datetime(STDATE).year) == True:
        ndaysYEAR = 366
    else:
        ndaysYEAR = 365
        
    temp_EF_df = EF_yr_spd.data.copy()
    nspd = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 89, 97, 105, 113, 121]
    EFyears = list(np.sort(temp_EF_df.years.unique()))
    ADyears = list(np.sort(AD_yr_vhc.data.manufacture_date.unique()))
    npol = list(EF_yr_spd.EF_polls.EF_poll)
    EF_AD_max_yr = max(np.max(ADyears), np.max(EFyears))
    EF_AD_min_yr = min(np.min(ADyears), np.min(EFyears))
    EF_yr_min = np.min(EFyears)
    nyears = list(np.arange(EF_AD_min_yr,EF_AD_max_yr+1))
    ADvhc = list(np.sort(AD_yr_vhc.fullname.vhc_name))
    EFvhc = list(np.sort(EF_yr_spd.EF_fullname.EF_fullname))
    missVHC = []
    for j in ADvhc:
        if j not in EFvhc:
            missVHC.append(j)
    missVHC = list(dict.fromkeys(missVHC))
    
    ADvhc_notin_EFvhc = 'vehicles_with_AD_without_EF.txt'
    with open(output_dir+'/LOGS/'+ADvhc_notin_EFvhc, 'w') as log:
        log.write('# ***** Vehicles with Activity data without Emissions Factors *****\n')
        log.write('# \n')
        log.write('vehicle_engine_fuel\n')
        for i in missVHC:
            log.write('{0}\n'.format(i))

    final_EF_df = pd.DataFrame()
    for ispd in nspd:
        for ipol in list(temp_EF_df.pollutant.unique()):
            yr_df = pd.DataFrame({'years' : nyears,
                                  'spd'   : [ispd for i in range(0,len(nyears))],
                                  'pollutant' : [ipol for i in range(0,len(nyears))]})
            
            aux_ef = temp_EF_df.loc[(temp_EF_df.pollutant == ipol) & (temp_EF_df.spd == ispd)]
            aux_ef = (aux_ef.pivot(index = 'years', columns='fullname', values='emis_fact')).reset_index(drop=False)
            
            aux_ef = pd.merge(yr_df,aux_ef, left_on = 'years', 
                              right_on = 'years', how='left').reset_index(drop=True)
            aux_ef = aux_ef.ffill(axis = 0).bfill(axis = 0) 

            final_EF_df = final_EF_df.append(aux_ef, ignore_index=True, sort=False).fillna(0.0)

    miss_EF_df = final_EF_df.loc[:,['years','spd','pollutant']]
    for mvhc in missVHC:
        auxdf = pd.DataFrame({mvhc :  (miss_EF_df.spd * 0.0).values})
        miss_EF_df = pd.concat([miss_EF_df, auxdf], axis=1, sort=False)

    final_EF_df = pd.merge(final_EF_df,miss_EF_df, left_on = ['years','spd','pollutant'], 
                              right_on = ['years','spd','pollutant'], how='left').reset_index(drop=True)
    
    EFvhc = list(final_EF_df.columns)
    EFvhc.remove('years')
    EFvhc.remove('spd')
    EFvhc.remove('pollutant')
    
    road_weight = roads_DF.roads_df.groupby(['region_cd', 'road_type']).sum().reset_index(drop=False)
    road_weight = road_weight.sort_values(by=['region_cd', 'road_type'])
    road_weight = road_weight.loc[:,['region_cd', 'road_type','vkt_split_county']]
    auxroad_weight = pd.DataFrame(np.ones((road_weight.shape[0],len(ADvhc))), columns=ADvhc)
    road_weight = pd.concat([road_weight, auxroad_weight], axis=1, sort=False)
    road_weight.loc[:,ADvhc] = road_weight.loc[:,ADvhc].apply(lambda x: x * road_weight.vkt_split_county.values)

    if (process_road_restriction.lower() == 'yes') or (process_road_restriction.lower() == 'y'):

        inroad = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',road_restriction)
        if os.path.exists(inroad) == True:
            print ('')
            print ('Reading Road restriction input file ...')
            print (inroad)
            
            print ('')
            print ('Processing Roads Restriction...')

    
        with open(inroad, 'r') as csvfile:
            sep = csv.Sniffer().sniff(csvfile.read(4096)).delimiter

        roadsRes = pd.read_csv(inroad, sep = sep).fillna(0)
        roadsRes.columns    = roadsRes.columns.str.lower()
        roadsRes.iloc[:,3:] = roadsRes.iloc[:,3:].astype(int)
        roadsRes.loc[:,'vehicle'] = roadsRes.loc[:,'vehicle'].str.lower()
        roadsRes.loc[:,'engine']   = roadsRes.loc[:,'engine'].str.lower()
        roadsRes.loc[:,'fuel']    = roadsRes.loc[:,'fuel'].str.lower()
        roadsRes['fullname']      = roadsRes.vehicle.str.cat(roadsRes[['engine','fuel']], sep='_')

        roadsRes = roadsRes.loc[~roadsRes.road1.isin([0.0])].reset_index(drop=True)
        roadsRes = roadsRes.drop(columns=['vehicle','engine','fuel'])
        roadsRes = roadsRes.astype(str)
        roadsRes['aggroads'] = roadsRes['road1'].str.cat(roadsRes[['road2','road3',
                                                                   'road4','road5',
                                                                   'road6','road7',
                                                                   'road8']],sep=",")

        roadsRes = roadsRes.groupby('aggroads')['fullname'].apply(list).reset_index(drop=False)
        
        def cleanlist(listOfnum):
            new = []
            for elem in listOfnum:
                if elem != '0':
                    new.append(elem)
            return new 
        road_weightEX = road_weight.loc[:,['region_cd','road_type']]
        road_EXaux = pd.DataFrame()
        for index, row in roadsRes.iterrows():
            prfl = row.aggroads
            vhcEX = row.fullname
            EXRoads = cleanlist(prfl.split(','))
            for irgn in road_weight.region_cd.unique():
                auxEX = road_weight.loc[road_weight.region_cd == irgn]
                auxEX['prfl'] = prfl
                shp   = auxEX.loc[~auxEX.road_type.isin(EXRoads)].shape[0]
                if shp > 0:
                    ival = auxEX.vkt_split_county.loc[auxEX.road_type.isin(EXRoads)].sum() / shp
                    auxEX.vkt_split_county = auxEX.vkt_split_county.add(ival)
                    auxEX.vkt_split_county.loc[auxEX.road_type.isin(EXRoads)] = 0
                    road_EXaux = road_EXaux.append(auxEX)
                else:
                    road_EXaux = road_EXaux.append(auxEX)
            prfl_val = road_EXaux.loc[road_EXaux.prfl == prfl, 'vkt_split_county']
            auxdfEX = pd.DataFrame(np.ones((prfl_val.shape[0],len(vhcEX))), columns=vhcEX)
            auxdfEX.loc[:,vhcEX] = auxdfEX.loc[:,vhcEX].apply(lambda x: x * prfl_val.values)
            auxdfEX = pd.concat([road_EXaux.loc[road_EXaux.prfl == prfl, ['region_cd','road_type']], auxdfEX], axis=1, sort=False)
            
            road_weightEX = pd.merge(road_weightEX, auxdfEX, 
                                     left_on=['region_cd','road_type'],
                                     right_on=['region_cd','road_type'])

        #changing the vkt_split by the vehicles with restrictions
        road_weight.loc[:,vhcEX] = road_weightEX.loc[:,vhcEX]
        road_weight = road_weight.drop(columns=['vkt_split_county'])

    final_AD_df = AD_yr_vhc.data.copy()
    inter_AD = final_AD_df.loc[final_AD_df.manufacture_date <= EF_yr_min]
    drop_list_AD = list(inter_AD.index)
    inter_AD = inter_AD.groupby(['region_cd']).sum().reset_index(drop=False)
    inter_AD.loc[:,'manufacture_date'] = EF_yr_min
    final_AD_df = final_AD_df.drop(index=drop_list_AD).reset_index(drop=True)
    final_AD_df = final_AD_df.append(inter_AD,ignore_index=True, sort=False)
    final_AD_df = final_AD_df.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)

    final_AD_df.loc[:,:] = final_AD_df.loc[:,:].astype(float)
    final_AD_df.loc[:,'region_cd'] = final_AD_df.loc[:,'region_cd'].astype(int)
    final_AD_df.loc[:,'manufacture_date'] = final_AD_df.loc[:,'manufacture_date'].astype(int)

    AD_byRoad = pd.merge(final_AD_df, road_weight.loc[:,['region_cd', 'road_type']], 
                         left_on = ['region_cd'], 
                         right_on = ['region_cd'], how='left').fillna(0.0)
    
    road_weight = pd.merge(final_AD_df.loc[:,['manufacture_date', 'region_cd']],
                           road_weight, left_on = ['region_cd'], 
                           right_on = ['region_cd'], how='left').fillna(0.0)


    AD_byRoad.loc[:,ADvhc] = AD_byRoad.loc[:,ADvhc] * road_weight.loc[:,ADvhc]
    ISzeroRD = AD_byRoad['road_type'].loc[AD_byRoad['road_type'] == 0].index 
    AD_byRoad['road_type'].iloc[ISzeroRD] = 105
    AD_byRoad['road_type'] = AD_byRoad['road_type'].astype(int)
    
    road_by_county = AD_byRoad.groupby(['region_cd','road_type']).sum().reset_index(drop=False)
    road_by_county = road_by_county.drop(columns=['manufacture_date'])
    road_by_county = road_by_county.sort_values(by=['region_cd','road_type'])
    aux_RbC        = road_by_county.groupby(['region_cd']).sum().reset_index(drop=False)
    aux_RbC        = aux_RbC.drop(columns=['road_type'])
    aux_RbC        = pd.merge(aux_RbC, road_by_county.loc[:,['region_cd','road_type']],
                              left_on = 'region_cd', right_on='region_cd', how='left')
    road_by_county = road_by_county.set_index(['region_cd','road_type'])
    aux_RbC = aux_RbC.set_index(['region_cd','road_type'])
    road_by_county.loc[:,ADvhc] = road_by_county.loc[:,ADvhc].div(aux_RbC.loc[:,ADvhc], axis=0)
    road_by_county = road_by_county.reset_index(drop=False).fillna(1.0)
    road_by_county.to_csv(inter_dir+'/'+'road_by_district_{0}.csv'.format(case_name))
    EF_avgSpd = pd.DataFrame()
    roads = ['101','102','103','104','105','106','107','108']
    for iroad in roads:
        aux_avgSpd = avgSpeedDist.data.loc[:,['speed',str(iroad)]]
        
        aux_EFavgSpd = pd.merge(final_EF_df,aux_avgSpd, left_on = 'spd', 
                                  right_on = 'speed', how='right').reset_index(drop=True)
        
        aux_EFavgSpd.loc[:,EFvhc] = aux_EFavgSpd.loc[:,EFvhc].apply(lambda x: np.asarray(x) * aux_EFavgSpd[iroad].values)
        aux_EFavgSpd = aux_EFavgSpd.groupby(['years','pollutant']).sum().reset_index(drop=False)
        aux_EFavgSpd = aux_EFavgSpd.drop(columns=['spd', 'speed', iroad])
        aux_EFavgSpd['road_type'] = aux_EFavgSpd['years'] * 0 + int(iroad)
        EF_avgSpd = EF_avgSpd.append(aux_EFavgSpd, ignore_index=True, sort=False)

    Emissions_by_yr_county = pd.DataFrame()
    for ipol in npol:
        EF_by_AD_df = pd.merge(AD_byRoad.loc[:,['manufacture_date','region_cd','road_type']], 
                       EF_avgSpd.loc[EF_avgSpd['pollutant'] == ipol,['years','pollutant','road_type']+ADvhc],
                       left_on = ['manufacture_date','road_type'],
                       right_on = ['years','road_type'], how='left').reset_index(drop=True)
                       

        EF_by_AD_df.loc[:,ADvhc] = EF_by_AD_df.loc[:,ADvhc].mul(AD_byRoad.loc[:,ADvhc], axis=1)
        EF_by_AD_df = EF_by_AD_df.groupby(['manufacture_date','region_cd','pollutant','road_type']).sum().reset_index(drop=False)
        EF_by_AD_df = EF_by_AD_df.drop(columns=['years']).sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
        
        Emissions_by_yr_county = Emissions_by_yr_county.append(EF_by_AD_df, ignore_index=True, sort=False)
    
    Emissions_by_yr_county = Emissions_by_yr_county.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
    
    Emissions_by_yr_county2csv = Emissions_by_yr_county.groupby(['region_cd','pollutant']).sum().reset_index(drop=False)
    Emissions_by_yr_county2csv.loc[:,ADvhc] = Emissions_by_yr_county2csv.loc[:,ADvhc] * (ndaysYEAR/1000000)
    Emissions_by_yr_county2csv = Emissions_by_yr_county2csv.drop(columns=['manufacture_date','road_type'])
    Emissions_by_yr_county2csv.to_csv(output_dir+'/hot_exhaust_emissions_Tons_per_Year.csv', sep=',', index=False)
    Emissions_by_yr_county2csv = []
    '''
    Calculating the Cold Start Emissions
    beta = 0.647 - (0.025 * 1Trip) - (0.00974 - 0.000385 * 1Trip)* temp_mean
    1Trip = 12.35km - This value came from South Korea survey and it means that
    each vehicle run, on average, 12.35 km per day
    '''
        
    for ifile in range(0,len(cold_start)):
#        print(ef_list[ifile])
        name = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',cold_start[ifile])

        if os.path.exists(name) == True:
            print('******************************************************************')
            print('*****   Reading Cold Start vehicles and roads tables ...     *****')
            print('*****                    Please wait ...                     *****')
            print('******************************************************************')
            print('')
            print (name)
            with open(name, 'r') as csvfile:
                sep = csv.Sniffer().sniff(csvfile.read(4096)).delimiter

            cold_start_df = (pd.read_csv(name, sep = sep)).fillna(0)
            cold_start_df.columns = cold_start_df.columns.str.lower()
            cold_start_df.loc[:,'vehicle'] = cold_start_df.loc[:,'vehicle'].str.lower()
            cold_start_df.loc[:,'engine'] = cold_start_df.loc[:,'engine'].str.lower()
            cold_start_df.loc[:,'fuel'] = cold_start_df.loc[:,'fuel'].str.lower()
            cold_start_df['fullname'] = cold_start_df.vehicle.str.cat(cold_start_df[['engine','fuel']], sep='_')
        else:
            print('')
            print('***** ERROR *****: There is no Cold Start csv file')
            print('***** ERROR *****: Please check the Cold start input file')
            sys.exit('')


    ColdStart_roads = list(cold_start_df.road_type.unique())
    ColdStart_fuels = ['gasoline','diesel','lpg']
    aux_list_CS = list(cold_start_df['fullname'].unique())
    CS_vhc_list = list(cold_start_df['fullname'].unique())
    EFvhc = EF_yr_spd.EF_fullname.EF_fullname.tolist()
    ADvhc = AD_yr_vhc.fullname.vhc_name.tolist()
    CSlog_name = 'LOG_Cold_Start_Vehicle_checking.txt'
    miss_CS = []
    with open(output_dir+'/LOGS/'+CSlog_name, 'w') as CS_log:
        CS_log.write('*****   List of vehicle without Cold Start Emissions    *****\n')
        for ivhc in aux_list_CS:
            if (ivhc not in ADvhc) or (ivhc not in EFvhc):
                miss_CS.append(ivhc)
                CS_vhc_list.remove(ivhc)
                CS_log.write('**** WARNING ***:There is no EF for vehicle {0} or the vehicle do not exist\n'.format(ivhc))
    
    if len(CS_vhc_list) == 0:
        print('')
        print('***** ERROR *****: There is no vehicle in the Cold Start list which match with AD or EF tables')
        print('***** ERROR *****: Please check the Cold start input list')
        sys.exit('')
    
    cold_start_aux = Emissions_by_yr_county.loc[(Emissions_by_yr_county.road_type.isin(ColdStart_roads)) &
                                                Emissions_by_yr_county.pollutant.isin(['CO','VOC','NOX','PM2.5', 'PM10'])].reset_index(drop=True)
    
    cold_start_aux = cold_start_aux.loc[:,['manufacture_date','region_cd','road_type','pollutant']+CS_vhc_list]

    beta = 0.647 - (0.025 * 12.35) - (0.00974 - 0.000385 * 12.35) * temp_mean
    
    EcEh_Gasoline_CO  =  beta * ((9.040 - 0.090 * temp_mean) - 1)
    EcEh_Gasoline_VOC =  beta * ((12.59 - 0.060 * temp_mean) - 1)
    EcEh_Gasoline_NOX =  beta * ((3.660 - 0.006 * temp_mean) - 1)

    EcEh_Diesel_CO    =  beta * ((1.900 - 0.030 * temp_mean) - 1)
    if temp_mean > 29:
        EcEh_Diesel_VOC   =  beta * 0.5
    else:
        EcEh_Diesel_VOC   =  beta * ((3.100 - 0.090 * temp_mean) - 1)

    EcEh_Diesel_NOX   =  beta * ((1.300 - 0.013 * temp_mean) - 1)

    if temp_mean > 26:
        EcEh_Diesel_PM    =  beta * 0.5
    else:
        EcEh_Diesel_PM    =  beta * ((3.100 - 0.100 * temp_mean) - 1)


    EcEh_LPG_CO       =  beta * ((3.660 - 0.090 * temp_mean) - 1)

    if temp_mean > 29:
        EcEh_LPG_VOC      =  beta * 0.5
    else:
        EcEh_LPG_VOC      =  beta * ((2.240 - 0.060 * temp_mean) - 1)
        
    EcEh_LPG_NOX      =  beta * ((0.980 - 0.006 * temp_mean) - 1)
    
    aux_dict = {'pollutant' : ['CO','VOC','NOX','PM2.5', 'PM10'],
                'gasoline'  : [EcEh_Gasoline_CO, EcEh_Gasoline_VOC, EcEh_Gasoline_NOX, 0.0,            0.0],
                'diesel'    : [EcEh_Diesel_CO  , EcEh_Diesel_VOC  , EcEh_Diesel_NOX  , EcEh_Diesel_PM, EcEh_Diesel_PM ],
                'lpg'       : [EcEh_LPG_CO     , EcEh_LPG_VOC     , EcEh_LPG_NOX     , 0.0,            0.0]}    
    EcEhBeta_df = pd.DataFrame(aux_dict)
    EcEhBeta_df.loc[:,ColdStart_fuels] = EcEhBeta_df.loc[:,ColdStart_fuels].clip(lower=0)

    cold_start_emissions = pd.DataFrame()
    for ifuel in ColdStart_fuels:
        aux_CS_emis = pd.DataFrame()
        for ipol in ['CO','VOC','NOX','PM2.5', 'PM10']:
            aux_CS_vhc = list(cold_start_aux.loc[cold_start_aux.pollutant.isin([ipol]),
                               cold_start_aux.columns.str.contains(ifuel)].columns)
            aux_CS_cols = ['manufacture_date', 'region_cd', 'road_type', 'pollutant']

            aux_fuel = cold_start_aux.loc[cold_start_aux.pollutant.isin([ipol]),
                                           aux_CS_cols + aux_CS_vhc]
            
            cs_ef_aux = np.float(EcEhBeta_df.loc[EcEhBeta_df.pollutant == ipol, ifuel].to_numpy())
            aux_fuel.loc[:,aux_CS_vhc] = aux_fuel.loc[:,aux_CS_vhc] * cs_ef_aux
            
            aux_CS_emis = aux_CS_emis.append(aux_fuel)
        
        if cold_start_emissions.shape[0] == 0:
            cold_start_emissions = cold_start_emissions.append(aux_CS_emis, ignore_index=True)
        else:
            cold_start_emissions = pd.merge(cold_start_emissions,aux_CS_emis, 
                                            left_on = ['manufacture_date', 'region_cd', 'road_type', 'pollutant'],
                                            right_on = ['manufacture_date', 'region_cd', 'road_type', 'pollutant'],
                                            how='left')

    cols_2merge = ['manufacture_date', 'region_cd', 'road_type', 'pollutant']
    merge_CS_HE = pd.merge(Emissions_by_yr_county.loc[:,cols_2merge],cold_start_emissions,
                   left_on = cols_2merge,
                   right_on = cols_2merge,
                   how='left').fillna(0.0)
    
    #Adding Cold Start Emissions into the emissions DF
    Emissions_by_yr_county.loc[:,CS_vhc_list] = Emissions_by_yr_county.loc[:,CS_vhc_list].add \
                                                (merge_CS_HE.loc[:,CS_vhc_list], axis=1)
                                                
    #converting Cold Start emissions to Tons/yr to output it
    cold_start_emissions.loc[:,CS_vhc_list] = cold_start_emissions.loc[:,CS_vhc_list] * (ndaysYEAR/1000000)
    cold_start_emissions = cold_start_emissions.groupby(['region_cd','pollutant']).sum().reset_index(drop=False)
    cold_start_emissions = cold_start_emissions.drop(columns=['road_type','manufacture_date'])
    cold_start_emissions.to_csv(output_dir+'/cold_start_emissions_Tons_per_Year.csv', sep=',', index=False)
    
    print('')
    print('******************************************************************')
    print('*****      Calculating the Evaporative Emissions ...        *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    '''
    Calculating the Evaporative Emissions
    beta = 0.647 - (0.025 * 1Trip) - (0.00974 - 0.000385 * 1Trip)* temp_mean
    1Trip = 12.35km - This value came from South Korea survey and it means that
    each vehicle run, on average, 12.35 km per day
    '''   
    RVP    = 82.0
    beta   = 0.647 - (0.025 * 12.35) - (0.00974 - 0.000385 * 12.35) * temp_mean
    #The new emissions factors were updated because all vehicles are now considered
    # as Carbon canister controlled vehicle
    ErHot  = 0.1 * 0.136 * (np.exp(-5.967 + 0.04259 * RVP + 0.1773 * temp_mean ))
    ErWarm = 0.1 * 0.100 * (np.exp(-5.967 + 0.04259 * RVP + 0.1773 * temp_mean ))
    Ed_aux = 0.2 * 9.100 * (np.exp(0.0158*(RVP-61.2) + 0.0574 *(temp_min - 22.5) + 0.0614 * \
                           (temp_max - temp_min - 11.7)))

    aux_EE_AD = AD_byRoad.copy()
    list_vhc_no_gas = list(aux_EE_AD.loc[0,~aux_EE_AD.columns.str.contains('gasoline')].index)
    list_vhc_no_gas.remove('manufacture_date')
    list_vhc_no_gas.remove('region_cd')
    list_vhc_no_gas.remove('road_type')
    # list_vhc_no_gas.remove('vkt_split_county')
    aux_EE_AD.loc[:,list_vhc_no_gas] = aux_EE_AD.loc[:,list_vhc_no_gas] * 0.0
    
    # Ed emissions is based on the factor Ed_aux, so it is just multiply the vhc count by the Ed_aux
    Ed = AD_yr_vhc.vhc_count.copy()
    Ed.loc[:,list_vhc_no_gas] = Ed.loc[:,list_vhc_no_gas] * 0.0
    Ed.loc[:,ADvhc] = Ed.loc[:,ADvhc].mul(Ed_aux)
    Ed['region_cd'] = Ed['region_cd'].astype(int)
    Ed = pd.merge(aux_EE_AD.loc[:,['manufacture_date','region_cd','road_type']],
                  Ed,
                  left_on  = ['manufacture_date','region_cd'],
                  right_on = ['manufacture_date','region_cd'], how='left')
    
    Ed.loc[:,ADvhc] = (Ed.loc[:,ADvhc] * road_weight.loc[:,ADvhc]).fillna(0.0)

    #Calculating R based on South Korea equation
    R_aux = ((1 - beta) * ErHot) + (beta * ErWarm)
    R     = aux_EE_AD.copy()
    R.loc[:,ADvhc] = R.loc[:,ADvhc].mul(R_aux)
    R.set_index(['manufacture_date','region_cd','road_type'], inplace=True)
    R = R.sort_index()
    R = R.reset_index(drop=False)

    # Applying conversion factor to reflect engine size of motocycle(0.075)
    Ed.loc[:,Ed.columns.str.contains('motorcycle')] = Ed.loc[:,Ed.columns.str.contains('motorcycle')] * 0.075
    R.loc[:,R.columns.str.contains('motorcycle')]   = R.loc[:,R.columns.str.contains('motorcycle')]   * 0.075

    #Calculating Sfi based on South Korea equation
    # Sfi = aux_EE_AD.copy()
    # Sfi.loc[:,ADvhc] = Sfi.loc[:,ADvhc].mul(0.7 / 12.35)
    # Sfi.set_index(['manufacture_date','region_cd','road_type'], inplace=True)
    # Sfi = Sfi.sort_index()
    # Sfi = Sfi.reset_index(drop=False)

    Ed  = Ed.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
    R   = R.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
    # Sfi = Sfi.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
    

    #Copying Ed emissions to the final dataframe
    evaporative_emissions = Ed.copy()
    #Adding R emissions to final dataframe
    evaporative_emissions.loc[:,ADvhc] = evaporative_emissions.loc[:,ADvhc].add \
                                                (R.loc[:,ADvhc], axis=1)
    # #Adding Sfi emissions to final dataframe
    # evaporative_emissions.loc[:,ADvhc] = evaporative_emissions.loc[:,ADvhc].add \
    #                                             (Sfi.loc[:,ADvhc], axis=1)
        
    evaporative_emissions['pollutant'] = 'VOC'
    evaporative_emissions = evaporative_emissions.loc[:,
    ['manufacture_date','region_cd','pollutant','road_type'] + ADvhc]

    cols_2merge = ['manufacture_date', 'region_cd', 'road_type', 'pollutant']
    merge_EE_HE = pd.merge(Emissions_by_yr_county.loc[:,cols_2merge],evaporative_emissions,
                   left_on = cols_2merge,
                   right_on = cols_2merge,
                   how='left').fillna(0.0)
    
    Emissions_by_yr_county.loc[:,ADvhc] = Emissions_by_yr_county.loc[:,ADvhc].add \
                                                (merge_EE_HE.loc[:,ADvhc], axis=1)
   
    #converting Evaporative emissions to Tons/yr to output it
    evaporative_emissions.loc[:,ADvhc] = evaporative_emissions.loc[:,ADvhc] * (ndaysYEAR/1000000)
    evaporative_emissions = evaporative_emissions.groupby(['region_cd','pollutant']).sum().reset_index(drop=False)
    evaporative_emissions.to_csv(output_dir+'/evaporative_emissions_Tons_per_Year.csv', sep=',', index=False)


    total_county = Emissions_by_yr_county.groupby(['region_cd','pollutant']).sum().reset_index(drop=False)
    total_county = total_county.drop(columns=['manufacture_date','road_type'])
    
    Emissions_by_road = Emissions_by_yr_county.groupby(['pollutant','road_type']).sum().reset_index(drop=False)
    Emissions_by_road = Emissions_by_road.drop(columns=['manufacture_date','region_cd'])
    Emissions_by_road.loc[:,ADvhc] = Emissions_by_road.loc[:,ADvhc] * (ndaysYEAR/1000000)
    Emissions_by_road['total_emis'] = Emissions_by_road.loc[:,ADvhc].sum(axis=1)
    Emissions_by_road.to_csv(output_dir+'/'+'Pollutant_Total_Emis_by_Road_{0}.csv'.format(case_name), sep=',', index=False)
    total_poll = []


    #Calculation the total emissions for county and link level
    Emissions_by_yr_county = Emissions_by_yr_county.groupby(['region_cd','pollutant','manufacture_date']).sum().reset_index(drop=False)
    Emissions_by_yr_county['total_emis'] = Emissions_by_yr_county.loc[:,ADvhc].sum(axis=1)
    total_county['total_emis'] = total_county.loc[:,ADvhc].sum(axis=1)
    county_table = pd.DataFrame({'region_cd' : Emissions_by_yr_county.region_cd.unique()})
    county_table.to_csv(inter_dir+'/'+'district_table_{0}.csv'.format(case_name))
    
    years_table  = pd.DataFrame({'years' : Emissions_by_yr_county.manufacture_date.unique()})
    years_table.to_csv(inter_dir+'/'+'years_table_{0}.csv'.format(case_name))
    
    vhc_table    = pd.DataFrame({'fullname' : ADvhc})
    vhc_table.to_csv(inter_dir+'/'+'vhc_table_{0}.csv'.format(case_name))
    
    #Inserting the Geometry in the county emissions
    total_county = pd.merge(total_county,county_df, left_on = ['region_cd'], 
                        right_on = ['region_cd'], how='left')
    nanGeometry = list(total_county.loc[total_county['geometry'].isna()].index)
    total_county = total_county.drop(index=nanGeometry).reset_index(drop=True)
    total_county = total_county.loc[:,['region_cd','pollutant','geometry']+ADvhc]
    #adding georeference in the emissions
    total_county = gpd.GeoDataFrame(total_county, crs=county_df.crs, geometry='geometry')
    
    #Calculating th weight factor to apply the normalization by VKT
    aux_df_wgt = roads_DF.roads_df.groupby(['region_cd']).sum()
    aux_df_wgt = aux_df_wgt.reset_index(drop=False)
    aux_df_wgt = pd.merge(county_table, aux_df_wgt.loc[:,['region_cd','vkt_avg']], 
                          left_on=['region_cd'], right_on=['region_cd'], how='left')
    
    emislog_name = 'LOG_Calculation_emissions_by_County.txt'
    with open(output_dir+'/LOGS/'+emislog_name, 'w') as EF_log:
        EF_log.write('*****   WARNING   *****')
        for icd in aux_df_wgt.region_cd.loc[aux_df_wgt.vkt_avg.isna()]:
            EF_log.write('There is no VKT data for the county {0}! Please, review your input link shapefile\n'.format(icd))
    aux_df_wgt = aux_df_wgt.fillna(0)
    aux_df_wgt.loc[:,'region_cd'] = aux_df_wgt.loc[:,'region_cd'].astype(str)
    aux_df_wgt['region_state']    = aux_df_wgt.loc[:,'region_cd'].str.slice(stop=2).astype(int)
    aux_df_wgt.loc[:,'region_cd'] = aux_df_wgt.loc[:,'region_cd'].astype(int)

    rgn_st = aux_df_wgt.groupby(['region_state']).sum().reset_index(drop=False)
    rgn_st = rgn_st.rename(columns={'vkt_avg' : 'vkt_avg_sum'})
    aux_df_wgt = pd.merge(aux_df_wgt, rgn_st.loc[:,['region_state','vkt_avg_sum']],
                          left_on='region_state', right_on='region_state', how='left')
    aux_df_wgt['vkt_weight'] = aux_df_wgt.loc[:,'vkt_avg'] /  aux_df_wgt.loc[:,'vkt_avg_sum']
    aux_df_wgt = aux_df_wgt.sort_values(by=['region_cd']).reset_index(drop=True).fillna(0.0)
    
    total_county.loc[:,'region_cd'] = total_county.loc[:,'region_cd'].astype(str)
    total_county['region_state']    = total_county.loc[:,'region_cd'].str.slice(stop=2).astype(int)
    total_county.loc[:,'region_cd'] = total_county.loc[:,'region_cd'].astype(int)

    #Apply the weight factor to normalize the emissions by county VKT
    total_county_WGT = total_county.groupby(['region_state','pollutant']).sum().reset_index(drop=False)
    total_county_WGT = total_county_WGT.drop(columns=['region_cd'])
    total_county_WGT = pd.merge(total_county_WGT, aux_df_wgt.loc[:,['region_cd','region_state','vkt_weight']], 
                                left_on=['region_state'], right_on=['region_state'], how='left')
    total_county_WGT.loc[:,ADvhc] = total_county_WGT.loc[:,ADvhc].apply(lambda x: np.asarray(x) * total_county_WGT.vkt_weight.values)
    total_county_WGT = total_county_WGT.loc[:,['region_cd','pollutant', 'vkt_weight']+ADvhc]
    total_county_WGT = pd.merge(total_county_WGT, total_county.loc[:,['region_cd', 'pollutant', 'geometry']],
                                left_on  = ['region_cd', 'pollutant'],
                                right_on = ['region_cd', 'pollutant'], how='left')
    total_county_WGT = gpd.GeoDataFrame(total_county_WGT, crs=county_df.crs, geometry='geometry')

    #Calculationg the link emissions based on Normalized County total
    Emissions_by_link = pd.merge(roads_DF.roads_df,total_county_WGT, left_on = ['region_cd'], 
                        right_on = ['region_cd'], how='left')
    Emissions_by_link.loc[:,ADvhc] = Emissions_by_link.loc[:,ADvhc].apply(lambda x: np.asarray(x) * \
                         Emissions_by_link.vkt_split_county.values)    

    #Apply the weight factor to normalize the emissions by county VKT and year
    Emissions_by_yr_county.loc[:,'region_cd'] = Emissions_by_yr_county.loc[:,'region_cd'].astype(str)
    Emissions_by_yr_county['region_state']    = Emissions_by_yr_county.loc[:,'region_cd'].str.slice(stop=2).astype(int)
    Emissions_by_yr_county.loc[:,'region_cd'] = Emissions_by_yr_county.loc[:,'region_cd'].astype(int)

    Emissions_by_yr_county_WGT = Emissions_by_yr_county.groupby(['region_state','manufacture_date','pollutant']).sum().reset_index(drop=False)
    Emissions_by_yr_county_WGT = Emissions_by_yr_county_WGT.drop(columns=['region_cd', 'road_type'])
    Emissions_by_yr_county_WGT = pd.merge(Emissions_by_yr_county_WGT, aux_df_wgt.loc[:,['region_cd','region_state','vkt_weight']], 
                                left_on=['region_state'], right_on=['region_state'], how='left')

    Emissions_by_yr_county_WGT.loc[:,ADvhc] = Emissions_by_yr_county_WGT.loc[:,ADvhc].apply(lambda x: np.asarray(x) * Emissions_by_yr_county_WGT.vkt_weight.values)
    Emissions_by_yr_county_WGT = Emissions_by_yr_county_WGT.loc[:,['region_cd','manufacture_date','pollutant']+ADvhc+['total_emis']]
    
    #Calculating the annual total emissions
    total_county.loc[:,ADvhc]               = total_county.loc[:,ADvhc] * (ndaysYEAR/1000000) # tons/year
    total_county_WGT.loc[:,ADvhc]           = total_county_WGT.loc[:,ADvhc] * (ndaysYEAR/1000000) # tons/year
    Emissions_by_yr_county.loc[:,ADvhc]     = Emissions_by_yr_county.loc[:,ADvhc] * (ndaysYEAR/1000000) # tons/year
    Emissions_by_yr_county_WGT.loc[:,ADvhc] = Emissions_by_yr_county_WGT.loc[:,ADvhc] * (ndaysYEAR/1000000) # tons/year

    total_county['total_emis']               = total_county.loc[:,ADvhc].sum(axis=1)
    total_county_WGT['total_emis']           = total_county_WGT.loc[:,ADvhc].sum(axis=1)
    Emissions_by_yr_county['total_emis']     = Emissions_by_yr_county.loc[:,ADvhc].sum(axis=1)
    Emissions_by_yr_county_WGT['total_emis'] = Emissions_by_yr_county_WGT.loc[:,ADvhc].sum(axis=1)

    print('')
    print('******************************************************************')
    print('*****      Saving the results in the .csv files   ...        *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')
    
    #Exporting the CSV files with the emissions 
    Emis_by_yr_county_2CVS = Emissions_by_yr_county.groupby(['region_cd', 'pollutant', 'manufacture_date']).sum().reset_index(drop=False)
    Emis_by_yr_county_2CVS = Emis_by_yr_county_2CVS.drop(columns=['road_type','region_state'])
    Emis_by_yr_county_2CVS.to_csv(output_dir+'/'+'District_Total_Emissions_by_Year_{0}.csv'.format(case_name), sep=',', index=False)
    Emis_by_yr_county_2CVS = []
    
    Emis_by_yr_county_WGT_2CVS = Emissions_by_yr_county_WGT.groupby(['region_cd', 'pollutant', 'manufacture_date']).sum().reset_index(drop=False)
#    Emis_by_yr_county_WGT_2CVS = Emis_by_yr_county_WGT_2CVS.drop(columns=['road_type','region_state'])
    Emis_by_yr_county_WGT_2CVS.to_csv(output_dir+'/'+'District_Total_Normalized_Emissions_by_Year_{0}.csv'.format(case_name), sep=',', index=False)
    Emis_by_yr_county_WGT_2CVS = []
    
    total_poll = total_county.groupby(['pollutant']).sum().reset_index(drop=False)
    total_poll = total_poll.drop(columns=['region_cd', 'region_state'])
    total_poll.to_csv(output_dir+'/'+'Pollutant_Total_Emissions_Tons_per_Year_{0}.csv'.format(case_name), sep=',', index=False)
    total_poll = []
    
    total_poll_WGT = total_county_WGT.groupby(['pollutant']).sum().reset_index(drop=False)
    total_poll_WGT = total_poll_WGT.drop(columns=['region_cd', 'vkt_weight'])
    total_poll_WGT.to_csv(output_dir+'/'+'Pollutant_Total_Normalized_Emissions_Tons_per_Year_{0}.csv'.format(case_name), sep=',', index=False)
    total_poll_WGT = []
    
    total_county_2CVS = total_county.drop(columns = ['geometry','region_state'])
    total_county_2CVS.to_csv(output_dir+'/'+'District_Total_Emissions_Tons_per_Year_{0}.csv'.format(case_name), sep=',', index=False)
    total_county_2CVS = []
    
    total_county_WGT_2CVS = total_county_WGT.drop(columns = ['geometry'])
    total_county_WGT_2CVS.to_csv(output_dir+'/'+'District_Total_Normalized_Emissions_Tons_per_Year_{0}.csv'.format(case_name), sep=',', index=False)
    total_county_WGT_2CVS = []

    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))
    print('')
    print('******************************************************************')
    print('*****             Emissions calculation is done              *****')
    print('******************************************************************')
    print('')
       
    return Emissions_table(Emissions_by_yr_county_WGT, total_county, total_county_WGT, 
                           county_table, years_table, vhc_table, road_by_county)





def aplly_control(Control_list_file, County_Emissions_df, County_SHP_df):
    start_time = time.time()
    print('')
    print('******************************************************************')
    print('*****        Applying control factor in the emissions        *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    control_list = Control_list_file    #Control_list      #
    cnt_df       = County_Emissions_df  # County_Emissions  #
    county_df    = County_SHP_df        #county_SHP        #
    ADvhc        = cnt_df.fullname.fullname.to_list()
    for ifile in range(0,len(control_list)):
#        print(ef_list[ifile])
        name = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',control_list[ifile])

        if os.path.exists(name) == True:
            print('******************************************************************')
            print('*****   Reading control factor for the vehicles ...          *****')
            print('*****                    Please wait ...                     *****')
            print('******************************************************************')
            print('')
            print (name)
            with open(name, 'r') as csvfile:
                sep = csv.Sniffer().sniff(csvfile.read(4096)).delimiter

            control_df = (pd.read_csv(name, sep = sep)).fillna(0)
            control_df.columns = control_df.columns.str.lower()
            control_df.loc[:,'vehicle'] = control_df.loc[:,'vehicle'].str.lower()
            control_df.loc[:,'engine'] = control_df.loc[:,'engine'].str.lower()
            control_df.loc[:,'fuel'] = control_df.loc[:,'fuel'].str.lower()
            control_df['fullname'] = control_df.vehicle.str.cat(control_df[['engine','fuel']], sep='_')
        else:
            print('')
            print('***** ERROR: There is no control factor file                 *****')
            print('***** Please, check the input file names                     *****')
            sys.exit()

    # control_vhc = list(control_df.fullname.unique())
    # control_region = list(control_df.region_cd.unique())
    Controllog_name = 'LOG_Control_factor_checking.txt'
    with open(output_dir+'/LOGS/'+Controllog_name, 'w') as CF_log:
        CF_log.write('*****   List of vehicle with more than one control factor per year    *****\n')
        
        duplicate = control_df.drop(columns = 'control_factor')
        aux_dup = duplicate.loc[control_df.duplicated()]
        if aux_dup.shape[0] != 0:
            print('')
            print('***** WARNING *****: Please check the LOG_Control_factor_checking.txt file')
            print('***** WARNING *****: There are vehicles with more than one control factor per year')
            print('')
            for i in aux_dup.index:
                CF_log.write('Vehicle {0} and region {1} has more than one control factor per year and/or per region.\n'.format(duplicate.loc[i,'fullname'], duplicate.loc[i,'region_cd']))
        
        duplicate.loc[:,'region_cd'] = duplicate.loc[:,'region_cd'].astype(str)
        duplicate.loc[:,'region_cd'] = duplicate.loc[:,'region_cd'].str.slice(stop=2).astype(int)
        duplicate.loc[:,'region_cd'] = duplicate.loc[:,'region_cd'].astype(int)
        aux_dup = duplicate.loc[duplicate.duplicated()]
        if aux_dup.shape[0] != 0:
            print('')
            print('***** WARNING *****: Please check the LOG_Control_factor_checking.txt file')
            print('***** WARNING *****: There are vehicles with more than one control factor per year')
            print('')
            for i in aux_dup.index:
                CF_log.write('Vehicle {0} and region {1} has more than one control factor per year and/or per region.\n'.format(duplicate.loc[i,'fullname'], duplicate.loc[i,'region_cd']))


    control_df.loc[:,'control_factor'] = 1 - (control_df.loc[:,'control_factor'] / 100)
    control_df.loc[:,'region_cd']      = control_df.loc[:,'region_cd'].astype(str)
    control_df['region_state']         = control_df.loc[:,'region_cd'].str.slice(stop=2).astype(int)
    control_df.loc[:,'region_cd']      = control_df.loc[:,'region_cd'].astype(int)
    
    cnt_df.county_by_yr.loc[:,'region_cd']      = cnt_df.county_by_yr.loc[:,'region_cd'].astype(str)
    cnt_df.county_by_yr['region_state']         = cnt_df.county_by_yr.loc[:,'region_cd'].str.slice(stop=2).astype(int)
    cnt_df.county_by_yr.loc[:,'region_cd']      = cnt_df.county_by_yr.loc[:,'region_cd'].astype(int)
    
    drop_all_pol_list = control_df.loc[control_df.data == 'ALL'].index
    control_df_all = control_df.loc[control_df.data == 'ALL'].reset_index(drop=True)
    control_df_all = control_df_all.drop(columns=['vehicle','engine','fuel'])
    control_df = control_df.drop(index=drop_all_pol_list).reset_index(drop=True)
    
    apply_factor_sta = cnt_df.county_by_yr.loc[:,['region_cd','region_state','pollutant', 'manufacture_date']].reset_index(drop=True)
    
    sta_vhc = []
    cnt_vhc = []
    
    control_df.loc[:,'region_cd']      = control_df.loc[:,'region_cd'].astype(int)
    if len(control_df.region_cd.unique()) > 0:
        CF_pol_sta = pd.DataFrame()
        CF_pol_cnt = pd.DataFrame()
        for irgn in control_df.region_cd.unique():
            aux_irg = control_df.loc[(control_df.region_cd == irgn)]
            if (len(str(irgn)) == 2) or (len(str(irgn)) == 1):
                sta_vhc.extend(list(aux_irg.fullname.unique()))
                for ipol in aux_irg.data.unique():
                    aux_CF = aux_irg.loc[(aux_irg.data == ipol)]
                    aux_CF = aux_CF.pivot(index='year', columns='fullname', values='control_factor').fillna(1).reset_index(drop=False)
                    aux_CF['data']      = ipol
                    aux_CF['region_cd'] = irgn
                    aux_CF = aux_CF.rename(columns={'region_cd' : 'region_general'})
                    CF_pol_sta = CF_pol_sta.append(aux_CF, ignore_index=True, sort=False)

            elif (len(str(irgn)) == 8) or (len(str(irgn)) == 7):
                cnt_vhc.extend(list(aux_irg.fullname.unique()))
                for ipol in aux_irg.data.unique():
                    aux_CF = aux_irg.loc[(aux_irg.data == ipol)]
                    aux_CF = aux_CF.pivot(index='year', columns='fullname', values='control_factor').fillna(1).reset_index(drop=False)
                    aux_CF['data']      = ipol
                    aux_CF['region_cd'] = irgn
                    aux_CF = aux_CF.rename(columns={'region_cd' : 'region_general'})
                    CF_pol_cnt = CF_pol_cnt.append(aux_CF, ignore_index=True, sort=False)
            else:
                print('')
                print('*****ERROR: Control factor is only applied in State level (2 digits state code e.g: 11) or')
                print('county level (8 digitis county code e.g: 11110129)')
                sys.exit()
    else:
        CF_pol_sta = None
        CF_pol_cnt = None
        print('')
        print('*** WARNING: There is no control factor for county level')
        print('')
    

    if len(control_df_all.region_cd.unique()) > 0:
        CF_all_sta = pd.DataFrame()
        CF_all_cnt = pd.DataFrame()
        for irgn in control_df_all.region_cd.unique():
            if (len(str(irgn)) == 2) or (len(str(irgn)) == 1):
                aux_all_cf = control_df_all.loc[control_df_all.region_cd == irgn]
                sta_vhc.extend(list(aux_all_cf.fullname.unique()))
                aux_all_cf = aux_all_cf.pivot(index='year', columns='fullname', values='control_factor').fillna(1).reset_index(drop=False)
                aux_all_cf['region_cd'] = irgn
                aux_all_cf = aux_all_cf.rename(columns={'region_cd' : 'region_general'})
                CF_all_sta = CF_all_sta.append(aux_all_cf, ignore_index=True)
            elif (len(str(irgn)) == 8) or (len(str(irgn)) == 7):
                aux_all_cf = control_df_all.loc[control_df_all.region_cd == irgn]
                cnt_vhc.extend(list(aux_all_cf.fullname.unique()))
                aux_all_cf = aux_all_cf.pivot(index='year', columns='fullname', values='control_factor').fillna(1).reset_index(drop=False)
                aux_all_cf['region_cd'] = irgn
                aux_all_cf = aux_all_cf.rename(columns={'region_cd' : 'region_general'})
                CF_all_cnt = CF_all_cnt.append(aux_all_cf, ignore_index=True)
            else:
                print('*****ERROR: Control factor is only applied in State level (2 digits state code e.g: 11) or')
                print('county level (8 digitis county code e.g: 11110129)')
                sys.exit()
    else:
        CF_all_sta = None
        CF_all_cnt = None
        print('')
        print('*** WARNING: There is no control factor for state level')
        print('')
    
    # Getting vehicles names which have CF to apply.
    sta_vhc = sta_vhc + cnt_vhc
    sta_vhc = list(dict.fromkeys(sta_vhc))
    
    if (CF_all_sta is not None) and (CF_all_sta.shape[0] > 0):
        # Inserting the State control factor for ALL pollutants in the final Control factor DF.
        apply_factor_sta = pd.merge(apply_factor_sta,CF_all_sta, 
                                    left_on  =['region_state', 'manufacture_date'],
                                    right_on =['region_general', 'year'], how='left').fillna(1.0)
        apply_factor_sta = apply_factor_sta.drop(columns=['region_general','year'])

    # Inserting the State control factor by Pollutant in the final Control factor DF.
    if (CF_pol_sta is not None) and (CF_pol_sta.shape[0] > 0):
        apply_factor_aux = cnt_df.county_by_yr.loc[:,['region_cd','region_state','pollutant', 'manufacture_date']].reset_index(drop=True)
        apply_factor_aux = pd.merge(apply_factor_aux,CF_pol_sta, 
                                    left_on  =['region_state', 'manufacture_date'],
                                    right_on =['region_general', 'year'], how='left')
        apply_factor_aux = apply_factor_aux.drop(columns=['region_general','year'])
    
        for ivhc in CF_pol_sta.columns.tolist():
            if (ivhc == 'year') or (ivhc == 'data') or (ivhc == 'region_general'):
                continue
            else:
                apply_factor_sta.loc[:,ivhc] = apply_factor_aux.loc[:,ivhc].fillna(apply_factor_sta.loc[:,ivhc])

    # Inserting the County control factor for ALL pollutants in the final Control factor DF.
    if (CF_all_cnt is not None) and (CF_all_cnt.shape[0] > 0):
        apply_factor_aux = cnt_df.county_by_yr.loc[:,['region_cd','region_state','pollutant', 'manufacture_date']].reset_index(drop=True)
        apply_factor_aux = pd.merge(apply_factor_aux,CF_all_cnt, 
                                    left_on  =['region_cd', 'manufacture_date'],
                                    right_on =['region_general', 'year'], how='left')
        apply_factor_aux = apply_factor_aux.drop(columns=['region_general','year'])
    
        for ivhc in CF_all_cnt.columns.tolist():
            if (ivhc == 'year') or (ivhc == 'data') or (ivhc == 'region_general'):
                continue
            else:
                apply_factor_sta.loc[:,ivhc] = apply_factor_aux.loc[:,ivhc].fillna(apply_factor_sta.loc[:,ivhc])

    # Inserting the County control factor by Pollutants in the final Control factor DF.
    if (CF_pol_cnt is not None) and (CF_pol_cnt.shape[0] > 0):
        apply_factor_aux = cnt_df.county_by_yr.loc[:,['region_cd','region_state','pollutant', 'manufacture_date']].reset_index(drop=True)
        apply_factor_aux = pd.merge(apply_factor_aux,CF_pol_cnt, 
                                    left_on  =['region_cd', 'manufacture_date'],
                                    right_on =['region_general', 'year'], how='left')
        apply_factor_aux = apply_factor_aux.drop(columns=['region_general','year'])
    
        for ivhc in CF_pol_cnt.columns.tolist():
            if (ivhc == 'year') or (ivhc == 'data') or (ivhc == 'region_general'):
                continue
            else:
                apply_factor_sta.loc[:,ivhc] = apply_factor_aux.loc[:,ivhc].fillna(apply_factor_sta.loc[:,ivhc])

    apply_factor_sta = apply_factor_sta.fillna(1.0)

    print('******************************************************************')
    print('*****   Applying Control Factor into emissions ...           *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')
    
    for ivhc in sta_vhc:
        cnt_df.county_by_yr.loc[:,ivhc] = cnt_df.county_by_yr.loc[:,ivhc] * \
                                          apply_factor_sta.loc[:,ivhc]
                    
    cnt_df.county_by_yr['total_emis'] = cnt_df.county_by_yr.loc[:,ADvhc].sum(axis=1)
    cnt_df.county_by_yr = cnt_df.county_by_yr.sort_values(by=['region_cd','pollutant'])
    
    
    cnt_df.county_total = cnt_df.county_by_yr.groupby(['region_cd','pollutant']).sum().reset_index(drop=False)
    cnt_df.county_total = cnt_df.county_total.drop(columns=['manufacture_date'])
    cnt_df.county_total = cnt_df.county_total.sort_values(by=['pollutant','region_cd']).reset_index(drop=True)
    

    cnt_df.county_total.loc[:,'region_cd'] = cnt_df.county_total.loc[:,'region_cd'].astype(str)
    cnt_df.county_total['region_state']    = cnt_df.county_total.loc[:,'region_cd'].str.slice(stop=2).astype(int)
    cnt_df.county_total.loc[:,'region_cd'] = cnt_df.county_total.loc[:,'region_cd'].astype(int)

    aux_pol = cnt_df.county_total.pollutant.unique()[0]
    aux_wgt = cnt_df.county_total_WGT.loc[(cnt_df.county_total_WGT.pollutant == aux_pol),
                                          ['region_cd','geometry','region_nm','vkt_weight']].reset_index(drop=True)

    cnt_df.county_total = pd.merge(cnt_df.county_total, aux_wgt, 
                                   left_on = ['region_cd'],
                                   right_on = ['region_cd'], how='left').reset_index(drop=True)

    #inserting the geometry and transforming it to geodataframe
    cnt_df.county_total = gpd.GeoDataFrame(cnt_df.county_total, crs=county_df.crs, geometry='geometry')

    cnt_df.county_total_WGT = cnt_df.county_total.copy()
    cnt_df.county_total_WGT.loc[:,ADvhc] = cnt_df.county_total_WGT.loc[:,ADvhc].apply(lambda x: np.asarray(x) * \
                               cnt_df.county_total_WGT.vkt_weight.values)

    #inserting the geometry and transforming it to geodataframe
    cnt_df.county_total_WGT = gpd.GeoDataFrame(cnt_df.county_total_WGT, crs=county_df.crs, geometry='geometry')

    total_county_YR_2CVS = cnt_df.county_by_yr.groupby(['region_cd', 'pollutant', 'manufacture_date']).sum().reset_index(drop=False)
    total_county_YR_2CVS.to_csv(output_dir+'/'+'District_Total_by_Year_Control_Norm_Emis_Tons_Year_{0}.csv'.format(case_name), sep=',', index=False)
    total_county_YR_2CVS = []

    total_county_2CVS = cnt_df.county_by_yr.groupby(['region_cd', 'pollutant']).sum().reset_index(drop=False)
    total_county_2CVS = total_county_2CVS.drop(columns = ['region_state','manufacture_date'])
    total_county_2CVS.to_csv(output_dir+'/'+'District_Total_Poll_Control_Norm_Emis_Tons_Year_{0}.csv'.format(case_name), sep=',', index=False)
    total_county_2CVS = []

    total_county_WGT_2CSV = cnt_df.county_total_WGT.drop(columns = ['geometry','region_state','vkt_weight'])
    total_county_WGT_2CSV.to_csv(output_dir+'/'+'District_Total_Contrl_Norm_Emis_Tons_Year_{0}.csv'.format(case_name), sep=',', index=False)
    total_county_WGT_2CSV = []
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))
    print('')
    print('******************************************************************')
    print('*****          The control factor has been applied           *****')
    print('******************************************************************')
    print('')

    return




# =============================================================================
# Function to apply chemical into emissions of PM2.5, VOC and NOX
# =============================================================================
def chemical_speciation(Input_dir, County_Emissions_DataFrame, Activity_Data_DataFrame,
                        Chemical_Speciation_Table):
    start_time = time.time()
    print('******************************************************************')
    print('*****        Starting chemical speciation of emissions       *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')
    
    input_dir = Input_dir
    cnty_emis = County_Emissions_DataFrame  #
    act_data  = Activity_Data_DataFrame     #
    ChemSpec  = Chemical_Speciation_Table   #
    
#     cnty_emis = County_Emissions    #County_Emissions_DataFrame  #
#     act_data  = AD_SK               #Activity_Data_DataFrame     #
#     ChemSpec  = Chemical_Spec_Table  #Chemical_Speciation_Table   #
# #
    nvhc = list(np.sort(act_data.fullname.vhc_name)) #vehicles from activity data

    
    pol_inv      = list(cnty_emis.county_total.pollutant.unique())
    pol_crossref = list(ChemSpec.crossref.columns[3:])
    pol_profile  = list(ChemSpec.chempro.pollutant.unique())
    
    crossref = ChemSpec.crossref.copy()
    crossref['fullname'] = crossref.vehicle.str.cat(crossref[['engine','fuel']], sep='_')
    vhc_CR = pd.DataFrame({'vhc_CR' : list(crossref['fullname'].unique())})
    cnty_emis_WGT = cnty_emis.county_total_WGT.copy()
    
    if False in list(act_data.fullname['vhc_name'].isin(vhc_CR['vhc_CR'])):
        vhc_CR = pd.merge(vhc_CR, act_data.fullname['vhc_name'], 
                          left_on='vhc_CR', right_on='vhc_name', how='right').fillna(-999)
        missing_VHC = vhc_CR['vhc_name'].loc[vhc_CR['vhc_CR'] == -999].reset_index(drop=True)
        missing_log_name = 'missing_vehicles_in_cross_refence_table.csv'
        with open(output_dir+'/LOGS/'+missing_log_name, 'w') as MCF_log:
            MCF_log.write('*****   List of missing vehicle in Cross Reference table    *****\n')
            missing_VHC.to_csv(MCF_log, index=False, header=False, line_terminator='\n')
        print('')
        print('*** ERROR: There are missing vehicles in the Cross Refence file')
        print('*** ERROR: Please, check the LOG file: {0}'.format(missing_log_name))
        print('')
        sys.exit()

    if ('PM10' in pol_inv) and ('PM2.5' in pol_inv):
        print('')
        print('***: Calculating the PMC emissions based on difference PM10 - PM2.5')
        print('')
        pmc_emis = cnty_emis_WGT.loc[cnty_emis_WGT.pollutant == 'PM10'].reset_index(drop=True)
        pmc_emis.loc[:,nvhc] = pmc_emis.loc[:,nvhc] - \
        cnty_emis_WGT.loc[cnty_emis_WGT.pollutant == 'PM2.5', nvhc].reset_index(drop=True)
        pmc_emis.loc[:,nvhc] = pmc_emis.loc[:,nvhc].clip(lower=0)
        pmc_emis['total_emis'] = pmc_emis.loc[:,nvhc].sum(axis=1)
        pmc_emis['pollutant']  = 'PMC'
        cnty_emis_WGT = cnty_emis_WGT.append(pmc_emis, ignore_index=True)
        pol_inv.extend(['PMC'])
    elif ('PM10' not in pol_inv) or ('PM2.5' not in pol_inv):
        print('')
        print('***: PM10 or PM2.5 is missing in the emissions dataframe')
        print('***: The PMC can not be calculated, please check your EF input table')
        print('***: Setting PMC to zero (0)')
        print('')
        ipmc = pol_inv[0]
        pmc_emis = cnty_emis_WGT.loc[cnty_emis_WGT.pollutant == ipmc].reset_index(drop=True)
        pmc_emis.loc[:,nvhc] = pmc_emis.loc[:,nvhc] * 0
        pmc_emis['pollutant']  = 'PMC'
        pmc_emis['total_emis'] = 0.0
        cnty_emis_WGT = cnty_emis_WGT.append(pmc_emis, ignore_index=True)
        pol_inv.extend(['PMC'])
        
    grams_sec_pol = list(ChemSpec.chempro.species.loc[ChemSpec.chempro.mw == 1].unique())
    grams_sec_pol.extend(['PM2.5', 'PM10', 'PMC'])
    
#    ipol = 'CO'  #'PM2.5'
    chem_factor_df = pd.DataFrame()
    for ipol in pol_inv:
        if (ipol in pol_crossref) and (ipol in pol_profile):
            aux_chempro = ChemSpec.chempro.loc[ChemSpec.chempro.pollutant.isin([ipol])].reset_index(drop=True)
            aux_chempro['mw_factor'] = aux_chempro.loc[:,'fraction'] / aux_chempro.loc[:,'mw']
            if False in list(crossref[ipol].isin(aux_chempro.profile)):
                miss = pd.DataFrame({'miss' : list(crossref[ipol].isin(aux_chempro.profile))})
                miss_pfl = list(crossref[ipol].loc[miss['miss'].loc[miss['miss'].isin([False])].index].unique())
                for ipfl in miss_pfl:
                    print('')
                    print('*** ERROR: Cross reference profile {0} for {1} is missing in the chemical profile'.format(ipfl,ipol))
                    print('*** ERROR: Please, add the profile {0} in the chemical profile file'.format(ipfl))
                    print('')
                    sys.exit()
            else:
                aux_chempro = aux_chempro.loc[aux_chempro.profile.isin(crossref[ipol])].reset_index(drop=True)
                aux_chempro = pd.merge(crossref.loc[:,['fullname', ipol]], aux_chempro, 
                               left_on=ipol, right_on='profile', how='left')
                aux_chempro = aux_chempro.pivot(index='species', columns='fullname', values='mw_factor').reset_index(drop=False).fillna(0.0)
                aux_chempro['pollutant'] = ipol
        elif (ipol in pol_crossref) and (ipol not in pol_profile):
            print('')
            print('*** ERROR: There is no pollutant {0} is the chemical profile file'.format(ipol))
            print('*** ERROR: Please, add a profile for the pollutant {0} in the chemical profile file'.format(ipol))
            print('')
            sys.exit()
        else:
            if ipol in list(['PM10', 'PMC']):
                print('')
                print('*** WARNING: Pollutant {0} is not in the Cross reference file'.format(ipol))
                print('*** WARNING: Molar mass for pollutant {0} is 1, because there is no need to change it to molar basis'.format(ipol))
                print('*** WARNING: Setting Molar weight to 1 (one) for the pollutant {0}'.format(ipol))
                print('')
            else:
                print('')
                print('*** WARNING: Pollutant {0} is not in the Cross reference file'.format(ipol))
                print('*** WARNING: Please, add a new column for the pollutant {0} in the cross refenrece'.format(ipol))
                print('*** WARNING: Setting Molar weight to 1 (one) for the pollutant {0}'.format(ipol))
                print('')          
            aux_chempro = crossref.loc[:,['fullname']]
            aux_chempro['species']   = ipol
            aux_chempro['mw_factor'] = 1.0
            aux_chempro = aux_chempro.pivot(index='species', columns='fullname', values='mw_factor').reset_index(drop=False).fillna(0.0)
            aux_chempro['pollutant'] = ipol

        chem_factor_df = chem_factor_df.append(aux_chempro, ignore_index=True)        


    County_Emissions_ChemSpec = pd.merge(cnty_emis_WGT,
                                         chem_factor_df.loc[:,['species', 'pollutant']], 
                                         left_on='pollutant', right_on='pollutant', how='left')

    chem_factor_df = pd.merge(cnty_emis_WGT.loc[:,['region_cd', 'pollutant']],
                              chem_factor_df,                             
                              left_on='pollutant', right_on='pollutant', how='left')

    County_Emissions_ChemSpec.loc[:,nvhc] = County_Emissions_ChemSpec.loc[:,nvhc] * \
                                            chem_factor_df.loc[:,nvhc]
    '''
    # The emissions calculated in the  previous steps is outputed as Tons/year
    # And the molar Mass is mols/grams. So, we need to bring it to grams
    # to do it, we just need to multiply by 10**6, and our unit is now mols/year
    # Later, in temporal variation, this value will converted into mols/hour
    '''
    County_Emissions_ChemSpec.loc[:,nvhc] = County_Emissions_ChemSpec.loc[:,nvhc] * 1000000
    County_Emissions_ChemSpec['total_emis'] = County_Emissions_ChemSpec.loc[:,nvhc].sum(axis=1)
   
    moles_sec_pol = pd.DataFrame({ 'species' : list(County_Emissions_ChemSpec.species.unique())})
    moles_sec_pol = moles_sec_pol.loc[~moles_sec_pol.species.isin(grams_sec_pol)].reset_index(drop=True)
    grams_sec_pol = pd.DataFrame({ 'species' : grams_sec_pol})
    
    County_Emissions_ChemSpec['total_emis'] = County_Emissions_ChemSpec.loc[:,nvhc].sum(axis=1)
    County_Emissions_ChemSpec.loc[:,County_Emissions_ChemSpec.columns != 'geometry'] = \
        County_Emissions_ChemSpec.loc[:,County_Emissions_ChemSpec.columns != 'geometry'].fillna(0.0)
    County_Emissions_ChemSpec = County_Emissions_ChemSpec.sort_values(by=['pollutant','region_cd']).reset_index(drop=True)
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))
    print('')
    print('******************************************************************')
    print('*****              Chemical Speciation is done               *****')
    print('******************************************************************')
    print('')
    return Emissions_ChemSpec(County_Emissions_ChemSpec, grams_sec_pol, moles_sec_pol)


# =============================================================================
# Function to generate hourly emissions and Gridding
# =============================================================================

def gridded_emis_NC(Activity_data_DF, Temporal_profile, County_Emissions_ChemSpec,
                    roads_RGS, run_period, GRID_info, IOAPI_out = 'N', grid4AQM = 'N'):
    
    start_time = time.time()
    print('****************************************************************')
    print('*****             Generating gridded emissions             *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    grd_AD        = Activity_data_DF
    grd_TempPro   = Temporal_profile
    ChemSpec_emis = County_Emissions_ChemSpec.chemspec_emissions.copy()
    run_period = run_period
    roads_RGS = roads_RGS
    GRID_info = GRID_info

#    grd_AD        = AD_SK
#    grd_TempPro   = TempPro
#    ChemSpec_emis = County_Emissions_ChemSpec.chemspec_emissions.copy()

    
    nvhc = list(np.sort(grd_AD.fullname.vhc_name))
    auxTP = grd_TempPro.diurnalPro.groupby(['Day','fullname']).mean().reset_index(drop=False)
    '''
    Dividing the emissions by 3600 to bring it from grams or mols per hour to
    grams or mols per second. This is necessary because Air quality models
    expects gridded emissions in this unit.
    '''
    ChemSpec_emis.loc[:,nvhc] = ChemSpec_emis.loc[:,nvhc] / 3600 
    
    diurnalPro_out = pd.DataFrame(columns=['DateTime']+list(np.sort(auxTP.fullname.unique())))
    for iday in run_period.Day.unique():
        auxDP = auxTP.loc[auxTP.Day == iday].T
        auxDP.columns = list(auxDP.loc['fullname'])
        auxDP = auxDP.drop(index=['Day','fullname']).reset_index(drop=True)
        auxDP = auxDP.sort_index(axis=1)
        auxDP.insert(loc=0, column='DateTime', value=pd.date_range(start=iday, freq='H', periods=24))
        
        diurnalPro_out = diurnalPro_out.append(auxDP,ignore_index=True, sort=False)
    
    diurnalPro_vhc = diurnalPro_out.columns.tolist()
    diurnalPro_vhc.remove('DateTime')
    missing_log_name = 'LOG_Missing_Vehicle_in_Temporal_Cross_reference.txt'
    with open(output_dir+'/LOGS/'+missing_log_name, 'w') as MCF_log:
        MCF_log.write('*****   List of missing vehicle in Temporal Cross Reference table    *****\n')
        cnt = 0
        for i in nvhc:
            if i not in diurnalPro_vhc:
                cnt += 1
                MCF_log.write('{0}\n'.format(i))
        if cnt > 0:
            print('')
            print('*** ERROR: There are missing vehicles in the Temporal Cross Refence file')
            print('*** ERROR: Please, check the LOG file: {0}'.format(missing_log_name))
            print('')
            sys.exit()
        else:
            MCF_log.write('*****  End of List of missing vehicle in Temporal Cross Reference table *****\n')


    
    aux_hourly = ChemSpec_emis.loc[:,['region_cd','species']+list(grd_AD.fullname.vhc_name)]
    hourly_emissions = ChemSpec_emis.loc[:,['region_cd','species']]
    for idt in diurnalPro_out.DateTime:
        auxDT = diurnalPro_out.loc[diurnalPro_out.DateTime == idt, nvhc].reset_index(drop=True)
        auxEmis = aux_hourly.copy()
        for icol in auxDT.columns:
            auxEmis[icol] = aux_hourly[icol] * auxDT[icol].values
        hourly_emissions.insert(loc=hourly_emissions.shape[1], column=idt, value=auxEmis.loc[:,nvhc].sum(axis=1))
    
    gridded_emissions = pd.DataFrame()
    for ipol in hourly_emissions.species.unique():
    
        grd_all = pd.merge(roads_RGS.surrogate.loc[:,['grid_id','region_cd','weight_factor']],
                       hourly_emissions.loc[hourly_emissions.species == ipol], left_on='region_cd', right_on='region_cd', how='left').fillna(0)
        
        for icol in run_period.DateTime:
            grd_all.loc[:,icol] = grd_all.loc[:,icol] * grd_all.weight_factor.values
        
        grd_all.loc[:,'species'] = ipol
        grd_all = grd_all.groupby(['grid_id','species']).sum().drop(columns=['weight_factor','region_cd']).reset_index(drop=False)
        grd_all = pd.merge(grd_all,roads_RGS.grid.loc[:,['geometry','grid_id','col','row']], left_on='grid_id', right_on='grid_id', how='left')
        grd_all = gpd.GeoDataFrame(grd_all, crs=roads_RGS.grid.crs, geometry='geometry')
        
        grid_id_list = list(grd_all.grid_id.unique())
        gridded_emissions = gridded_emissions.append(grd_all, ignore_index=True, sort=False)
    
    drop_zero = gridded_emissions.loc[gridded_emissions.species == 0].index
    gridded_emissions = gridded_emissions.drop(index=list(drop_zero)).reset_index(drop=True)
    #releasing memory of gridding dataframe
    grd_all = []
    
    if (IOAPI_out.lower() == 'yes') or (IOAPI_out.lower() == 'y')  and \
       (grid4AQM.lower() == 'yes')  or (grid4AQM.lower() == 'y'):
        
        print('')
        print('******************************************************************')
        print('*****     Generating chemical speciation IOAPI files         *****')
        print('******************************************************************')
        print('')
 
        #getting the number of vars and the VAR-LIST
        NVARS = len(gridded_emissions.species.unique())
        VAR_LIST = []
        for kpol in gridded_emissions.species.unique():
            VAR_LIST.append('{:<16}'.format(kpol))
        s = pd.Series(VAR_LIST)
        VAR_LIST = s.str.cat(sep='')
        df_times = run_period.DateTime
        ntflag = run_period.TFLAG.to_list()
        TFLAG_list = np.zeros((df_times.shape[0], NVARS, 2))
        
        for i in range(0,NVARS):
            TFLAG_list[:,i,:] = ntflag
        TFLAG_list = TFLAG_list.astype(np.int32)
        
        #getting the FILEDESC
        filedesc = ['{:<80.80}'.format('Vehicle emissions inventory'),
                    '{:<80.80}'.format('/FROM/ CARS - Comprehensive Automobile Emissions Research Simulator'),
                    '{:<160.160}'.format('/VERSION/ CARSv1'),
                    '{:<80.80}'.format('/BASE YEAR/     {}'.format(run_period.DateTime[0].year)),
                    '{:<80.80}'.format('/NUMBER OF FILES/   1'),
                    '{:<80.80}'.format('/FILE POSITION/   1'),
                    '{:<4240.4240}'.format('/NUMBER OF VARIABLES/  {}'.format(NVARS))]
        f = pd.Series(filedesc)
        filedesc = f.str.cat(sep='')
        
        now = dt.datetime.now()
        col = gridded_emissions.col.max()
        row = gridded_emissions.row.max()
        
        # openign the netcdf file
        ncname = (output_dir+'/vehicle_emissions_{0}_n{1}_{2}.nc'.format(run_period.jul_day.iloc[0],run_period.shape[0],case_name))
        
        ncfile = Dataset(ncname, 'w', format='NETCDF3_CLASSIC')
        ncfile.createDimension('TSTEP' , None)
        ncfile.createDimension('DATE-TIME', 2)
        ncfile.createDimension('LAY' , 1)
        ncfile.createDimension('VAR'   , NVARS)
        ncfile.createDimension('ROW'   , row)
        ncfile.createDimension('COL'   , col)
        
        TFLAG_var = ncfile.createVariable('TFLAG' , np.int32 , ('TSTEP', 'VAR', 'DATE-TIME'))
        TFLAG_var.units     = '{:<16.16}'.format("<YYYYDDD,HHMMSS>")
        TFLAG_var.long_name = '{:<16.16}'.format("TFLAG")
        TFLAG_var.var_desc  = '{:<80.80}'.format("Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS")
        TFLAG_var[:] = TFLAG_list
        
        for ipol in gridded_emissions.species.unique():
        
            pol_grd = gridded_emissions.loc[gridded_emissions.species == ipol].reset_index(drop=True)
            pol_grd = pol_grd.drop(columns=['grid_id','species', 'geometry'])
        
            pol_array = np.zeros((df_times.shape[0], 1, row,col))
            
            for idx ,idt in zip(range(0,df_times.shape[0]), df_times):
                aux_grd = pol_grd.loc[:,idt].to_numpy()
                aux_grd = aux_grd.reshape((row,col)) #np.flipud(aux_grd.reshape((row,col)))
                pol_array[idx,0,:,:] = aux_grd
            
            pol_array = pol_array.astype(np.float32)
            
            nc_var    = ncfile.createVariable(ipol , np.float32 , ('TSTEP','LAY', 'ROW','COL'))
            nc_var[:] = pol_array
            if ipol in County_Emissions_ChemSpec.grams_pol.species.to_list():
                nc_var.long_name = '{:<16.16}'.format('{}'.format(ipol))
                nc_var.units     = '{:<16.16}'.format('g/s')
                nc_var.var_desc  = '{:<80.80}'.format('Model species {}'.format(ipol))
        
            else:
                nc_var.long_name = '{:<16.16}'.format('{}'.format(ipol))
                nc_var.units     = '{:<16.16}'.format('moles/s')
                nc_var.var_desc  = '{:<80.80}'.format('Model species {}'.format(ipol))
        
        
        ncfile.IOAPI_VERSION = '{:<80.80}'.format("$Id: @(#) ioapi library version 3.0 $") 
        ncfile.EXEC_ID       = '{:<80.80}'.format("????????????????")
        ncfile.FTYPE         = 1 
        ncfile.CDATE         = int('{}{:0>3}'.format(str(now.year),str(now.day))) 
        ncfile.CTIME         = int('{}{}{}'.format(now.hour,now.minute,now.second)) 
        ncfile.WDATE         = int('{}{:0>3}'.format(str(now.year),str(now.day))) 
        ncfile.WTIME         = int('{}{}{}'.format(now.hour,now.minute,now.second)) 
        ncfile.SDATE         = run_period.jul_day[0] 
        ncfile.STIME         = run_period.jul_hour[0] 
        ncfile.TSTEP         = 10000 
        ncfile.NTHIK         = GRID_info.NTHIK
        ncfile.NCOLS         = GRID_info.NCOLS
        ncfile.NROWS         = GRID_info.NROWS
        ncfile.NLAYS         = GRID_info.NLAYS
        ncfile.NVARS         = NVARS
        ncfile.GDTYP         = GRID_info.GDTYP
        ncfile.P_ALP         = GRID_info.P_ALP
        ncfile.P_BET         = GRID_info.P_BET
        ncfile.P_GAM         = GRID_info.P_GAM 
        ncfile.XCENT         = GRID_info.XCENT 
        ncfile.YCENT         = GRID_info.YCENT 
        ncfile.XORIG         = GRID_info.XORIG 
        ncfile.YORIG         = GRID_info.YORIG 
        ncfile.XCELL         = GRID_info.XCELL
        ncfile.YCELL         = GRID_info.YCELL
        ncfile.VGTYP         = 7  #-1 ?
        ncfile.VGTOP         = np.float32(0.) #0.f ??
        ncfile.VGLVLS        = (np.float32(0.), np.float32(0.))
        ncfile.GDNAM         = "{:<16}".format(GRID_info.GDNAM)
        ncfile.UPNAM         = "{:<16}".format("OPENSET")
        ncfile.VAR_LIST      = VAR_LIST
        ncfile.FILEDESC      = filedesc
        ncfile.HISTORY       = now.ctime() 
        ncfile.renameAttribute('VAR_LIST','VAR-LIST')
        
        #closing netcdf
        ncfile.close()
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))
    print('')
    print('******************************************************************')
    print('*****              Gridded emissions is done               *****')
    print('******************************************************************')
    print('') 
    
    return gridded_emissions




def plot_figures_chart(County_Emissions, County_Emissions_ChemSpec, Activity_data_DF, 
                       gridded_emissions, roads_RGS, avgSpeedDist, run_period,
                       plot_figures = 'no', plot_24 = 'no', pol_list = [], species_list = [],
                       adj_scale = 0.4):
    print('')
    print('******************************************************************')
    print('*****              The Plot option is YES.                   *****')
    print('***** The total emissions of CO, NOX, PM10, PM2.5, VOC and SOX ***')
    print('******************************************************************')
    print('')  
    start_time_plot = time.time()
    
    County_Emissions = County_Emissions
    grams_pol = County_Emissions_ChemSpec.grams_pol.species.to_list()
    moles_pol = County_Emissions_ChemSpec.moles_pol.species.to_list()
    AD_SK = Activity_data_DF
    gridded_emissions = gridded_emissions
    roads_RGS = roads_RGS
    run_period = run_period
    
    # =============================================================================
    # Function to create the GrayJet colormap
    # =============================================================================
    def createGrayJet():
        from matplotlib import cm
        from matplotlib.colors import ListedColormap
        import numpy as np
        jet = cm.get_cmap('jet', 256)
        newcolors = jet(np.linspace(0, 1, 256))
        newcolors[0, :] = np.array([256/256, 256/256, 256/256, 1])    #white
        newcolors[1:4, :] = np.array([200/256, 200/256, 200/256, 1])  #lightgray
        GrayJet = ListedColormap(newcolors)
        return GrayJet

    GrayJet = createGrayJet()
    
    if len(pol_list) == 0:
        polls   = list(County_Emissions.county_total_WGT.pollutant.unique())
    else:
        polls = pol_list
    if len(species_list) == 0:
        species = list(gridded_emissions.species.unique())
    else:
        species = species_list
    
    if plot_figures == 'yes':
        
        for ipol in polls:
            nvhc = list(np.sort(AD_SK.fullname.vhc_name))
            total_by_county = County_Emissions.county_total_WGT.loc[County_Emissions.county_total_WGT.pollutant == ipol]
            total_by_county = total_by_county.reset_index(drop=True)
            total_by_county = total_by_county.drop(columns = ['geometry'])
            Emissions_by_link = pd.merge(roads_RGS.roads_df,total_by_county, left_on = ['region_cd'], 
                                right_on = ['region_cd'], how='left')
            
            Emissions_by_link.loc[:,nvhc] = Emissions_by_link.loc[:,nvhc].apply(lambda x: np.asarray(x) * \
                                 Emissions_by_link.vkt_split_county.values) 
            Emissions_by_link['total_emis']      = Emissions_by_link.loc[:,nvhc].sum(axis=1)    
            ccc = Emissions_by_link.loc[(Emissions_by_link.pollutant == ipol)]
            ccc['total_emis'] = ccc.loc[:,AD_SK.fullname.vhc_name].sum(axis=1)
            ccc.loc[:,'total_emis'] = ccc.loc[:,'total_emis']
            fig, ax = plt.subplots(1, figsize = (13,8))
            ccc.plot(column='total_emis', cmap=GrayJet, linewidth=1.0, ax=ax,
                               vmax = ccc.total_emis.max()* adj_scale)# ,edgecolor='0.8')
            
            ax.set_title('Total Emissions of {0}-{1} [ton/yr]'.format(case_name,ipol), fontsize = 22)
            sm = plt.cm.ScalarMappable(cmap=GrayJet, 
                                       norm=plt.Normalize(vmin=ccc.total_emis.min(),
                                                          vmax=ccc.total_emis.max() * adj_scale))
            sm._A = []
            cbar = fig.colorbar(sm)
            #    ax.axis('off')
            name = '{0}_total_link_emissions_{1}.png'.format(case_name,ipol)
            fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
            plt.close()
        
        for ipol in polls:
            ddd = County_Emissions.county_total_WGT.loc[County_Emissions.county_total_WGT.pollutant == ipol].reset_index(drop=True)
            fig, ax = plt.subplots(1, figsize = (13,8))
            ddd.plot(column='total_emis', cmap=GrayJet, linewidth=0.8, ax=ax,
                               vmax = ddd.total_emis.max()* adj_scale)# ,edgecolor='0.8')
            
            ax.set_title('Total Emissions of {0}-{1} [ton/yr]'.format(case_name,ipol), fontsize = 22)
            sm = plt.cm.ScalarMappable(cmap=GrayJet, 
                                       norm=plt.Normalize(vmin=ddd.total_emis.min(),
                                                          vmax=ddd.total_emis.max() * adj_scale))
            sm._A = []
            cbar = fig.colorbar(sm)
            #    ax.axis('off')
            name = '{0}_total_county_emissions_{1}.png'.format(case_name,ipol)
            fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
            plt.close()
    
        #Piechart Plots
        fuels = list(AD_SK.fuels.fuels)
        for ipol in polls:
            pol_df = County_Emissions.county_total_WGT.loc[County_Emissions.county_total_WGT.pollutant == ipol].reset_index(drop=True)
            plot_df_fuel = pd.DataFrame()
            for ifuel in fuels:  #make a correction on the columns extraction to separate h-gasoline from gasoline
                fuel_df = pol_df.loc[:,pol_df.columns.str.contains(ifuel)]
                aux_df = pd.DataFrame({ifuel : [fuel_df.sum().sum()]})
                plot_df_fuel = pd.concat([plot_df_fuel,aux_df],axis=1,  sort=False)
                
            plot_df_fuel = plot_df_fuel.T
            plot_df_fuel = plot_df_fuel.rename(columns={ plot_df_fuel.columns[0]: "Fuels" })
            fig, ax = plt.subplots(1, figsize = (13,8))
    
            plot_df_fuel.plot(y='Fuels', kind='pie',autopct='%.0f%%', ax=ax, legend=True, fontsize=14)
            ax.legend(loc='center left', bbox_to_anchor=(-0.3, 0.25, 0.5, 0.5), fontsize='large')            
            ax.set_title('Total Emissions by Fuel - {0}'.format(ipol), fontsize = 22)
    
            name = 'Pie_Chart_Emissions_by_Fuel_{0}.png'.format(ipol)
            fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
            plt.close()

        #Average Speed profile
        from matplotlib.pyplot import cm
        fntsz = 16
        roads = list(avgSpeedDist.data.columns)
        roads.remove('speed')
        plt.style.use('bmh')
        maker = ['*','^','o','+','h','x','s','>']
        color=iter(cm.rainbow(np.linspace(0,1,len(roads))))
        
        fig, ax = plt.subplots(1, figsize = (10,4))
        rlbl =  {'101': '101 : National Expressway',
                  '102': '102 : Urban Expressway',
                  '103': '103 : Highway',
                  '104': '104 : Urban Expressway',
                  '105': '105 : National Local Road',
                  '106': '106 : Rural Local Road',
                  '107': '107 : Urban Local Road',
                  '108': '108 : Local Road'}
        
        for idx in np.arange(0, len(roads)):
            print(idx, roads[idx])
            c=next(color)
            m=maker[idx]
            p2 = plt.plot(avgSpeedDist.data['speed'], 
                          avgSpeedDist.data[roads[idx]], color = c, marker=m, label=rlbl[roads[idx]])
        plt.ylabel('Fraction')
        plt.title('Average Speed Profile by road')
        plt.yticks(fontsize=fntsz)#
        plt.legend()
        name = '{0}_Average_Speed_Profile_by_road.png'.format(case_name)
        fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
        plt.close()

        
        #Gridded emissions Plots
        for ipol in species:
            plot_pol = gridded_emissions[gridded_emissions.species == ipol].reset_index(drop=True)
            plot_pol['total_emis'] = plot_pol.loc[:,run_period.DateTime.loc[0:24]].sum(axis=1)
            vmin = 0
            vmax = adj_scale * plot_pol['total_emis'].max()
    
            fig, ax = plt.subplots(1, figsize = (13,8))
            plot_pol.plot(column='total_emis', cmap=GrayJet, linewidth=0.8, ax=ax,
                               vmax = vmax)# ,edgecolor='0.8')
            
            ax.set_title('{0} Grid Emissions - {1}'.format(case_name,ipol), fontsize = 22)
            sm = plt.cm.ScalarMappable(cmap=GrayJet, 
                                       norm=plt.Normalize(vmin=vmin,
                                                          vmax=vmax))
            sm._A = []
            cbar = fig.colorbar(sm)
            particles = grams_pol
            if ipol in particles:
                cbar.ax.set_ylabel(r'$grams{\cdot}day^{-1}$', fontsize=14)
            else:
                cbar.ax.set_ylabel(r'$mols{\cdot}day^{-1}$', fontsize=14)
    
            name = '{0}_total_grid_emissions_{1}.png'.format(case_name,ipol)
            fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
            plt.close()
    
    #Hourly gridded emissions Plots
    if plot_24 == 'yes':
        for ipol in species:
            plot_pol = gridded_emissions[gridded_emissions.species == ipol].reset_index(drop=True)
            vmin = 0
            vmax = adj_scale * plot_pol.loc[:,run_period.DateTime.loc[0:24]].max().max()
            for itime in run_period.DateTime.loc[0:24]:
                fig, ax = plt.subplots(1, figsize = (13,8))
                plot_pol.plot(column=plot_pol[itime], cmap=GrayJet, linewidth=0.8, ax=ax,
                              vmax = vmax)# ,edgecolor='0.8')
    			
                ax.set_title('Grid Emissions of {0} - {1}'.format(case_name,ipol), fontsize = 22)
                sm = plt.cm.ScalarMappable(cmap=GrayJet, 
                                           norm=plt.Normalize(vmin=vmin, vmax=vmax))
                xpos = roads_RGS.surrogate.total_bounds[0] + 2000
                ypos = roads_RGS.surrogate.total_bounds[3] - 2000
                ax.text(xpos, ypos, '{0}:00H'.format(itime.hour), color='red', bbox=dict(facecolor='white', alpha=1))
                sm._A = []
                cbar = fig.colorbar(sm)
                particles = grams_pol
                
                if ipol in particles:
                    cbar.ax.set_ylabel(r'$grams{\cdot}s^{-1}$', fontsize=14)
                else:
                    cbar.ax.set_ylabel(r'$mols{\cdot}s^{-1}$', fontsize=14)

                if itime.hour < 10:
                    name = '{0}_hourly_grid_emissions_{1}_hour_0{2}.png'.format(case_name,ipol,itime.hour)
                else:
                    name = '{0}_hourly_grid_emissions_{1}_hour_{2}.png'.format(case_name,ipol,itime.hour)
                fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
                plt.close()
        
    total_plot = ((time.time() - start_time_plot))
    print('')  
    print('---  CARS Plot Elapsed time in seconds = {0}     ---'.format(total_plot))
    print('---  CARS Plot Elapsed time in minutes = {0}     ---'.format(total_plot/60))
    print('---  CARS Plot Elapsed time in   hours = {0}     ---'.format(total_plot/3600))    
    print('')
    print('******************************************************************')
    print('*****                CARS finished to plot                   *****')
    print('******************************************************************')
    print('')  

    return 






def main():
    print('')
    print('******************************************************************')
    print('** Comprehensive Automobile Emissions Research Simulator - CARS **')
    print('******************************************************************')
    print('')
    print('*****    Loading classes, functions and primary input data   *****')
    print('*****                      Please wait                       *****')
    print('')

    GrayJet = createGrayJet()
    run_period = run_period_times(STDATE, ENDATE, STTIME, RUNLEN)
    
    #Reading chemical profile csv file. Calling read_chemical_info function
    avgSpeedDist = read_avgSpeedDistribution(input_dir, avg_SPD_Dist_file)


    #Reading Activity data e csv file. Calling read_activity_data_csv_SK function
    if (process_activity_data.lower() == 'yes') or (process_activity_data.lower() == 'y'):
        AD_SK = read_activity_data_csv_SK(input_dir, activity_file, ENDATE)
        
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
        AD_SK = read_activity_data_csv_SK(input_dir, activity_file, ENDATE)
    
    #Reading and calculating Emissions factor. Calling read_emissions_factor_SK and calculate_EmisFact function
    if (process_emissions_factor.lower() == 'yes') or (process_emissions_factor.lower() == 'y'):
        EF_All_fuel_function = read_emissions_factor_SK(input_dir, Emis_Factor_list)
        EmisFactor_yr_spd    = calculate_EmisFact(input_dir, EF_All_fuel_function, output_dir)
        
        if (Deterioration_factor.lower() == 'yes') or (Deterioration_factor.lower() == 'y'):
            EmisFactor_yr_spd = apply_deterioration_emissions_factor_SK(input_dir, 
                                                                    Deterioration_list, 
                                                                    EmisFactor_yr_spd, output_dir)
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
        EmisFactor_yr_spd    = calculate_EmisFact(input_dir, EF_All_fuel_function, output_dir)
        
        if (Deterioration_factor.lower() == 'yes') or (Deterioration_factor.lower() == 'y'):
            EmisFactor_yr_spd = apply_deterioration_emissions_factor_SK(input_dir, 
                                                                    Deterioration_list, 
                                                                    EmisFactor_yr_spd, output_dir)

    #Comparing vehicles between Activity data Emissions factor. Calling check_VHC_AD_EF
    check_VHC_AD_EF(EmisFactor_yr_spd, AD_SK, output_dir)


    #Reading and processinf Roads and county shape file. Calling roads_RGS and county_SHP
    if (process_shapes.lower() == 'yes') or (process_shapes.lower() == 'y'):
        roads_RGS = roads_grid_surrogate_inf(input_dir, link_shape, link_shape_att,
                                     gridfile_name, Radius_Earth, Unit_meters = True ) 
        county_SHP = processing_County_shape(input_dir, county_shape, county_shape_att, gridfile_name, Radius_Earth,) 
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
            county_SHP = processing_County_shape(input_dir, county_shape, county_shape_att, gridfile_name, Radius_Earth,) 
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
        roads_RGS = roads_grid_surrogate_inf(input_dir, link_shape, link_shape_att,
                             gridfile_name, Radius_Earth, Unit_meters = True ) 
        county_SHP = processing_County_shape(input_dir, county_shape, county_shape_att, gridfile_name, Radius_Earth,) 

    # Calculating the county emissions
    if (calculate_county_emissions.lower() == 'yes') or (calculate_county_emissions.lower() == 'y'):
        County_Emissions = calc_County_Emissions(input_dir, EmisFactor_yr_spd,
                                           AD_SK, roads_RGS, county_SHP,      
                                           Cold_Start_list, avgSpeedDist, output_dir)
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
        County_Emissions = calc_County_Emissions(input_dir, EmisFactor_yr_spd,
                                   AD_SK, roads_RGS, county_SHP,      
                                   Cold_Start_list, avgSpeedDist, output_dir)

    # Applying control strategy
    if (control_strategy.lower() == 'yes') or (control_strategy.lower() == 'y'):
        aplly_control(control_list, County_Emissions, county_SHP)
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

        gridded_emissions = gridded_emis_NC(AD_SK, TempPro, County_Emissions_ChemSpec,
                              roads_RGS, run_period, GRID_info, IOAPI_out = 'Y', grid4AQM = 'Y')
        
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

        gridded_emissions = gridded_emis_NC(AD_SK, TempPro, County_Emissions_ChemSpec,
                              roads_RGS, run_period, GRID_info, IOAPI_out = 'N', grid4AQM = 'N')

    if (plot_figures.lower() == 'yes') or (plot_figures.lower() == 'y'):
        plot_figures_chart(County_Emissions, County_Emissions_ChemSpec, AD_SK, gridded_emissions, 
                           roads_RGS, avgSpeedDist, run_period,
                           plot_figures = plot_figures, plot_24 = plot_24, 
                           pol_list = [], species_list = ['NO2', 'CO'])
    else:
        print('')
        print('******************************************************************')
        print('*****            plot_figures is set to NO                   *****')
        print('******************************************************************')
        print('')
    
   
if __name__== "__main__":
   main()

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






# from matplotlib.pyplot import cm

# fntsz = 16

# roads = list(avgSpeedDist.data.columns)
# roads.remove('speed')
# plt.style.use('bmh')
# maker = ['*','^','o','+','h','x','s','>']
# color=iter(cm.rainbow(np.linspace(0,1,len(roads))))

# fig, ax = plt.subplots(1, figsize = (10,4))
# rlbl =  {'101': '101 : National Expressway',
#           '102': '102 : Urban Expressway',
#           '103': '103 : Highway',
#           '104': '104 : Urban Expressway',
#           '105': '105 : National Local Road',
#           '106': '106 : Rural Local Road',
#           '107': '107 : Urban Local Road',
#           '108': '108 : Local Road'}

# for idx in np.arange(0, len(roads)):
#     print(idx, roads[idx])
#     c=next(color)
#     m=maker[idx]
#     p2 = plt.plot(avgSpeedDist.data['speed'], 
#                   avgSpeedDist.data[roads[idx]], color = c, marker=m, label=rlbl[roads[idx]])

# plt.ylabel('weight')
# plt.title('Average Speed Profile by road')
# plt.yticks(fontsize=fntsz)#
# plt.legend()

# name = '{0}_Average_Speed_Profile_by_road.png'.format(case_name)
# fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
# plt.close()








# fuels = list(AD_SK.fuels.fuels)

# pol_df = County_Emissions.county_total.copy()
# plot_df_fuel = pd.DataFrame()
# for ifuel in fuels:  #make a correction on the columns extraction to separate h-gasoline from gasoline
#     fuel_df = pd.concat([pol_df.loc[:,['pollutant']]
#                         ,pol_df.loc[:,pol_df.columns.str.contains('_'+ifuel)]], axis=1)
#     fuel_df = fuel_df.groupby('pollutant').sum().reset_index(drop=False)
#     fuel_df = pd.DataFrame({'pollutant' : fuel_df.pollutant, 
#                            ifuel : fuel_df.loc[:,fuel_df.columns.str.contains(ifuel)].sum(axis=1)})
#     if plot_df_fuel.shape[1] == 0:
#         plot_df_fuel = pd.concat([plot_df_fuel,fuel_df], sort=False)
#     else:
#         plot_df_fuel = pd.merge(plot_df_fuel,fuel_df, left_on='pollutant', right_on='pollutant', how='outer')

# plot_df_fuel['hybrid'] = plot_df_fuel.loc[:,plot_df_fuel.columns.str.contains('h-')].sum(axis=1)
# plot_df_fuel = plot_df_fuel.loc[:,~plot_df_fuel.columns.str.contains('h-')]

# plot_df_fuel =  plot_df_fuel.set_index('pollutant')
# plot_df_fuel = plot_df_fuel.T

# Xcols = plot_df_fuel.columns.to_list()
# Ycols = plot_df_fuel.index.tolist()
# width = 0.5      # the width of the bars: can also be len(x) sequence
# fntsz = 12
# plt.style.use('bmh')
# fig, ax = plt.subplots(1, figsize = (10,4))
# p1 = plt.bar(Xcols, plot_df_fuel.loc[Ycols[0],:], width, label=Ycols[0])
# bot = plot_df_fuel.loc[Ycols[0],:].tolist()

# for idx in np.arange(1, len(Ycols)):
#     print(idx, Ycols[idx])
#     p2 = plt.bar(Xcols, plot_df_fuel.loc[Ycols[idx],:], width, bottom=bot, label=Ycols[idx])
#     bot = np.add(bot, plot_df_fuel.loc[Ycols[idx],:]).tolist()

# y = plot_df_fuel.sum().to_list()
# for i, v in enumerate(y):
#     plt.text(i -0.25, v + 1, str(int(np.round(v))))


# # plt.yscale('log')
# plt.ylabel('Tons/year')
# plt.title('Emissions by fuel type')
# plt.yticks(fontsize=fntsz)#

# plt.legend()

# name = '{0}_emissions_by_fuel_pollutant.png'.format(case_name)
# fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
# plt.close()








# groupY = list(np.arange(1985, 2017,2))
# groupY.append(2017)
# AAA = AD_SK.vhc_count.groupby('manufacture_date').sum().reset_index(drop=False)
# AAA = pd.merge(AAA, pd.DataFrame({'group' : groupY}), left_on='manufacture_date', 
#                right_on='group', how='left').bfill()
# AAA = AAA.groupby('group').sum()
# AAA = AAA.drop(columns='manufacture_date')

# BBB = AAA.stack().reset_index(drop=False)
# BBB['fuel'] = [x.split('_')[-1] for x in BBB.level_1]
# BBB['fuel'] = [x.split('-')[0].replace('h','hybrid') for x in BBB.fuel]
# BBB = BBB.groupby(['group','fuel']).sum().reset_index(drop=False)
# BBB = BBB.rename(columns={0: 'count'})
# BBB = BBB.pivot(index='fuel', columns='group', values='count')

# N = len(groupY)
# fuels = list(BBB.index.unique())
# width = 1.5       # the width of the bars: can also be len(x) sequence
# fntsz = 12
# plt.style.use('bmh')
# fig, ax = plt.subplots(1, figsize = (10,4))
# p1 = plt.bar(groupY, BBB.loc[fuels[0],:], width, label=fuels[0])
# bot = BBB.loc[fuels[0],:].tolist()

# for idx in np.arange(1, len(fuels)):
#     print(idx, fuels[idx])
#     p2 = plt.bar(groupY, BBB.loc[fuels[idx],:], width, bottom=bot, label=fuels[idx])
#     bot = np.add(bot, BBB.loc[fuels[idx],:]).tolist()

# plt.ylabel('#')
# plt.title('Vehicles age distribuition by fuel')
# plt.xticks(groupY, tuple([str(x) for x in groupY]), fontsize=fntsz, rotation=90)
# plt.yticks(fontsize=fntsz)#

# plt.legend()

# name = '{0}_vehicle_age_distribuition_fuel.png'.format(case_name)
# fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
# plt.close()


# CCC = AAA.stack().reset_index(drop=False)
# CCC['vhc'] = [x.split('_')[0] for x in CCC.level_1]
# CCC = CCC.groupby(['group','vhc']).sum().reset_index(drop=False)
# CCC = CCC.rename(columns={0: 'count'})
# CCC = CCC.pivot(index='vhc', columns='group', values='count')

# N = len(groupY)
# vhc = list(CCC.index.unique())
# width = 1.5       # the width of the bars: can also be len(x) sequence
# fntsz = 12
# plt.style.use('bmh')
# fig, ax = plt.subplots(1, figsize = (10,4))
# p1 = plt.bar(groupY, CCC.loc[vhc[0],:], width, label=vhc[0])
# bot = CCC.loc[vhc[0],:].tolist()

# for idx in np.arange(1, len(vhc)):
#     print(idx, vhc[idx])
#     p2 = plt.bar(groupY, CCC.loc[vhc[idx],:], width, bottom=bot, label=vhc[idx])
#     bot = np.add(bot, CCC.loc[vhc[idx],:]).tolist()

# plt.ylabel('#')
# plt.title('Vehicles age distribuition by type')
# plt.xticks(groupY, tuple([str(x) for x in groupY]), fontsize=fntsz, rotation=90)
# plt.yticks(fontsize=fntsz)#

# plt.legend()

# name = '{0}_vehicle_age_distribuition_type.png'.format(case_name)
# fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
# plt.close()




# #Average Speed profile
# from matplotlib.pyplot import cm
# fntsz = 16
# roads = list(ASD.columns)
# roads.remove('speed')
# plt.style.use('bmh')
# maker = ['*','^','o','+','h','x','s','>']
# color=iter(cm.rainbow(np.linspace(0,1,len(roads))))

# fig, ax = plt.subplots(1, figsize = (10,4))
# rlbl =  {'101': '101 : National Expressway',
#           '102': '102 : Urban Expressway',
#           '103': '103 : Highway',
#           '104': '104 : Urban Expressway',
#           '105': '105 : National Local Road',
#           '106': '106 : Rural Local Road',
#           '107': '107 : Urban Local Road',
#           '108': '108 : Local Road'}

# for idx in np.arange(0, len(roads)):
#     print(idx, roads[idx])
#     c=next(color)
#     m=maker[idx]
#     p2 = plt.plot(ASD['speed'], 
#                   ASD[roads[idx]], color = c, marker=m, label=rlbl[roads[idx]])
# plt.ylabel('Fraction')
# plt.title('Average Speed Profile by road')
# plt.yticks(fontsize=fntsz)#
# plt.legend()
# name = '{0}_Average_Speed_Profile_by_road_RAW.png'.format(case_name)
# fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=110)
# plt.close()

