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
import pyproj, csv
import calendar
import matplotlib.pyplot as plt
from netCDF4 import Dataset
#matplotlib.use('Agg')
total_start_time = time.time()

'''
*****      Case name and path to the folders - Please, set for all directories      *****
'''
case_name = 'Country_KNU_09_01' #'_Seoul_Old_AD'#
home_dir   = r'C:/Users/pedruzzi/OneDrive - University of North Carolina at Chapel Hill/0000_EUA_IE/001_mobile_source_inventory/CARS_source_code/CARS_test_Case_Country_KNU_09_01'
src_dir    = home_dir+'/src'
input_dir  = home_dir+'/input_country'
inter_dir  = home_dir+'/intermediate'#+case_name
output_dir = home_dir+'/output_country'#+case_name
met_dir    = input_dir+'/metdata'

'''
*****      Dates and Runlenght      *****
'''
STDATE = '2017-01-01'  # start date
ENDATE = '2017-01-02'  # end date
STTIME =  00            # start time 
RUNLEN =  24            # run length


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


activity_file          =  'cars_input_v3.1.csv' #'seoul_2017_OLD.csv' #

Emis_Factor_list       = ['cng_v3.csv','diesel_v3.csv','gasoline_v3.csv','lpg_v3.csv',
                          'H-CNG_v3.csv','H-Diesel_v3.csv','H-Gasoline_v3.csv' ,'H-LPG_v3.csv'] #['gasoline.csv','diesel.csv','cng.csv','lpg.csv']

Cold_Start_list        = ['cold_start_vehicles.csv']

avg_SPD_Dist_file      = 'avgSpeedDistribution_rev_01.csv'

link_shape             = '/shapes'+'/country_road_by_county_GRS80_Avg_VKT.shp'
link_shape_att         = ['link_id'  , 'EMD_CD' , 'EMD_ENG_NM', 'EMD_KOR_NM',
                           'link_type', 'speed', 'length', 'Avg_VKT']
                       
county_shape           = '/shapes/TL_SCCO_EMD.shp'
county_shape_att       = ['EMD_CD', 'EMD_ENG_NM', 'EMD_KOR_NM']

'''
Do you want to apply Deterioration factor into emissions factors?
If YES, set Deterioration_CALC     = 'yes' and insert the deterioration
list file. Remember that vehciles names should match with emissions factor
and Activity data
'''
Deterioration_CALC     = 'yes'
Deterioration_list     = ['degradation_rate_Diesel.csv','degradation_rate_Gasoline.csv','degradation_rate_LPG.csv']


'''
Do you want to apply Control factors into emissions?
If YES, set Control_CALC     = 'yes' and insert the control
list file. Remember that vehciles names should match with emissions factor
and Activity data
'''
Control_CALC     = 'no' #'yes'
Control_list     = ['control_factors.csv']


''' IT IS NOT IMPLEMENTED
Do you want to estimate emissions for the future ?
If YES, set future_case     = 'yes' and make sure there are activity data and
emissions factor for the vehicles you want to estimate.
Remember that vehciles names should match with emissions factor and Activity data
'''
future_case      = 'no' #'yes' or 'no'


'''
***** PLOTTING *****
Do you want to create the plots?
*** WARNING ***: The plotting takes a bit longer to finish
'''
plot_figures = 'yes'   #'yes' or 'no'
plot_24 = 'no'        #'yes' or 'no' #24 hours animated plot


'''
'**********      Gridding options      **********'
If the user DO NOT want to generate the grid based on the GRIDDESC,
set the grid_size in meters (e.g 1000) and the grid2CMAQ = 'no'
'''
grid_size = 1000

'''
CARS has the option to create grid to air quality modeling.
To set it, change set  grid2CMAQ to 'yes'
Set the gridfile_name pointing to GRIDDESC file.
Set the Radius_Earth (default=6370000.0)
Set the Datum
'''
grid2CMAQ = 'yes'     #'yes' or 'no'
gridfile_name = met_dir+'/GRIDDESC_NIER_09_01'
Radius_Earth = 6370000.0
Datum = 'WGS84'

'''
Do you want to CARS outpout the GRID shapefile?
If yes, set the outGridShape = 'yes'
This option was added because sometimes generate the shapefile
can take couple minutes.
'''
outGridShape = 'yes'


'''
'*****      Temporal profiles names      *****'
'''

temporal_profile_folder = input_dir+'/temporal_profile'
temporal_monthly_file  = 'monthly_profile.csv'
temporal_week_file     = 'week_profile.csv'
temporal_weekday_file  = 'weekday_profile.csv'
temporal_weekend_file  = 'weekend_profile.csv'
temporal_CrossRef      = 'temporal_profile_CrossRef.csv'


'''
'*****      Chemical speciation for PM2.5, VOC and NOx      *****'
'''
chemical_profile_folder = input_dir+'/chemical_speciation'
pm25_speciation_file  = 'chem_profile_PM2.5.csv'
voc_speciation_file   = 'chem_profile_VOC.csv'
nox_speciation_file   = 'chem_profile_NOx.csv'
speciation_CrossRef   = 'chem_profile_CrossRef.csv'

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
    def __init__(self, pm25_speciation, voc_speciation, nox_speciation, speciation_CrossRef):
        self.PM2_5_spec = pm25_speciation
        self.VOC_spec   = voc_speciation
        self.NOx_spec   = nox_speciation
        self.crossref   = speciation_CrossRef

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
GrayJet = createGrayJet()

        
# =============================================================================
# Function to read Temporal information
# =============================================================================
def read_griddesc(input_dir, gridfile_name, grid2CMAQ = 'no'):
    if ((grid2CMAQ == 'yes') or (grid2CMAQ == 'y') or (grid2CMAQ == 'Y') or \
    (grid2CMAQ == 'YES')) and (gridfile_name != ''):

        start_time = time.time()
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
        
        return GRID_info_table(NTHIK, NCOLS, NROWS, NLAYS, GDTYP, P_ALP, P_BET, P_GAM, 
                               XCENT, YCENT, XORIG, YORIG, XCELL, YCELL, GDNAM)

    else:
        return 

GRID_info = read_griddesc(input_dir, gridfile_name, grid2CMAQ = grid2CMAQ)




# =============================================================================
# Function to read Temporal information
# =============================================================================
def read_temporal_info(Temporal_Profile_Folder, Temporal_Monthly, Temporal_Week,
                       Temporal_WeekDay, Temporal_WeekEnd, Temporal_CrossRef, Run_Period, sep = ';'):
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
#    Temporal_CrossRef  = 'temporal_profile_CrossRef.csv'
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
    
    month_out    = pd.read_csv(month, sep = sep).fillna(np.nan)
    week_out     = pd.read_csv(week, sep = sep).fillna(np.nan)
    weekday_out  = pd.read_csv(weekday, sep = sep).fillna(np.nan)
    weekend_out  = pd.read_csv(weekend, sep = sep).fillna(np.nan)
    crossref_out = pd.read_csv(crossref, sep = sep).fillna(np.nan)

    month_out.columns    = map(str.lower, month_out.columns)
    week_out.columns     = map(str.lower, week_out.columns)
    weekday_out.columns  = map(str.lower, weekday_out.columns)
    weekend_out.columns  = map(str.lower, weekend_out.columns)
    crossref_out.columns = map(str.lower, crossref_out.columns)
    crossref_out['fullname'] = crossref_out.vehicle.str.cat(crossref_out[['types','fuel']], sep='_')

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
#    cols = list(diurnal_temp.columns)
#    if ('weekday') and ('weekend') in list(diurnal_temp.columns):
#        diurnal_temp = diurnal_temp.drop(columns=['weekday','weekend'])
#    elif ('weekday' in cols) and ('weekend' not in cols):
#        diurnal_temp = diurnal_temp.drop(columns=['weekday'])
#    elif ('weekday' not in cols) and ('weekend' in cols):
#    
        
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


TempPro = read_temporal_info(temporal_profile_folder, temporal_monthly_file, temporal_week_file,
                           temporal_weekday_file, temporal_weekend_file, temporal_CrossRef, run_period)



# =============================================================================
# Function to read Chemical speciation
# =============================================================================
def read_chemical_info(Chemical_Speciation_Folder, PM25_Speciation, VOC_Speciation,
                          NOx_Speciation, Speciation_CrossRef, sep = ';'):
    start_time = time.time()
    spec_dir = Chemical_Speciation_Folder
    pm25     = '{0}{1}'.format(spec_dir+'/',PM25_Speciation)
    voc      = '{0}{1}'.format(spec_dir+'/',VOC_Speciation)
    nox      = '{0}{1}'.format(spec_dir+'/',NOx_Speciation)
    crossref = '{0}{1}'.format(spec_dir+'/',Speciation_CrossRef)
    
    files = [pm25, voc, nox, crossref]
    for ifl in files:
        if os.path.exists(ifl) == False:
            print ('')
            print('*** ERROR ABORT ***:  Chemical Speciation file ', ifl, ' "" does not exist!')
            sys.exit('CARS preProcessor can not read Chemical Speciation file')
        else:
            print ('')
            print ('*** Reading Chemical Speciation information ... ***')
            print (ifl)
    
    pm25_out     = pd.read_csv(pm25, sep = sep).fillna(0)
    voc_out      = pd.read_csv(voc, sep = sep).fillna(0)
    nox_out      = pd.read_csv(nox, sep = sep).fillna(0)
    crossref_out = pd.read_csv(crossref, sep = sep).fillna(0)

    pm25_out.columns = map(str.lower, pm25_out.columns)
    voc_out.columns  = map(str.lower, voc_out.columns)
    nox_out.columns  = map(str.lower, nox_out.columns)
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))    
    return Chemical_Speciation_table(pm25_out, voc_out, nox_out, crossref_out)


Chemical_Spec_Table = read_chemical_info(chemical_profile_folder, pm25_speciation_file, voc_speciation_file,
                           nox_speciation_file, speciation_CrossRef)


# =============================================================================
# Function to read the Speed average distribution
# =============================================================================
def read_avgSpeedDistribution(input_dir, avg_Speed_Distribution_file):
    start_time = time.time()
    input_dir = input_dir
    spd_file = avg_SPD_Dist_file #avg_Speed_Distribution_file
    name = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',spd_file)
    if os.path.exists(name) == True:
        print ('')
        print ('Reading Average Speed Distribution table ...')
        print (name)

        with open(name, 'r') as csvfile:
            sep = csv.Sniffer().sniff(csvfile.read(4096)).delimiter

        ASD = pd.read_csv(name, sep = sep).fillna(np.nan)
        
        for ird  in [101, 102, 103, 104, 105, 106, 107, 108]:
            ASD.loc[:,str(ird)] = ASD.loc[:,str(ird)] / ASD.loc[:,str(ird)].sum()
    else:
        print ('')
        print('*** ERROR ABORT ***:  Average speed Distribution ', spd_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read Average speed Distribution')

    out_spd_bins = pd.DataFrame({'spd_bins': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]})
    out_spd      = pd.DataFrame({'spd': list(ASD.Speed)})

    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))
    return EF_Speed_Distribution(ASD, out_spd, out_spd_bins)

avgSpeedDist = read_avgSpeedDistribution(input_dir, avg_SPD_Dist_file)



# =============================================================================
# Function to read the Activity Data
# =============================================================================
def read_activity_data_csv_SK(input_dir, ad_file, End_date, future_case = 'NO'):
    start_time = time.time() 
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
        
    outdf = pd.DataFrame(columns=['fuel','vehicle','types','daily_vkt','region_cd',
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
            Fuel,vehicle,Types,Daily_VKT,Region_code,Manufacture_date,
            '''

            activity_data.columns = ['fuel','vehicle','types','daily_vkt','region_cd',
                                     'manufacture_date']

            activity_data.loc[:,'vehicle'] = activity_data.loc[:,'vehicle'].str.lower()
            activity_data.loc[:,'fuel']    = activity_data.loc[:,'fuel'].str.lower()
            activity_data.loc[:,'types']   = activity_data.loc[:,'types'].str.lower()
            activity_data.loc[:,'manufacture_date']   = (activity_data.loc[:,'manufacture_date'] / 10000).astype(int)
            activity_data['fullname'] = activity_data.vehicle.str.cat(activity_data[['types','fuel']], sep='_')
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
            
            sedan_correction = (activity_data.loc[(activity_data.vehicle == 'sedan'),'daily_vkt']) * 3
            activity_data.loc[(activity_data.vehicle == 'sedan'),'daily_vkt'] = sedan_correction
            suv_correction = (activity_data.loc[(activity_data.vehicle == 'suv'), 'daily_vkt']) * 3
            activity_data.loc[(activity_data.vehicle == 'suv'), 'daily_vkt'] = suv_correction
            
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


    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))
    return Activity_Data_table(outdf, vhc_names, vhc_years, vhc_fuels, vhc_cars, count_outdf)
## =============================================================================

AD_SK = read_activity_data_csv_SK(input_dir, activity_file, ENDATE)



# =============================================================================
# Function to read link level shapefile
# =============================================================================
def roads_grid_surrogate_inf(input_dir, File_Name, Link_Shape_att, GridDesc_file='', Radius_Earth='', outGridShape = 'yes', Unit_meters = True):
    start_time = time.time()
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
    if ((grid2CMAQ == 'yes') or (grid2CMAQ == 'y') or (grid2CMAQ == 'Y') or \
        (grid2CMAQ == 'YES')) and (gridfile_name == ''):
        print('')
        print('*** ERROR ***: Error on processing grid to CMAQ')
        print('*** ERROR ***: Please set the gridfile_name and grid2CMAQ the correctly')
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
        rt = {101 : 80 *reduc, 102 : 60 *reduc, 103 : 60 *reduc, 104 : 50 *reduc,   #60% reduction as BH asked
              105 : 30 *reduc, 106 : 30 *reduc, 107 : 30 *reduc, 108 : 30 *reduc}

        for key, ispd in rt.items():
            out_roads.loc[(out_roads.loc[:,'road_type'] == key),'max_speed'] = ispd
        # Creating grid
        if ((grid2CMAQ == 'yes') or (grid2CMAQ == 'y') or (grid2CMAQ == 'Y') or \
            (grid2CMAQ == 'YES')) and (gridfile_name != ''):
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
            if (outGridShape == 'yes') or (outGridShape == 'YES') or \
               (outGridShape == 'y')   or (outGridShape == 'Y'):
                grid.to_file(filename = output_dir+'/grid_{0}.shp'.format(case_name), driver='ESRI Shapefile',crs_wkt=prj)

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
        
        surrogate = pd.merge(grid, surrogate, how='left', on=['grid_id']).fillna(0)
        out_roads = out_roads.drop(columns=['geometry'])
        out_roads = out_roads.rename(columns={'geometry_BKP': 'geometry'}).set_geometry('geometry')
        out_roads['region_cd'] = out_roads['region_cd'].astype(int)
        surrogate['region_cd'] = surrogate['region_cd'].astype(int)
    else:
        print('*** ERROR ABORT ***:  Shapefile "" ', shp_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read link Shapefile file')


    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))

    return Roads_Grid_table(grid, surrogate, out_roads)
# =============================================================================
    
roads_RGS = roads_grid_surrogate_inf(input_dir, link_shape, link_shape_att,
                                     gridfile_name, Radius_Earth, Unit_meters = True ) 

# =============================================================================
# Function to read link level shapefile
# =============================================================================
def processing_County_shape(input_dir, file_name, County_Shape_att, GridDesc_file='', Radius_Earth=''):
    start_time = time.time()
    Region_Geocode       = County_Shape_att[0] #Region_CD
    Region_name_attr     = County_Shape_att[1] #Region_name_attr 
    Region_name_attr_SK  = County_Shape_att[2] #Region_name_attr_SK
#    Link_ID_attr       = 'LINK_ID'
#    RD_name_attr       = 'ROAD_NAME'
#    RD_type_attr       = 'ROAD_RANK'
#    Activity_data_attr = 'SHAPE_STLe'
#    Speed_attr         = 'MAX_SPD'
#    Link_length        = 'SHAPE_STLe'
#    file_name          = link_shape
    if Radius_Earth != '':
        radius = Radius_Earth
    else:
        radius = 6370000.0

    if ((grid2CMAQ == 'yes') or (grid2CMAQ == 'y') or (grid2CMAQ == 'Y') or \
        (grid2CMAQ == 'YES')) and (gridfile_name != ''):
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
                                           'region_name_SK']]
                
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
                                   'region_name_SK']]
        
        out_cnty['region_cd'] = out_cnty['region_cd'].astype(int)

    else:
        print('*** ERROR ABORT ***:  Shapefile "" ', cnty_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read county Shapefile file')         
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))

    return out_cnty
# =============================================================================
    
county_SHP = processing_County_shape(input_dir, county_shape, county_shape_att, gridfile_name, Radius_Earth,) 



# =============================================================================
# Function to read the emissions factor from South Korea* (*particular case)
# =============================================================================
def read_emissions_factor_SK(input_dir, EmisFactor_list):
    start_time = time.time()
    input_dir = input_dir
    ef_list = EmisFactor_list #ef_file #['gasoline.csv'] #
#    ef_list = ['diesel_v3.csv']
#    sep = ';'
    final_EF = pd.DataFrame()
    for ifile in range(0,len(ef_list)):
#        print(ef_list[ifile])
        name = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',ef_list[ifile])
        if os.path.exists(name) == True:
            print('******************************************************************')
            print('*****         Reading Emissions Factors tables ...           *****')
            print('*****                    Please wait ...                     *****')
            print('******************************************************************')
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
                aux = 0 in list(emissions_factor.loc[pos,['vehicle','types','scc','pollutant']])
                if aux == True:
                    drop_list.append(pos)
                iaux +=1
            emissions_factor = emissions_factor.drop(index=drop_list)
            emissions_factor['years'] =  emissions_factor['years'] / 10000
            emissions_factor['years'] =  emissions_factor.years.astype(int)
            emissions_factor.loc[:,'vehicle'] = emissions_factor.loc[:,'vehicle'].str.lower()
            emissions_factor.loc[:,'types'] = emissions_factor.loc[:,'types'].str.lower()
            emissions_factor.loc[:,'fuel'] = emissions_factor.loc[:,'fuel'].str.lower()
            emissions_factor['fullname'] = emissions_factor.vehicle.str.cat(emissions_factor[['types','fuel']], sep='_')
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
    
    return EF_Grid_table(final_EF, EF_years, EF_names, EF_fuels, EF_polls)
# =============================================================================

EF_All_fuel_function = read_emissions_factor_SK(input_dir, Emis_Factor_list)




# =============================================================================
# Function to calculate the emissions factor from South Korea
# =============================================================================                              
def calculate_EmisFact(Input_dir, Emissions_Factor_Table, Out_dir ):
    start_time = time.time()
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
    dat_EF['spd']   = dat_EF['spd'].astype(int)
    dat_EF['years'] = dat_EF['years'].astype(int)
    
    #Extracting the NOx Diesel EF from final EF table
    Temp_EF_Diesel = dat_EF[(dat_EF.temperatures != 0) & 
                            (dat_EF.pollutant == 'NOx')]
    Temp_EF_Diesel_drop_list = list(Temp_EF_Diesel.index)
    dat_EF = dat_EF.drop(index=Temp_EF_Diesel_drop_list).reset_index(drop=True)
    Temp_EF_Diesel = Temp_EF_Diesel.reset_index(drop=True)
    
    #Filtering the NOx Diesel EF based on the ambient temperature
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
    
    #Filtering the duplicate EF of NOx Diesel - we are averaging them
    Temp_EF_Diesel = dat_EF[(dat_EF.temperatures != 0) & 
                            (dat_EF.pollutant == 'NOx')]
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
    
    EmisFact_yr_spd = EF_Grid_table(dat_EF, EF_years, EF_names, EF_fuels, EF_polls)
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))
    
    return EmisFact_yr_spd

EmisFactor_yr_spd = calculate_EmisFact(input_dir, EF_All_fuel_function, output_dir)

# =============================================================================
# Function to check any missing emissions factor
# =============================================================================
def check_VHC_AD_EF(Emissions_Factor_DataFrame, Activity_Data_DataFrame, Out_dir):
    EFdf = Emissions_Factor_DataFrame
    ADdf = Activity_Data_DataFrame
    OutD = Out_dir
    EFvhc = EFdf.EF_fullname.EF_fullname.tolist()
    ADvhc = ADdf.fullname.vhc_name.tolist()
    print('')
    print('******************************************************************')
    print('***** Checking if vehicle if all Vehicles have Emission Factor ***')
    print('*****            Please wait and check the LOG...            *****')
    print('******************************************************************')
    print('')
    log_name = 'LOG_Missing_Vehicle_EmissionFactor.txt'
    with open(OutD+'/LOGS/'+log_name, 'w') as log:
        log.write('# ***** Vehicles without Emissions Factors *****\n')
        log.write('# \n')
        log.write('vehicle_types_fuel\n')
        for ivhc in ADvhc:
            if ivhc not in EFvhc:
                log.write('{0}\n'.format(ivhc))
    return       


try:
    AD_SK
    a = True
except: a = False

try:
    EmisFactor_yr_spd
    b = True
except: b = False

if (a == True) and (b == True):
    check_VHC_AD_EF(EmisFactor_yr_spd, AD_SK, output_dir)

# =============================================================================
# Function to apply Deterioration rate into the emissions factor from South Korea
# =============================================================================                              

if (Deterioration_CALC == 'yes') or (Deterioration_CALC == 'YES') or (Deterioration_CALC == 'Yes') or (Deterioration_CALC == 'Y') or (Deterioration_CALC == 'y'):
    print('') 
    print('******************************************************************')
    print('*****         Deterioration_CALC set to YES ...              *****')
    print('*****        Applying Deterioration factor to EF             *****')
    print('******************************************************************')
    print('')
    print('***** WARNING: The pure emissions factor will be discarded and only \
          the deteriorated EF will be used, unless there is no deterioration for an \
          specific vehicle or/and pollutant *****')    
    def apply_deterioration_emissions_factor_SK(input_dir, Deterioration_list, EmisFactor_yr_spd, output_dir):
        start_time = time.time()
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
                    aux = 0 in list(det_factor.loc[pos,['vehicle','types','scc','pollutant']])
                    if aux == True:
                        drop_list.append(pos)
                    iaux +=1
                det_factor = det_factor.drop(index=drop_list)
                det_factor['manufacture_date'] =  det_factor['manufacture_date'] / 10000
                det_factor['manufacture_date'] =  det_factor['manufacture_date'].astype(int)
                det_factor.loc[:,'vehicle'] = det_factor.loc[:,'vehicle'].str.lower()
                det_factor.loc[:,'types'] = det_factor.loc[:,'types'].str.lower()
                det_factor.loc[:,'fuel'] = det_factor.loc[:,'fuel'].str.lower()
                det_factor['fullname'] = det_factor.vehicle.str.cat(det_factor[['types','fuel']], sep='_')
                final_DEF = final_DEF.append(det_factor, ignore_index=True, sort=False)
        
        final_DEF = final_DEF.sort_values(by=['manufacture_date','fullname']).reset_index(drop=True)
        final_DEF = final_DEF.set_index(['fullname','pollutant','manufacture_date'])
        final_DEF = final_DEF.drop(columns=['vehicle','types','fuel','scc'])
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
        EF_fuels = pd.DataFrame({'EF_fuel'    : list(EmisFactor_yr_spd.EF_fuels.EF_fuel)})
        EF_polls = pd.DataFrame({'EF_poll'    : list(np.sort(merge_EF_Final.pollutant.unique()))})
        
        merge_EF_Final.to_csv(output_dir+'/'+'EmisFact_by_YR_SPD_Deteriorated.csv',sep=',', index=False)
        
        EmisFact_yr_spd = EF_Grid_table(merge_EF_Final, EF_years, EF_names, EF_fuels, EF_polls)
        run_time = ((time.time() - start_time))
        print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
        print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
        print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))
        
        return EmisFact_yr_spd

    EmisFactor_yr_spd = apply_deterioration_emissions_factor_SK(input_dir, Deterioration_list, EmisFactor_yr_spd, output_dir)
else:
    print('') 
    print('******************************************************************')
    print('*****         Deterioration_CALC set to NO ...               *****')
    print('*****        Deterioration factor will not be applied        *****')
    print('*****        Using the standad emissions factors             *****')
    print('******************************************************************')
    print('')

      
# =============================================================================
# Function to calculate the county emissions 
# =============================================================================
def calc_County_Emissions(Input_dir, Emissions_Factor_DataFrame, 
                          Activity_Data_DataFrame, Roads_DataFrame, County_DataFrame,
                          Cold_Start_input, Out_dir ):
    print('')
    print('******************************************************************')
    print('*****           Starting calculation of emissions            *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')
    start_time = time.time()
    input_dir  = Input_dir
    output_dir = Out_dir                     #
    EF_yr_spd  = Emissions_Factor_DataFrame  #
    AD_yr_vhc  = Activity_Data_DataFrame     #
    roads_DF   = Roads_DataFrame             #
    county_df  = County_DataFrame            # 
    cold_start = Cold_Start_input

#    input_dir  = input_dir         
#    output_dir = output_dir        
#    EF_yr_spd  = EmisFactor_yr_spd 
#    AD_yr_vhc  = AD_SK             
#    roads_DF   = roads_RGS         
#    county_df  = county_SHP        
#    cold_start = Cold_Start_list  
    
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
        log.write('vehicle_types_fuel\n')
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
    

    AD_byRoad = pd.merge(final_AD_df,road_weight, left_on = ['region_cd'], 
                        right_on = ['region_cd'], how='left').fillna(0.0)

    AD_byRoad.loc[:,ADvhc] = AD_byRoad.loc[:,ADvhc].apply(lambda x: np.asarray(x) * AD_byRoad.vkt_split_county.values)
    ISzeroRD = AD_byRoad['road_type'].loc[AD_byRoad['road_type'] == 0].index 
    AD_byRoad['road_type'].iloc[ISzeroRD] = 105
    AD_byRoad['road_type'] = AD_byRoad['road_type'].astype(int)
    
    road_by_county = AD_byRoad.groupby(['region_cd','road_type']).sum().reset_index(drop=False)
    road_by_county = road_by_county.drop(columns=['manufacture_date', 'vkt_split_county'])
    road_by_county = road_by_county.sort_values(by=['region_cd','road_type'])
    aux_RbC        = road_by_county.groupby(['region_cd']).sum().reset_index(drop=False)
    aux_RbC        = aux_RbC.drop(columns=['road_type'])
    aux_RbC        = pd.merge(aux_RbC, road_by_county.loc[:,['region_cd','road_type']],
                              left_on = 'region_cd', right_on='region_cd', how='left')
    road_by_county = road_by_county.set_index(['region_cd','road_type'])
    aux_RbC = aux_RbC.set_index(['region_cd','road_type'])
    road_by_county.loc[:,ADvhc] = road_by_county.loc[:,ADvhc].div(aux_RbC.loc[:,ADvhc], axis=0)
    road_by_county = road_by_county.reset_index(drop=False).fillna(1.0)
    
    EF_avgSpd = pd.DataFrame()
    roads = ['101','102','103','104','105','106','107','108']
    for iroad in roads:
        aux_avgSpd = avgSpeedDist.data.loc[:,['Speed',str(iroad)]]
        
        aux_EFavgSpd = pd.merge(final_EF_df,aux_avgSpd, left_on = 'spd', 
                                  right_on = 'Speed', how='right').reset_index(drop=True)
        
        aux_EFavgSpd.loc[:,EFvhc] = aux_EFavgSpd.loc[:,EFvhc].apply(lambda x: np.asarray(x) * aux_EFavgSpd[iroad].values)
        aux_EFavgSpd = aux_EFavgSpd.groupby(['years','pollutant']).sum().reset_index(drop=False)
        aux_EFavgSpd = aux_EFavgSpd.drop(columns=['spd', 'Speed', iroad])
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
    Emissions_by_yr_county2csv = Emissions_by_yr_county2csv.drop(columns=['manufacture_date','region_cd','road_type'])
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
            cold_start_df.loc[:,'types'] = cold_start_df.loc[:,'types'].str.lower()
            cold_start_df.loc[:,'fuel'] = cold_start_df.loc[:,'fuel'].str.lower()
            cold_start_df['fullname'] = cold_start_df.vehicle.str.cat(cold_start_df[['types','fuel']], sep='_')


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
                                                Emissions_by_yr_county.pollutant.isin(['CO','VOC','NOx','PM2.5'])].reset_index(drop=True)
    
    cold_start_aux = cold_start_aux.loc[:,['manufacture_date','region_cd','road_type','pollutant']+CS_vhc_list]

    beta = 0.647 - (0.025 * 12.35) - (0.00974 - 0.000385 * 12.35) * temp_mean
    
    EcEh_Gasoline_CO  =  beta * ((9.040 - 0.090 * temp_mean) - 1)
    EcEh_Gasoline_VOC =  beta * ((12.59 - 0.060 * temp_mean) - 1)
    EcEh_Gasoline_NOx =  beta * ((3.660 - 0.006 * temp_mean) - 1)

    EcEh_Diesel_CO    =  beta * ((1.900 - 0.030 * temp_mean) - 1)
    if temp_mean > 29:
        EcEh_Diesel_VOC   =  beta * 0.5
    else:
        EcEh_Diesel_VOC   =  beta * ((3.100 - 0.090 * temp_mean) - 1)

    EcEh_Diesel_NOx   =  beta * ((1.300 - 0.013 * temp_mean) - 1)

    if temp_mean > 26:
        EcEh_Diesel_PM    =  beta * 0.5
    else:
        EcEh_Diesel_PM    =  beta * ((3.100 - 0.100 * temp_mean) - 1)


    EcEh_LPG_CO       =  beta * ((3.660 - 0.090 * temp_mean) - 1)

    if temp_mean > 29:
        EcEh_LPG_VOC      =  beta * 0.5
    else:
        EcEh_LPG_VOC      =  beta * ((2.240 - 0.060 * temp_mean) - 1)
        
    EcEh_LPG_NOx      =  beta * ((0.980 - 0.006 * temp_mean) - 1)
    
    aux_dict = {'pollutant' : ['CO','VOC','NOx','PM2.5'],
                'gasoline'  : [EcEh_Gasoline_CO, EcEh_Gasoline_VOC, EcEh_Gasoline_NOx, 0.0],
                'diesel'    : [EcEh_Diesel_CO  , EcEh_Diesel_VOC  , EcEh_Diesel_NOx  , EcEh_Diesel_PM],
                'lpg'       : [EcEh_LPG_CO     , EcEh_LPG_VOC     , EcEh_LPG_NOx     , 0.0]}
    
    EcEhBeta_df = pd.DataFrame(aux_dict)
    EcEhBeta_df.loc[:,ColdStart_fuels] = EcEhBeta_df.loc[:,ColdStart_fuels].clip(lower=0)

    cold_start_emissions = pd.DataFrame()
    for ifuel in ColdStart_fuels:
        aux_CS_emis = pd.DataFrame()
        for ipol in ['CO','VOC','NOx','PM2.5']:
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
    print('*****      Calculating the Evaporatiuve Emissions ...        *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    '''
    Calculating the Evaporatiuve Emissions
    beta = 0.647 - (0.025 * 1Trip) - (0.00974 - 0.000385 * 1Trip)* temp_mean
    1Trip = 12.35km - This value came from South Korea survey and it means that
    each vehicle run, on average, 12.35 km per day
    '''   
    RVP    = 82.0
    beta   = 0.647 - (0.025 * 12.35) - (0.00974 - 0.000385 * 12.35) * temp_mean
    ErHot  = 0.136 * (np.exp(-5.967 + 0.04259 * RVP + 0.1773 * temp_mean ))
    ErWarm = 0.100 * (np.exp(-5.967 + 0.04259 * RVP + 0.1773 * temp_mean ))
    Ed_aux = 9.100 * (np.exp(0.0158*(RVP-61.2) + 0.0574 *(temp_min - 22.5) + 0.0614 * \
                     (temp_max - temp_min - 11.7)))
    

    aux_EE_AD = AD_byRoad.copy()
    list_vhc_no_gas = list(aux_EE_AD.loc[0,~aux_EE_AD.columns.str.contains('gasoline')].index)
    list_vhc_no_gas.remove('manufacture_date')
    list_vhc_no_gas.remove('region_cd')
    list_vhc_no_gas.remove('road_type')
    list_vhc_no_gas.remove('vkt_split_county')
    aux_EE_AD.loc[:,list_vhc_no_gas] = aux_EE_AD.loc[:,list_vhc_no_gas] * 0.0
    
    # Ed emissions is based on the factor Ed_aux, so it is just multiply the vhc count by the Ed_aux
    Ed = AD_yr_vhc.vhc_count.copy()
    Ed.loc[:,list_vhc_no_gas] = Ed.loc[:,list_vhc_no_gas] * 0.0
    Ed.loc[:,ADvhc] = Ed.loc[:,ADvhc].mul(Ed_aux)
    Ed['region_cd'] = Ed['region_cd'].astype(int)
    Ed = pd.merge(aux_EE_AD.loc[:,['manufacture_date','region_cd','road_type','vkt_split_county']],
                  Ed,
                  left_on  = ['manufacture_date','region_cd'],
                  right_on = ['manufacture_date','region_cd'], how='left')
    Ed.loc[:,ADvhc] = Ed.loc[:,ADvhc].apply(lambda x: np.asarray(x) * \
                      Ed.vkt_split_county.values).fillna(0.0)
    
    #Calculating R based on South Korea equation
    R_aux = ((1 - beta) * ErHot) + (beta * ErWarm)
    R     = aux_EE_AD.copy()
    R.loc[:,ADvhc] = R.loc[:,ADvhc].mul(R_aux)
    R.set_index(['manufacture_date','region_cd','road_type'], inplace=True)
    R = R.sort_index()
    R = R.reset_index(drop=False)
   
    #Calculating Sfi based on South Korea equation 
    Sfi = aux_EE_AD.copy()
    Sfi.loc[:,ADvhc] = Sfi.loc[:,ADvhc].mul(0.7 / 12.35)
    Sfi.set_index(['manufacture_date','region_cd','road_type'], inplace=True)
    Sfi = Sfi.sort_index()
    Sfi = Sfi.reset_index(drop=False)

    Ed  = Ed.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True) 
    R   = R.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True) 
    Sfi = Sfi.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True) 
    

    #Copying Ed emissions to the final dataframe
    evaporative_emissions = Ed.copy()
    #Adding R emissions to final dataframe
    evaporative_emissions.loc[:,ADvhc] = evaporative_emissions.loc[:,ADvhc].add \
                                                (R.loc[:,ADvhc], axis=1)
    #Adding Sfi emissions to final dataframe
    evaporative_emissions.loc[:,ADvhc] = evaporative_emissions.loc[:,ADvhc].add \
                                                (Sfi.loc[:,ADvhc], axis=1)
        
    evaporative_emissions['pollutant'] = 'VOC'
    evaporative_emissions = evaporative_emissions.loc[:,
    ['manufacture_date','region_cd','pollutant','road_type','vkt_split_county'] + ADvhc]

    cols_2merge = ['manufacture_date', 'region_cd', 'road_type', 'pollutant']
    merge_EE_HE = pd.merge(Emissions_by_yr_county.loc[:,cols_2merge],evaporative_emissions,
                   left_on = cols_2merge,
                   right_on = cols_2merge,
                   how='left').fillna(0.0)
    
    Emissions_by_yr_county.loc[:,ADvhc] = Emissions_by_yr_county.loc[:,ADvhc].add \
                                                (merge_EE_HE.loc[:,ADvhc], axis=1)
   
    #converting Evaporative emissions to Tons/yr to output it
    evaporative_emissions.loc[:,ADvhc] = evaporative_emissions.loc[:,ADvhc] * (ndaysYEAR/1000000)
    evaporative_emissions.to_csv(output_dir+'/evaporative_emissions_Tons_per_Year.csv', sep=',')


    total_county = Emissions_by_yr_county.groupby(['region_cd','pollutant']).sum().reset_index(drop=False)
    total_county = total_county.drop(columns=['manufacture_date','road_type'])


    #Calculation the total emissions for county and link level
    Emissions_by_yr_county = Emissions_by_yr_county.groupby(['region_cd','pollutant','manufacture_date']).sum().reset_index(drop=False)
    Emissions_by_yr_county['total_emis'] = Emissions_by_yr_county.loc[:,ADvhc].sum(axis=1)
    total_county['total_emis'] = total_county.loc[:,ADvhc].sum(axis=1)
    county_table = pd.DataFrame({'region_cd' : Emissions_by_yr_county.region_cd.unique()})
    years_table  = pd.DataFrame({'years' : Emissions_by_yr_county.manufacture_date.unique()})
    vhc_table    = pd.DataFrame({'fullname' : ADvhc})
    
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
    total_county_WGT = total_county_WGT.loc[:,['region_cd','pollutant']+ADvhc]
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

    #Exporting the CSV files with the emissions 
    Emis_by_yr_county_2CVS = Emissions_by_yr_county.groupby(['region_cd', 'pollutant', 'manufacture_date']).sum().reset_index(drop=False)
    Emis_by_yr_county_2CVS = Emis_by_yr_county_2CVS.drop(columns=['road_type','region_state'])
    Emis_by_yr_county_2CVS.to_csv(output_dir+'/'+'County_Total_Emissions_by_Yearly_Activity_Data_Year.csv', sep=',', index=False)
    Emis_by_yr_county_2CVS = []
    
    Emis_by_yr_county_WGT_2CVS = Emissions_by_yr_county_WGT.groupby(['region_cd', 'pollutant', 'manufacture_date']).sum().reset_index(drop=False)
#    Emis_by_yr_county_WGT_2CVS = Emis_by_yr_county_WGT_2CVS.drop(columns=['road_type','region_state'])
    Emis_by_yr_county_WGT_2CVS.to_csv(output_dir+'/'+'County_Total_Normalized_Emissions_by_Yearly_Activity_Data_Year.csv', sep=',', index=False)
    Emis_by_yr_county_WGT_2CVS = []
    
    total_poll = total_county.groupby(['pollutant']).sum().reset_index(drop=False)
    total_poll = total_poll.drop(columns=['region_cd'])
    total_poll.to_csv(output_dir+'/'+'Pollutant_Total_Emissions_Tons_per_Year.csv', sep=',', index=False)
    total_poll = []
    
    total_poll_WGT = total_county_WGT.groupby(['pollutant']).sum().reset_index(drop=False)
    total_poll_WGT = total_poll_WGT.drop(columns=['region_cd'])
    total_poll_WGT.to_csv(output_dir+'/'+'Pollutant_Total_Normalized_Emissions_Tons_per_Year.csv', sep=',', index=False)
    total_poll_WGT = []
    
    total_county_2CVS = total_county.drop(columns = ['geometry','region_state'])
    total_county_2CVS.to_csv(output_dir+'/'+'County_Total_Emissions_Tons_per_Year.csv', sep=',', index=False)
    total_county_2CVS = []
    
    total_county_WGT_2CVS = total_county_WGT.drop(columns = ['geometry'])
    total_county_WGT_2CVS.to_csv(output_dir+'/'+'County_Total_Normalized_Emissions_Tons_per_Year.csv', sep=',', index=False)
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
       
    return Emissions_table(Emissions_by_yr_county_WGT, total_county, total_county_WGT, county_table, years_table, vhc_table, road_by_county)

County_Emissions = calc_County_Emissions(input_dir, EmisFactor_yr_spd,
                                       AD_SK, roads_RGS, county_SHP, Cold_Start_list, output_dir)




def aplly_control(Control_list_file, County_Emissions_df, County_SHP_df):
    start_time = time.time()
    control_list = Control_list      #Control_list_file    #
    cnt_df       = County_Emissions  #County_Emissions_df  #
    county_df    = county_SHP        #County_SHP_df        #
    ADvhc        = cnt_df.fullname.fullname.to_list()
    for ifile in range(0,len(Control_list)):
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
            control_df.loc[:,'types'] = control_df.loc[:,'types'].str.lower()
            control_df.loc[:,'fuel'] = control_df.loc[:,'fuel'].str.lower()
            control_df['fullname'] = control_df.vehicle.str.cat(control_df[['types','fuel']], sep='_')


    control_vhc = list(control_df.fullname.unique())
    control_region = list(control_df.region_cd.unique())
    Controllog_name = 'LOG_Control_factor_checking.txt'
    with open(output_dir+'/LOGS/'+Controllog_name, 'w') as CF_log:
        CF_log.write('*****   List of vehicle with more than one control factor per year    *****\n')
        
        duplicate = control_df.drop(columns = 'control_factor')
        aux_dup = duplicate.loc[control_df.duplicated()]
        if aux_dup.shape[0] != 0:
            print('')
            print('***** ERROR *****: Please check the LOG_Control_factor_checking.txt file')
            print('***** ERROR *****: There are vehicles with more than one control factor per year')
            print('')
            for i in aux_dup.index:
                CF_log.write('Vehicle {0} and region {1} has more than one control factor per year and/or per region.\n'.format(duplicate.loc[i,'fullname'], duplicate.loc[i,'region_cd']))
            sys.exit()
        
        duplicate.loc[:,'region_cd'] = duplicate.loc[:,'region_cd'].astype(str)
        duplicate.loc[:,'region_cd'] = duplicate.loc[:,'region_cd'].str.slice(stop=2).astype(int)
        duplicate.loc[:,'region_cd'] = duplicate.loc[:,'region_cd'].astype(int)
        aux_dup = duplicate.loc[duplicate.duplicated()]
        if aux_dup.shape[0] != 0:
            print('')
            print('***** ERROR *****: Please check the LOG_Control_factor_checking.txt file')
            print('***** ERROR *****: There are vehicles with more than one control factor per year')
            print('')
            for i in aux_dup.index:
                CF_log.write('Vehicle {0} and region {1} has more than one control factor per year and/or per region.\n'.format(duplicate.loc[i,'fullname'], duplicate.loc[i,'region_cd']))
            sys.exit()

    control_df.loc[:,'control_factor'] = 1 - (control_df.loc[:,'control_factor'] / 100)
    control_df.loc[:,'region_cd']      = control_df.loc[:,'region_cd'].astype(str)
    control_df['region_state']         = control_df.loc[:,'region_cd'].str.slice(stop=2).astype(int)
    control_df.loc[:,'region_cd']      = control_df.loc[:,'region_cd'].astype(int)
    
    cnt_df.county_by_yr.loc[:,'region_cd']      = cnt_df.county_by_yr.loc[:,'region_cd'].astype(str)
    cnt_df.county_by_yr['region_state']         = cnt_df.county_by_yr.loc[:,'region_cd'].str.slice(stop=2).astype(int)
    cnt_df.county_by_yr.loc[:,'region_cd']      = cnt_df.county_by_yr.loc[:,'region_cd'].astype(int)
    
    drop_all_pol_list = control_df.loc[control_df.data == 'ALL'].index
    control_df_all = control_df.loc[control_df.data == 'ALL'].reset_index(drop=True)
    control_df_all = control_df_all.drop(columns=['vehicle','types','fuel'])
    control_df = control_df.drop(index=drop_all_pol_list).reset_index(drop=True)
    
    apply_factor_sta = cnt_df.county_by_yr.loc[:,['region_cd','region_state','pollutant', 'manufacture_date']]
    apply_factor_cnt = cnt_df.county_by_yr.loc[:,['region_cd','region_state','pollutant', 'manufacture_date']]
    
    sta_vhc = []
    cnt_vhc = []
    
    control_df.loc[:,'region_cd']      = control_df.loc[:,'region_cd'].astype(int)
    for irgn in control_df.region_cd.unique():
        aux_irg = control_df.loc[(control_df.region_cd == irgn)]
        appd_CF = pd.DataFrame()
        for ipol in aux_irg.data.unique():
            aux_CF = aux_irg.loc[(aux_irg.data == ipol)]
            aux_CF = aux_CF.pivot(index='year', columns='fullname', values='control_factor').fillna(1).reset_index(drop=False)
            aux_CF['data']      = ipol
            aux_CF['region_cd'] = irgn
            aux_CF = aux_CF.rename(columns={'region_cd' : 'region_general'})
            appd_CF = appd_CF.append(aux_CF, ignore_index=True, sort=False)
        if (len(str(irgn)) == 2) or (len(str(irgn)) == 1):
            apply_factor_sta = pd.merge(apply_factor_sta, appd_CF, 
                                        left_on=['region_state','pollutant', 'manufacture_date'],
                                        right_on=['region_general','data','year'],
                                        how='left')
            apply_factor_sta = apply_factor_sta.drop(columns=['region_general','data','year'])
            sta_vhc.extend(list(aux_irg.fullname.unique()))
        elif (len(str(irgn)) == 8) or (len(str(irgn)) == 7):
            apply_factor_cnt = pd.merge(apply_factor_cnt, appd_CF, 
                                        left_on=['region_cd','pollutant', 'manufacture_date'],
                                        right_on=['region_general','data','year'],how='left')
            apply_factor_cnt = apply_factor_cnt.drop(columns=['region_general','data','year'])
            cnt_vhc.extend(list(aux_irg.fullname.unique()))
        else:
            print('*****ERROR: Control factor is only applied in State level (2 digits state code e.g: 11) or')
            print('county level (8 digitis county code e.g: 11110129)')
            sys.exit()

    for irgn in control_df_all.region_cd.unique():
        if (len(str(irgn)) == 2) or (len(str(irgn)) == 1):
            aux_all_cf = control_df_all.loc[control_df_all.region_cd == irgn]
            sta_vhc.extend(list(aux_all_cf.fullname.unique()))
            aux_all_cf = aux_all_cf.pivot(index='year', columns='fullname', values='control_factor').fillna(1).reset_index(drop=False)
            aux_all_cf['region_cd'] = irgn
            aux_all_cf = aux_all_cf.rename(columns={'region_cd' : 'region_general'})
            apply_factor_sta = pd.merge(apply_factor_sta,aux_all_cf, 
                                        left_on  =['region_state', 'manufacture_date'],
                                        right_on =['region_general', 'year'], how='left')
            apply_factor_sta = apply_factor_sta.drop(columns=['region_general','data','year'])
            
        elif (len(str(irgn)) == 8) or (len(str(irgn)) == 7):
            aux_all_cf = control_df_all.loc[control_df_all.region_cd == irgn]
            cnt_vhc.extend(list(aux_all_cf.fullname.unique()))
            aux_all_cf = aux_all_cf.pivot(index='year', columns='fullname', values='control_factor').fillna(1).reset_index(drop=False)
            aux_all_cf['region_cd'] = irgn
            aux_all_cf = aux_all_cf.rename(columns={'region_cd' : 'region_general'})

            apply_factor_cnt = pd.merge(apply_factor_cnt,aux_all_cf, 
                                        left_on  =['region_cd', 'manufacture_date'],
                                        right_on =['region_general', 'year'], how='left')

            apply_factor_cnt = apply_factor_cnt.drop(columns=['region_general','year'])
            
        else:
            print('*****ERROR: Control factor is only applied in State level (2 digits state code e.g: 11) or')
            print('county level (8 digitis county code e.g: 11110129)')
            sys.exit()

    apply_factor_sta = apply_factor_sta.fillna(1.0)
    apply_factor_cnt = apply_factor_cnt.fillna(1.0)
    
    for ivhc in sta_vhc:
        cnt_df.county_by_yr.loc[:,ivhc] = cnt_df.county_by_yr.loc[:,ivhc] * \
                                          apply_factor_sta.loc[:,ivhc]
    
    for ivhc in cnt_vhc:
        cnt_df.county_by_yr.loc[:,ivhc] = cnt_df.county_by_yr.loc[:,ivhc] * \
                                          apply_factor_cnt.loc[:,ivhc]
                    
    cnt_df.county_by_yr['total_emis'] = cnt_df.county_by_yr.loc[:,ADvhc].sum(axis=1)
    cnt_df.county_by_yr = cnt_df.county_by_yr.sort_values(by=['region_cd','pollutant'])
    
    cnt_df.county_total = cnt_df.county_by_yr.groupby(['region_cd','pollutant']).sum().reset_index(drop=False)
    cnt_df.county_total = cnt_df.county_total.drop(columns=['manufacture_date','road_type'])
    cnt_df.county_total = cnt_df.county_total.sort_values(by=['pollutant','region_cd']).reset_index(drop=True)
    
    aux_pol = cnt_df.county_total.pollutant.unique()[0]
    aux_wgt = cnt_df.county_total_WGT.loc[(cnt_df.county_total_WGT.pollutant == aux_pol),
                                          ['region_cd','geometry','region_nm','vkt_weight']].reset_index(drop=True)

    cnt_df.county_total = pd.merge(cnt_df.county_total, aux_wgt, 
                                   left_on = ['region_cd'],
                                   right_on = ['region_cd'], how='left')
    
    #inserting the geometry and transforming it to geodataframe
    cnt_df.county_total = gpd.GeoDataFrame(cnt_df.county_total, crs=county_df.crs, geometry='geometry')
    
    cnt_df.county_total_WGT = cnt_df.county_total.copy()
    cnt_df.county_total_WGT.loc[:,ADvhc] = cnt_df.county_total_WGT.loc[:,ADvhc].apply(lambda x: np.asarray(x) * \
                               cnt_df.county_total_WGT.vkt_weight.values)
    
    #inserting the geometry and transforming it to geodataframe
    cnt_df.county_total_WGT = gpd.GeoDataFrame(cnt_df.county_total_WGT, crs=county_df.crs, geometry='geometry')

    total_county_YR_2CVS = cnt_df.county_by_yr.groupby(['region_cd', 'pollutant', 'manufacture_date']).sum().reset_index(drop=False)
    total_county_YR_2CVS = total_county_YR_2CVS.drop(columns = ['road_type'])
    total_county_YR_2CVS.to_csv(output_dir+'/'+'County_Total_by_Year_Controled_Emissions_Tons_per_Year.csv', sep=',', index=False)

    total_county_2CVS = cnt_df.county_by_yr.groupby(['region_cd', 'pollutant']).sum().reset_index(drop=False)
    total_county_2CVS = total_county_2CVS.drop(columns = ['region_state','manufacture_date','road_type'])
    total_county_2CVS.to_csv(output_dir+'/'+'County_Total_Controled_Emissions_Tons_per_Year.csv', sep=',', index=False)
    
    total_county_WGT_2CSV = cnt_df.county_total_WGT.drop(columns = ['geometry','region_state','vkt_weight'])
    total_county_WGT_2CSV.to_csv(output_dir+'/'+'County_Total_Controled_Normalized_Emissions_Tons_per_Year.csv', sep=',', index=False)

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


if (Control_CALC == 'yes') or (Control_CALC == 'YES') or (Control_CALC == 'y') or (Control_CALC == 'Y'):
    aplly_control(Control_list, County_Emissions, county_SHP)
else:
    print('')





# =============================================================================
# Function to apply chemical into emissions of PM2.5, VOC and NOx
# =============================================================================
def chemical_speciation(Input_dir, County_Emissions_DataFrame, Activity_Data_DataFrame,
                        Chemical_Speciation_Table):
    print('******************************************************************')
    print('*****        Starting chemical speciation of emissions       *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')

    start_time = time.time()
    input_dir = Input_dir
    cnty_emis =  County_Emissions_DataFrame  #County_Emissions    #
    act_data  =  Activity_Data_DataFrame     #AD_SK               #
    ChemSpec  = Chemical_Speciation_Table    #Chemical_Spec_Table #
    nvhc = list(np.sort(act_data.fullname.vhc_name))
    pol_ChemSpec = ['PM2.5', 'VOC', 'NOx'] #list(ChemSpec.crossref.columns[3:])
    crossref = ChemSpec.crossref.copy()
    crossref['fullname'] = crossref.vehicle.str.cat(crossref[['types','fuel']], sep='_')
    County_Emissions_ChemSpec = cnty_emis.county_total_WGT.copy()
    grams_sec_pol = []
    
    for ipol in pol_ChemSpec:
        if ipol == 'PM2.5':
            print('')
            print('*****              Calculating for PM2.5             *****')
            print('')
            
            split_pol = list(ChemSpec.PM2_5_spec.pollutant.unique())
            ChemSpec.PM2_5_spec.loc[:,'mw_factor'] = ChemSpec.PM2_5_spec.loc[:,'fraction'] / \
                                                    ChemSpec.PM2_5_spec.loc[:,'mw']
            aux_chem  = ChemSpec.PM2_5_spec.pivot(index='profile', columns='pollutant',
                                                  values='mw_factor').fillna(0).reset_index(drop=False)
            aux_CR = crossref.loc[:,['fullname',ipol]]
            aux_CR = pd.merge(aux_CR,aux_chem, left_on = ipol, 
                                          right_on = ['profile'], how='left')
            grams_sec_pol.extend(ChemSpec.PM2_5_spec.pollutant.loc[ChemSpec.PM2_5_spec.mw == 1].drop_duplicates().to_list())
        elif ipol == 'NOx':
            print('')
            print('*****              Calculating for NOx               *****')
            print('')
            
            split_pol = list(ChemSpec.NOx_spec.pollutant.unique())
            ChemSpec.NOx_spec.loc[:,'mw_factor'] = ChemSpec.NOx_spec.loc[:,'fraction'] / \
                                                  ChemSpec.NOx_spec.loc[:,'mw']
            aux_chem  = ChemSpec.NOx_spec.pivot(index='profile', columns='pollutant',
                                                  values='mw_factor').fillna(0).reset_index(drop=False)
            aux_CR = crossref.loc[:,['fullname',ipol]]
            aux_CR = pd.merge(aux_CR,aux_chem, left_on = ipol, 
                                          right_on = ['profile'], how='left')
            grams_sec_pol.extend(ChemSpec.NOx_spec.pollutant.loc[ChemSpec.NOx_spec.mw == 1].drop_duplicates().to_list())

        elif ipol == 'VOC':
            print('')
            print('*****              Calculating for VOC               *****')
            print('')
            
            split_pol = list(ChemSpec.VOC_spec.pollutant.unique())
            ChemSpec.VOC_spec.loc[:,'mw_factor'] = ChemSpec.VOC_spec.loc[:,'fraction'] / \
                                                  ChemSpec.VOC_spec.loc[:,'mw']
            aux_chem  = ChemSpec.VOC_spec.pivot(index='profile', columns='pollutant',
                                                  values='mw_factor').fillna(0).reset_index(drop=False)
            aux_CR = crossref.loc[:,['fullname',ipol]]
            aux_CR = pd.merge(aux_CR,aux_chem, left_on = ipol, 
                                          right_on = ['profile'], how='left')
            grams_sec_pol.extend(ChemSpec.VOC_spec.pollutant.loc[ChemSpec.VOC_spec.mw == 1].drop_duplicates().to_list())

        else:
            print('')
            print('*** ERROR ABORT ***:')
            print('*****   There is no pollutant named as {0}     *****'.format(ipol))
            print('*****   Please review your input tables        *****')
            print('*****   Make sure that all tables have the same and equal pollutants names *****')
            print('')
            sys.exit('*** ERROR ABORT ***: Differeces in pollutants name. Check your input files')
        
        for splt_pol in split_pol:
            aux_ChemSpec         = aux_CR.loc[:,[splt_pol,'fullname']]
            aux_ChemSpec         = aux_ChemSpec.T
            aux_ChemSpec.columns = aux_ChemSpec.loc['fullname',:]
            aux_ChemSpec         =  aux_ChemSpec.drop(index='fullname')
            
            split_total = (cnty_emis.county_total_WGT.loc[(cnty_emis.county_total_WGT.pollutant == ipol)].reset_index(drop=True)).copy()
            split_total['pollutant'] = splt_pol
            apply_ChemSpec = (aux_ChemSpec.append([aux_ChemSpec] * (split_total.shape[0] - 1),ignore_index=True)).reset_index(drop=True)
            split_total.loc[:,nvhc] = (split_total.loc[:,nvhc]).mul(apply_ChemSpec, axis='columns')
            County_Emissions_ChemSpec = County_Emissions_ChemSpec.append(split_total,ignore_index=True)
    
    for jpol in ['CO','SOx']:
        if jpol in list(County_Emissions_ChemSpec.pollutant.unique()):
            print('')
            print('*****              Calculating for {0}                *****'.format(jpol))
            print('')
            
            if (jpol == 'CO') or (jpol == 'co'):
                MW = 28.01  #molar mass for CO = 28.01 g/mol
            elif (jpol == 'SOx') or (jpol == 'SO2'):
                print('*****              Coverting SOx to SO2              *****')
                MW = 64.066 #molar mass for SO2 =  64.066 g/mol
            
            drop_rows = list(County_Emissions_ChemSpec.loc[County_Emissions_ChemSpec.pollutant.isin([jpol]),nvhc].index)
            aux_jpol = County_Emissions_ChemSpec.loc[County_Emissions_ChemSpec.pollutant.isin([jpol])].reset_index(drop=True)
            aux_jpol.loc[:,nvhc] = aux_jpol.loc[:,nvhc].div(MW, axis='columns')  
            if jpol == 'SO2' : aux_jpol.loc[:,'pollutant'] = 'SO2'
            County_Emissions_ChemSpec = County_Emissions_ChemSpec.drop(index=drop_rows).reset_index(drop=True)
            County_Emissions_ChemSpec = County_Emissions_ChemSpec.append(aux_jpol, ignore_index=True, sort=False).reset_index(drop=True)
    
    if 'PM2.5' in list(County_Emissions_ChemSpec.pollutant.unique()): grams_sec_pol.extend(['PM2.5'])
    if 'PM10' in list(County_Emissions_ChemSpec.pollutant.unique()): grams_sec_pol.extend(['PM10'])
    moles_sec_pol = pd.DataFrame({ 'pollutant' : list(County_Emissions_ChemSpec.pollutant.unique())})
    moles_sec_pol = moles_sec_pol.loc[~moles_sec_pol.pollutant.isin(grams_sec_pol)].reset_index(drop=True)
    grams_sec_pol = pd.DataFrame({ 'pollutant' : grams_sec_pol})
    
    County_Emissions_ChemSpec['total_emis'] = County_Emissions_ChemSpec.loc[:,nvhc].sum(axis=1)
    County_Emissions_ChemSpec = County_Emissions_ChemSpec.fillna(0.0)
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

County_Emissions_ChemSpec = chemical_speciation(input_dir, County_Emissions,
                                                AD_SK, Chemical_Spec_Table)

# =============================================================================
# Function to generate hourly emissions and Gridding
# =============================================================================

nvhc = list(np.sort(AD_SK.fullname.vhc_name))

auxTP = TempPro.diurnalPro.groupby(['Day','fullname']).mean().reset_index(drop=False)

diurnalPro_out = pd.DataFrame(columns=['DateTime']+list(np.sort(auxTP.fullname.unique())))
for iday in run_period.Day.unique():
    auxDP = auxTP.loc[auxTP.Day == iday].T
    auxDP.columns = list(auxDP.loc['fullname'])
    auxDP = auxDP.drop(index=['Day','fullname']).reset_index(drop=True)
    auxDP = auxDP.sort_index(axis=1)
    auxDP.insert(loc=0, column='DateTime', value=pd.date_range(start=iday, freq='H', periods=24))
    
    diurnalPro_out = diurnalPro_out.append(auxDP,ignore_index=True, sort=False)


aux_hourly = County_Emissions_ChemSpec.chemspec_emissions.loc[:,['region_cd','pollutant']+list(AD_SK.fullname.vhc_name)]
hourly_emissions = County_Emissions_ChemSpec.chemspec_emissions.loc[:,['region_cd','pollutant']]
for idt in diurnalPro_out.DateTime:
    auxDT = diurnalPro_out.loc[diurnalPro_out.DateTime == idt, nvhc].reset_index(drop=True)
    auxEmis = aux_hourly.copy()
    for icol in auxDT.columns:
        auxEmis[icol] = aux_hourly[icol] * auxDT[icol].values
    hourly_emissions.insert(loc=hourly_emissions.shape[1], column=idt, value=auxEmis.loc[:,nvhc].sum(axis=1))

gridded_emissions = pd.DataFrame()
for ipol in hourly_emissions.pollutant.unique():

    grd_all = pd.merge(roads_RGS.surrogate.loc[:,['grid_id','region_cd','weight_factor']],
                   hourly_emissions.loc[hourly_emissions.pollutant == ipol], left_on='region_cd', right_on='region_cd', how='left').fillna(0)
    
    for icol in run_period.DateTime:
        grd_all.loc[:,icol] = grd_all.loc[:,icol] * grd_all.weight_factor.values
    
    grd_all.loc[:,'pollutant'] = ipol
    grd_all = grd_all.groupby(['grid_id','pollutant']).sum().drop(columns=['weight_factor','region_cd']).reset_index(drop=False)
    grd_all = pd.merge(grd_all,roads_RGS.grid.loc[:,['geometry','grid_id','col','row']], left_on='grid_id', right_on='grid_id', how='left')
    grd_all = gpd.GeoDataFrame(grd_all, crs=roads_RGS.grid.crs, geometry='geometry')
    
    grid_id_list = list(grd_all.grid_id.unique())
    gridded_emissions = gridded_emissions.append(grd_all, ignore_index=True, sort=False)

drop_zero = gridded_emissions.loc[gridded_emissions.pollutant == 0].index
gridded_emissions = gridded_emissions.drop(index=list(drop_zero)).reset_index(drop=True)
#releasing memory of gridding dataframe
grd_all = []

#getting the number of vars and the VAR-LIST
NVARS = len(gridded_emissions.pollutant.unique())
VAR_LIST = []
for kpol in gridded_emissions.pollutant.unique():
    VAR_LIST.append('{:<16}'.format(kpol))
s = pd.Series(VAR_LIST)
VAR_LIST = s.str.cat(sep='')

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

for ipol in gridded_emissions.pollutant.unique():

    pol_grd = gridded_emissions.loc[gridded_emissions.pollutant == ipol].reset_index(drop=True)
    pol_grd = pol_grd.drop(columns=['grid_id','pollutant', 'geometry'])

    pol_array = np.zeros((df_times.shape[0], 1, row,col))
    
    for idx ,idt in zip(range(0,df_times.shape[0]), df_times):
        aux_grd = pol_grd.loc[:,idt].to_numpy()
        aux_grd = aux_grd.reshape((row,col)) #np.flipud(aux_grd.reshape((row,col)))
        pol_array[idx,0,:,:] = aux_grd
    
    pol_array = pol_array.astype(np.float32)
    
    nc_var    = ncfile.createVariable(ipol , np.float32 , ('TSTEP','LAY', 'ROW','COL'))
    nc_var[:] = pol_array
    if ipol in County_Emissions_ChemSpec.grams_pol.pollutant.to_list():
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



def plot_figures_chart(County_Emissions, AD_SK, gridded_emissions, roads_RGS, 
                       plot_figures = 'no', plot_24 = 'no', pol_list = [], ):
    print('')
    print('******************************************************************')
    print('*****              The Plot option is YES.                   *****')
    print('***** The total emissions of CO, NOx, PM10, PM2.5, and VOC   ')
    print('******************************************************************')
    print('')  
    start_time_plot = time.time()
    
    County_Emissions = County_Emissions
    AD_SK = AD_SK
    gridded_emissions = gridded_emissions
    roads_RGS = roads_RGS
    
    if len(pol_list) == 0:
        polls = ['CO','NOx','PM10','PM2.5','VOC']
    else: 
        pol_list = pol_list
    
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
                               vmax = ccc.total_emis.max()*0.40)# ,edgecolor='0.8')
            
            ax.set_title('Total Emissions of {0}-{1} [ton/yr]'.format(case_name,ipol), fontsize = 22)
            sm = plt.cm.ScalarMappable(cmap=GrayJet, 
                                       norm=plt.Normalize(vmin=ccc.total_emis.min(),
                                                          vmax=ccc.total_emis.max() * 0.50))
            sm._A = []
            cbar = fig.colorbar(sm)
            #    ax.axis('off')
            name = '{0}_link_emis_SK_Total_emis_Fig_ax_{1}.png'.format(case_name,ipol)
            fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=300)
            plt.close()
        
        for ipol in polls:
            ddd = County_Emissions.county_total_WGT.loc[County_Emissions.county_total_WGT.pollutant == ipol]
            fig, ax = plt.subplots(1, figsize = (13,8))
            ddd.plot(column='total_emis', cmap=GrayJet, linewidth=0.8, ax=ax,
                               vmax = ddd.total_emis.max()*0.40)# ,edgecolor='0.8')
            
            ax.set_title('Total Emissions of {0}-{1} [ton/yr]'.format(case_name,ipol), fontsize = 22)
            sm = plt.cm.ScalarMappable(cmap=GrayJet, 
                                       norm=plt.Normalize(vmin=ddd.total_emis.min(),
                                                          vmax=ddd.total_emis.max() * 0.40))
            sm._A = []
            cbar = fig.colorbar(sm)
            #    ax.axis('off')
            name = '{0}_county_emis_SK_Total_emis_{1}.png'.format(case_name,ipol)
            fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=300)
            plt.close()
    
        #Piechart Plots
        fuels = list(AD_SK.fuels.fuels)
        for ipol in polls:
            pol_df = County_Emissions.county_total_WGT.loc[County_Emissions.county_total_WGT.pollutant == ipol]
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
            fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=300)
            plt.close()
        
        #Gridded emissions Plots
        for ipol in polls:
            plot_pol = gridded_emissions[gridded_emissions.pollutant == ipol]
            plot_pol['total_emis'] = plot_pol.loc[:,run_period.DateTime.loc[0:24]].sum(axis=1)
            vmin = 0
            vmax = 0.4 * plot_pol['total_emis'].max()
    
            fig, ax = plt.subplots(1, figsize = (13,8))
            plot_pol.plot(column='total_emis', cmap=GrayJet, linewidth=0.8, ax=ax,
                               vmax = vmax)# ,edgecolor='0.8')
            
            ax.set_title('{0} Grid Emissions - {1}'.format(case_name,ipol), fontsize = 22)
            sm = plt.cm.ScalarMappable(cmap=GrayJet, 
                                       norm=plt.Normalize(vmin=vmin,
                                                          vmax=vmax))
            sm._A = []
            cbar = fig.colorbar(sm)
            particles = ['PAL','PCA','PCL','PEC','PFE','PH2O','PK',
                         'PMC','PMG','PMN','PMOTHR','PNA','PNCOM',
                         'PNH4','PNO3','POC','PSI','PSO4','PTI']
            if ipol in particles:
                cbar.ax.set_ylabel(r'$grams{\cdot}day^{-1}$', fontsize=14)
            else:
                cbar.ax.set_ylabel(r'$mols{\cdot}day^{-1}$', fontsize=14)
    
            name = '{0}_Grid_emis_SK_Total_emis_{1}.png'.format(case_name,ipol)
            fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=200)
            plt.close()
    
    #Hourly gridded emissions Plots
    if plot_24 == 'yes':
    	for ipol in polls:
    		plot_pol = gridded_emissions[gridded_emissions.pollutant == ipol]
    		vmin = 0
    		vmax = 0.9 * plot_pol.loc[:,run_period.DateTime.loc[0:24]].max().max()
    		for itime in run_period.DateTime.loc[0:24]:
    			fig, ax = plt.subplots(1, figsize = (13,8))
    			plot_pol.plot(column=plot_pol[itime], cmap=GrayJet, linewidth=0.8, ax=ax,
    							   vmax = vmax)# ,edgecolor='0.8')
    			
    			ax.set_title('Grid Emissions of {0} - {1}'.format(case_name,ipol), fontsize = 22)
    			sm = plt.cm.ScalarMappable(cmap=GrayJet, 
    									   norm=plt.Normalize(vmin=vmin,
    														  vmax=vmax))
    			xpos = roads_RGS.surrogate.total_bounds[0] + 2000
    			ypos = roads_RGS.surrogate.total_bounds[3] - 2000
    			ax.text(xpos, ypos, '{0}:00H'.format(itime.hour), color='red', bbox=dict(facecolor='white', alpha=1))
    			sm._A = []
    			cbar = fig.colorbar(sm)
    			cbar.ax.set_ylabel(r'$mols{\cdot}s^{-1}$', fontsize=12)
    			#    ax.axis('off')
    			if itime.hour < 10:
    				name = 'Hourly_{0}_Grid_emis_SK_Total_emis_{1}_hour_0{2}.png'.format(case_name,ipol,itime.hour)
    			else:
    				name = 'Hourly_{0}_Grid_emis_SK_Total_emis_{1}_hour_{2}.png'.format(case_name,ipol,itime.hour)
    			fig.savefig(output_dir+'/plots/'+name, bbox_inches = 'tight',dpi=150)
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

plot_figures_chart(County_Emissions, AD_SK, gridded_emissions, roads_RGS, 
                       plot_figures = plot_figures, plot_24 = plot_24, pol_list = [])




