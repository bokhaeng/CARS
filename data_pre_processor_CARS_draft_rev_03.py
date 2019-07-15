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
import pyproj
import matplotlib.pyplot as plt
matplotlib.use('Agg')
total_start_time = time.time()

'''
*****      Case name and path to the folders - Please, set for all directories      *****
'''
case_name = 'Seoul_1x1km' 
home_dir   = r'C:\Users\pedruzzi\OneDrive - University of North Carolina at Chapel Hill\0000_EUA_IE\001_mobile_source_inventory\CARS_source_code\CARS_test_Case'
src_dir    = home_dir+'/src'
input_dir  = home_dir+'/input_seoul'
inter_dir  = home_dir+'/intermediate'
output_dir = home_dir+'/output_seoul'
met_dir    = input_dir+'/metdata'

'''
*****      Dates and Runlenght      *****
'''
STDATE = '2017-01-01'  # start date
ENDATE = '2017-01-01'  # end date
STTIME = 00            # start time (HH)
RUNLEN = 24            # run length (HHH)


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
ambient_temp           = 15.0
activity_file          =  'seoul_AD_v3.1.csv'

Emis_Factor_list       = ['cng_v3.csv','diesel_v3.csv','gasoline_v3.csv','lpg_v3.csv',
                          'H-CNG_v3.csv','H-Diesel_v3.csv','H-Gasoline_v3.csv' ,'H-LPG_v3.csv'] #['gasoline.csv','diesel.csv','cng.csv','lpg.csv']

avg_SPD_Dist_file      = 'avgSpeedDistribution_rev_01.csv'

link_shape             = '/shapes'+'/seoul_eup_links_Avg_VKT_UTM52N.shp'
link_shape_att         = ['link_id'  , 'EMD_CD' , 'EMD_ENG_NM', 'EMD_KOR_NM',
                           'road_type', 'speed', 'length_2', 'Avg_VKT']
                       
county_shape           = '/shapes/seoul_eup_UTM52N.shp'
county_shape_att       = ['EMD_CD', 'EMD_ENG_NM', 'EMD_KOR_NM']


'''
***** PLOTTING *****
Do you want to create the plots?
*** WARNING ***: The plotting takes a bit longer to finish
'''
plot_figures = 'no'   #'yes' or 'no'
plot_24 = 'no'        #'yes' or 'no'


'''
'*****      Gridding options      *****'
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
grid2CMAQ = 'no'     #'yes' or 'no'
gridfile_name = met_dir+'/GRIDDESC'
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
temporal_CrossRef      = 'temporal_profile_CrossRef_Seoul.csv'


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

if not os.path.exists(output_dir+'/plots'):
    print('*** intermediate directory does not exist: Creating a new one')
    os.makedirs(output_dir+'/plots/')
    
if not os.path.exists(output_dir+'/LOGS'):
        print('*** LOG directory does not exist: Creating a new one')
        os.makedirs(output_dir+'/LOGS')


df_times = pd.date_range(start=STDATE, freq='H', periods=RUNLEN )
run_period = pd.DataFrame({'DateTime' : df_times})
run_period['dayofweek'] = run_period.DateTime.dt.strftime('%a').str.lower()  #%a or %A
run_period['month'] = run_period.DateTime.dt.strftime('%b').str.lower()  #%a or %A
run_period[['Day','Hour']] = run_period['DateTime'].dt.strftime('%Y-%m-%d_%H:%M:%S').str.split('_',expand=True)
run_period['TFLAG']  = run_period.DateTime.dt.strftime('%j,%H0000')


class EmissionFactor_table:
    def __init__(self, dataframe, name ):
        self.dataframe = dataframe
        self.name      = name.split('.')[0]

class Activity_Data_table:
    def __init__(self, dataframe, vhc_name, years, fuels, vehicle):
        self.data     = dataframe
        self.fullname = vhc_name
        self.years    = years
        self.fuels    = fuels
        self.vhc     = vehicle

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
    def __init__(self, County_Emissions, County_Emissions_GeoRef, County_Emissions_GeoRef_WGT,County, Years):
        self.county_by_yr     = County_Emissions
        self.county_total     = County_Emissions_GeoRef
        self.county_total_WGT = County_Emissions_GeoRef_WGT
        self.counties         = County
        self.years            = Years

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
def read_avgSpeedDistribution(input_dir, avg_Speed_Distribution_file, sep = ';'):
    start_time = time.time()
    input_dir = input_dir
    spd_file = avg_Speed_Distribution_file
    name = '{0}{1}{2}'.format(input_dir,'/emissions_factor/',spd_file)
    if os.path.exists(name) == True:
        print ('')
        print ('Reading Average Speed Distribution table ...')
        print (name)

        ASD = pd.read_csv(name, sep = sep).fillna(np.nan)

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

avgSpeedDist = read_avgSpeedDistribution(input_dir, avg_SPD_Dist_file, sep = ';')


# =============================================================================
# Function to read the Activity Data
# =============================================================================
def read_activity_data_csv_SK(input_dir, ad_file, sep = ';'):
    start_time = time.time()
    ad_file = ad_file
    now = dt.datetime.now()
    current_yr = now.year
#    ad_file = activity_file
#    sep = ','
    name = '{0}{1}{2}'.format(input_dir,'/activity_data/',ad_file)
    outdf = pd.DataFrame(columns=['fuel','vehicle','types','daily_vkt','region_cd',
                                  'manufacture_date'])
    if os.path.exists(name) == True:
        print ('')
        print ('Reading Activity Data table ...')
        print (name)
        for activity_data in pd.read_csv(name, chunksize=1500000, sep = sep, usecols = [0,1,2,3,4,5]):
            print(activity_data.shape)
            activity_data = activity_data.fillna(0)
#        activity_data = (pd.read_csv(name, sep = sep, usecols = [0,1,2,3,4,5], encoding = 'utf-8')).fillna(0)
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
    #        activity_data = activity_data.rename(columns={'Region_code' : 'region_cd'})
            activity_data.loc[:,'manufacture_date']   = (activity_data.loc[:,'manufacture_date'] / 10000).astype(int)
            activity_data['fullname'] = activity_data.vehicle.str.cat(activity_data[['types','fuel']], sep='_')
            
            vhc_fuels = pd.DataFrame({'fuels'  : list(activity_data['fuel'].unique())})
            vhc_cars  = pd.DataFrame({'vhc'  : list(activity_data['vehicle'].unique())})
            if activity_data.loc[:,'manufacture_date'].max() > current_yr:
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
            
            grouped = activity_data.groupby(['manufacture_date','region_cd','fullname']).sum()
            grouped = grouped.unstack().fillna(0)
            grouped.columns = [x[1] for x in grouped.columns]
            grouped.reset_index(inplace=True)
            grouped.columns = grouped.columns.str.lower()
            grouped = grouped.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
            outdf = outdf.append(grouped, ignore_index=True, sort=False)
#            vhc_names = pd.DataFrame({'vhc_name'  : list(activity_data.FullName.unique())})
#            vhc_years = pd.DataFrame({'vhc_years' : list(grouped.manufacture_date.unique())})
    #        out_table = Activity_Data_table(grouped, vhc_names,vhc_years)
    else:
        print('*** ERROR ABORT ***:  Emissions Factor file "" ', ad_file, ' "" does not exist!')
        sys.exit('CARS preProcessor can not read Emissions Factor file')

    outdf = outdf.groupby(['manufacture_date','region_cd']).sum()
    outdf = outdf.reset_index(drop=False)
    outdf = outdf.sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
    nvhc = list(np.sort(outdf.columns))
    nvhc.remove('manufacture_date')
    nvhc.remove('region_cd')
    vhc_names = pd.DataFrame({'vhc_name'  : nvhc})
    vhc_years = pd.DataFrame({'vhc_years' : list(outdf.manufacture_date.unique())})
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))
    return Activity_Data_table(outdf, vhc_names, vhc_years, vhc_fuels, vhc_cars)
## =============================================================================

AD_SK = read_activity_data_csv_SK(input_dir, activity_file, sep = ',')


# =============================================================================
# Function to read link level shapefile
# =============================================================================
def roads_grid_surrogate_inf(input_dir, File_Name, Link_Shape_att, GridDesc_file='', Radius_Earth='', Unit_meters = True):
    start_time = time.time()
    Link_ID_attr       = Link_Shape_att[0]   #= Link_ID_attr
    Region_CD          = Link_Shape_att[1]   #= Region_Code
    Region_NM          = Link_Shape_att[2]   #= Region_Name
    RD_name_attr       = Link_Shape_att[3]   #= RD_name_attr 
    RD_type_attr       = Link_Shape_att[4]   #= RD_type_attr
    Speed_attr         = Link_Shape_att[5]   #= Speed_attr
    Link_length        = Link_Shape_att[6]   #= Link_length
    VKT_attr           = Link_Shape_att[7]   #= VKT_attr
    file_name          = File_Name
    gridfile_name      = GridDesc_file
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

#    Link_ID_attr       = link_shape_att[0]
#    Region_CD          = link_shape_att[1]
#    Region_NM          = link_shape_att[2]
#    RD_name_attr       = link_shape_att[3]
#    RD_type_attr       = link_shape_att[4]
#    Speed_attr         = link_shape_att[5]
#    Link_length        = link_shape_att[6]
#    VKT_attr           = link_shape_att[7]
#    file_name          = link_shape
#    gridfile_name      = gridfile_name
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
        out_roads['geometry']      = out_roads.buffer(0.1) #changing link to area to do the overlay
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
                print('***** There are {0} columns and {1} rows'.format(cols,rows))

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
                for j in range(0,rows-1):
                    yini = lat_list[j]
                    yfin = lat_list[j+1] 
                    for i in range(0,cols-1):
                        grid_row.append(j)
                        grid_col.append(i)
                        xini = lon_list[i]
                        xfin = lon_list[i+1]
                        polygons.append(Polygon([(xini, yini), (xfin, yini),
                                                 (xfin, yfin), (xini, yfin), (xini, yini)])) 

                grid_ID = [x for x in range (1,len(polygons)+1)]
                grid = gpd.GeoDataFrame({'geometry': polygons, 'grid_id': grid_ID,
                                         'row'     : grid_row, 'col'    : grid_col})
                grid.crs = proj_grid.srs
                
                # exporting grid as shapefile
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
            grid.to_file(filename = output_dir+'/grid_{0}.shp'.format(case_name), driver='ESRI Shapefile',crs_wkt=prj)
        
        
        out_roads['weight_factor'] = out_roads['vkt_avg'] / out_roads['total_area'] #dividing the wgt factor by the link "area" for consistency with overlay  
        #creating the surrogate
        surrogate = gpd.overlay(out_roads, grid, how='intersection').reset_index(drop=True)
        surrogate['split_area'] = surrogate.area
        surrogate.loc[:,'weight_factor'] = ((surrogate.loc[:,'weight_factor'] * \
                                             surrogate.loc[:,'split_area']) / \
                                             surrogate.loc[:,'total_area']).astype(float)
        

        for igeocd in out_roads.region_cd.unique():
            srgt_split_vkt = surrogate.weight_factor.loc[surrogate.region_cd   == igeocd].values / \
                 (surrogate.weight_factor.loc[surrogate.region_cd  == igeocd]).sum()
            surrogate.loc[surrogate.region_cd == igeocd, ['weight_factor']] = srgt_split_vkt
            
            vkt_split_county = out_roads.vkt_avg.loc[out_roads.region_cd   == igeocd].values / \
                 (out_roads.vkt_avg.loc[out_roads.region_cd  == igeocd]).sum()
            out_roads.loc[out_roads.region_cd == igeocd, ['vkt_split_county']] = vkt_split_county

            lnk_split_county = out_roads.link_length.loc[out_roads.region_cd   == igeocd].values / \
                 (out_roads.link_length.loc[out_roads.region_cd  == igeocd]).sum()
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
def read_emissions_factor_SK(input_dir, EmisFactor_list, sep = ';'):
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

EF_All_fuel_function = read_emissions_factor_SK(input_dir, Emis_Factor_list, sep = ';')


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
            print('*** WARNING *** Check your Emissions Factor input at EF temp depended ')
        elif type(i) != int:
            EF_spd_limit.append(i)
    EF_spd_limit = np.unique(EF_spd_limit)

    dat_EF = pd.DataFrame(columns=['fullname','pollutant','v','years','temperatures','spd','emis_fact'])
    for ispd in nspd:
        aux_EF = EF_All_fuel.loc[:,['fullname','pollutant','v','years','temperatures']]
        ef = EF_All_fuel.k * ((EF_All_fuel.a * (ispd ** EF_All_fuel.b)) + \
                              (EF_All_fuel.c * (ispd ** EF_All_fuel.d)) + \
                              (EF_All_fuel.f))

        aux_EF['spd']       = (aux_EF.years * 0) + ispd
        aux_EF['emis_fact'] = ef
        
        dat_EF = dat_EF.append(aux_EF, ignore_index=True, sort=False)
    
    dat_EF['ambient_temp'] = ((dat_EF.spd * 0 ) + ambient_temp).astype(int)
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
            if (tlim[:2] == 'GT') and (ambient_temp < int(tlim[2:])):
                idx = list(Temp_EF_Diesel[(Temp_EF_Diesel.temperatures  == tlim)].index)
                filter_temp_diesel.extend(idx)
            elif (tlim[:2] == 'LE') and (ambient_temp > int(tlim[2:])):
                idx = list(Temp_EF_Diesel[(Temp_EF_Diesel.temperatures  == tlim)].index)
                filter_temp_diesel.extend(idx)
        elif len(tlim) > 4:
            tlim2 = tlim.split('_')
            if (ambient_temp < int(tlim2[0][2:])) or (ambient_temp > int(tlim2[1][2:])):
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
    log_name = 'LOG_Missing_Vehicle_EmissionFactor.txt'
    with open(OutD+'/LOGS/'+log_name, 'w') as log:
        log.write('# ***** Vehicles without Emissions Factors *****\n')
        log.write('# \n')
        log.write('vehicle_types_fuel\n')
        for ivhc in ADvhc:
            if ivhc not in EFvhc:
                log.write('{0}\n'.format(ivhc))
    return       
check_VHC_AD_EF(EmisFactor_yr_spd, AD_SK, output_dir)


# =============================================================================
# Function to calculate the county emissions 
# =============================================================================
def calc_County_Emissions(Input_dir, Emissions_Factor_DataFrame, 
                          Activity_Data_DataFrame, Roads_DataFrame, County_DataFrame,Out_dir ):
    print('******************************************************************')
    print('*****           Starting calculation of emissions            *****')
    print('*****                    Please wait ...                     *****')
    print('******************************************************************')
    print('')
    start_time = time.time()
    input_dir  = Input_dir                   #Input_dir
    output_dir = Out_dir
    EF_yr_spd  = Emissions_Factor_DataFrame #EmisFactor_yr_spd #
    AD_yr_vhc  = Activity_Data_DataFrame    #AD_SK             #
    roads_DF   = Roads_DataFrame            #roads_RGS         #
    county_df  = County_DataFrame           #county_SHP         #
    
  
    temp_EF_df = EF_yr_spd.data.copy()
    nspd = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 89, 97, 105, 113, 121]
    EFyears = list(np.sort(temp_EF_df.years.unique()))
    ADyears = list(np.sort(AD_yr_vhc.data.manufacture_date.unique()))
    EF_yr_min = np.min(EFyears)
    AD_yr_max = np.max(ADyears)
    npol = list(EF_yr_spd.EF_polls.EF_poll)
    nyears = list(np.arange(EF_yr_min,AD_yr_max+1))
    ADvhc = list(np.sort(AD_yr_vhc.fullname.vhc_name))
    EFvhc = list(np.sort(EF_yr_spd.EF_fullname.EF_fullname))
    missVHC = []
    for j in ADvhc:
        if j not in EFvhc:
            missVHC.append(j)
    missVHC = list(dict.fromkeys(missVHC))
    
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
    inter_AD = inter_AD.groupby(['region_cd']).sum()
    inter_AD = inter_AD.reset_index(drop=False)
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
    AD_byRoad['road_type'].iloc[ISzeroRD] = 101
    AD_byRoad['road_type'] = AD_byRoad['road_type'].astype(int)
    
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
        AD_by_EF_df = pd.merge(AD_byRoad.loc[:,['manufacture_date','region_cd','road_type']], 
                       EF_avgSpd.loc[EF_avgSpd['pollutant'] == ipol,['years','pollutant','road_type']+ADvhc],
                       left_on = ['manufacture_date','road_type'],
                       right_on = ['years','road_type'], how='left').reset_index(drop=True)
                       

        AD_by_EF_df.loc[:,ADvhc] = AD_by_EF_df.loc[:,ADvhc].mul(AD_byRoad.loc[:,ADvhc])
        AD_by_EF_df = AD_by_EF_df.groupby(['manufacture_date','region_cd','pollutant']).sum().reset_index(drop=False)
        AD_by_EF_df = AD_by_EF_df.drop(columns=['road_type','years']).sort_values(by=['region_cd','manufacture_date']).reset_index(drop=True)
        
        Emissions_by_yr_county = Emissions_by_yr_county.append(AD_by_EF_df, ignore_index=True, sort=False)
    

    total_county = Emissions_by_yr_county.groupby(['region_cd','pollutant']).sum()
    total_county = total_county.reset_index(drop=False)

    Emissions_by_link = pd.merge(roads_DF.roads_df,total_county, left_on = ['region_cd'], 
                        right_on = ['region_cd'], how='left')
    
    Emissions_by_link.loc[:,ADvhc] = Emissions_by_link.loc[:,ADvhc].apply(lambda x: np.asarray(x) * \
                         Emissions_by_link.vkt_split_county.values)    

    #Calculation the total emissions for county and link level
    Emissions_by_yr_county['total_emis'] = Emissions_by_yr_county.loc[:,ADvhc].sum(axis=1)
    total_county['total_emis'] = total_county.loc[:,ADvhc].sum(axis=1)
    county_table = pd.DataFrame({'region_cd' : Emissions_by_yr_county.region_cd.unique()})
    years_table  = pd.DataFrame({'years' : Emissions_by_yr_county.manufacture_date.unique()})

#    cnt_info = county_df.loc[:,['geometry','region_cd']]
    total_county = pd.merge(total_county,county_df, left_on = ['region_cd'], 
                        right_on = ['region_cd'], how='left')
    nanGeometry = list(total_county.loc[total_county['geometry'].isna()].index)
    total_county = total_county.drop(index=nanGeometry).reset_index(drop=True)
    
    total_county = gpd.GeoDataFrame(total_county, crs=county_df.crs, geometry='geometry')

    
    aux_df_wgt = roads_DF.roads_df.groupby(['region_cd']).sum()
    aux_df_wgt = aux_df_wgt.reset_index(drop=False)
    aux_df_wgt = pd.merge(county_table, aux_df_wgt.loc[:,['region_cd','vkt_avg']], 
                          left_on=['region_cd'], right_on=['region_cd'], how='left')
    
    emislog_name = 'LOG_Calculation_emissions_by_County.txt'
    with open(output_dir+'/LOGS/'+emislog_name, 'w') as EF_log:
        for icd in aux_df_wgt.region_cd.loc[aux_df_wgt.vkt_avg.isna()]:
            EF_log.write('*****   WARNING   *****')
            EF_log.write('There is no VKT data for the county {0}! Please, review \
                         your input link shapefile'.format(icd))
    aux_df_wgt = aux_df_wgt.fillna(0)
    aux_df_wgt['vkt_weight'] = aux_df_wgt.loc[:,'vkt_avg'] /  aux_df_wgt.loc[:,'vkt_avg'].sum()
    aux_df_wgt = aux_df_wgt.sort_values(by=['region_cd']).reset_index(drop=True)
    
    total_county_WGT = pd.DataFrame(columns=list(total_county.columns))

    for ipol in total_county.pollutant.unique():
        aux_df = total_county.loc[total_county.pollutant == ipol]
        for ivhc in ADvhc:
            aux_df.loc[:,ivhc] = aux_df.loc[:,ivhc].sum()
        aux_df = pd.merge(aux_df, aux_df_wgt, 
                          left_on=['region_cd'], right_on=['region_cd'], how='left')
        aux_df.loc[:,ADvhc] = aux_df.loc[:,ADvhc].apply(lambda x: np.asarray(x) * aux_df.vkt_weight.values)
        
        total_county_WGT = total_county_WGT.append(aux_df, ignore_index=True, sort=False)
    
    total_county_WGT = gpd.GeoDataFrame(total_county_WGT, crs=county_df.crs, geometry='geometry')

    total_county.loc[:,ADvhc] = total_county.loc[:,ADvhc] * (365/1000000) # tons/year
    total_county_WGT.loc[:,ADvhc] = total_county_WGT.loc[:,ADvhc] * (365/1000000) # tons/year
    
    total_county['total_emis'] = total_county.loc[:,ADvhc].sum(axis=1)
    total_county_WGT['total_emis'] = total_county_WGT.loc[:,ADvhc].sum(axis=1)
    
    total_county_WGT = total_county_WGT.drop(columns = ['geometry','manufacture_date','region_name_SK','vkt_avg','vkt_weight'])
    total_county_WGT.to_csv(output_dir+'/'+'County_Total_emissions.csv', sep=',', index=False)
    
    run_time = ((time.time() - start_time))
    print('---     Elapsed time in seconds = {0}     ---'.format(run_time))
    print('---     Elapsed time in minutes = {0}     ---'.format(run_time/60))
    print('---     Elapsed time in   hours = {0}     ---'.format(run_time/3600))
    print('')
    print('******************************************************************')
    print('*****             Emissions calculation is done              *****')
    print('******************************************************************')
    print('')
       
    return Emissions_table(Emissions_by_yr_county, total_county, total_county_WGT, county_table, years_table)

County_Emissions = calc_County_Emissions(input_dir, EmisFactor_yr_spd,
                                       AD_SK, roads_RGS, county_SHP, output_dir)

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
#    input_dir = Input_dir
    cnty_emis = County_Emissions_DataFrame  #County_Emissions   # 
    act_data  = Activity_Data_DataFrame     # AD_SK               #
    ChemSpec  = Chemical_Speciation_Table   #Chemical_Spec_Table #
    nvhc = list(np.sort(act_data.fullname.vhc_name))
    pol_ChemSpec = ['PM2.5', 'VOC', 'NOx'] #list(ChemSpec.crossref.columns[3:])
    crossref = ChemSpec.crossref.copy()
    crossref['fullname'] = crossref.vehicle.str.cat(crossref[['types','fuel']], sep='_')
    County_Emissions_ChemSpec = cnty_emis.county_total_WGT.copy()

    for ipol in pol_ChemSpec:
        if ipol == 'PM2.5':
            print('')
            print('*****              Calculating for PM2.5             *****')
            print('')
            
            split_pol = list(ChemSpec.PM2_5_spec.pollutant.unique())
            ChemSpec.PM2_5_spec.loc[:,'fraction'] = ChemSpec.PM2_5_spec.loc[:,'fraction'] / \
                                                    ChemSpec.PM2_5_spec.loc[:,'mw']
            aux_chem  = ChemSpec.PM2_5_spec.pivot(index='profile', columns='pollutant',
                                                  values='fraction').fillna(0).reset_index(drop=False)
            aux_CR = crossref.loc[:,['fullname',ipol]]
            aux_CR = pd.merge(aux_CR,aux_chem, left_on = ipol, 
                                          right_on = ['profile'], how='left')
        elif ipol == 'NOx':
            print('')
            print('*****              Calculating for NOx               *****')
            print('')
            
            split_pol = list(ChemSpec.NOx_spec.pollutant.unique())
            ChemSpec.NOx_spec.loc[:,'fraction'] = ChemSpec.NOx_spec.loc[:,'fraction'] / \
                                                  ChemSpec.NOx_spec.loc[:,'mw']
            aux_chem  = ChemSpec.NOx_spec.pivot(index='profile', columns='pollutant',
                                                  values='fraction').fillna(0).reset_index(drop=False)
            aux_CR = crossref.loc[:,['fullname',ipol]]
            aux_CR = pd.merge(aux_CR,aux_chem, left_on = ipol, 
                                          right_on = ['profile'], how='left')
        elif ipol == 'VOC':
            print('')
            print('*****              Calculating for VOC               *****')
            print('')
            
            split_pol = list(ChemSpec.VOC_spec.pollutant.unique())
            ChemSpec.VOC_spec.loc[:,'fraction'] = ChemSpec.VOC_spec.loc[:,'fraction'] / \
                                                  ChemSpec.VOC_spec.loc[:,'mw']
            aux_chem  = ChemSpec.VOC_spec.pivot(index='profile', columns='pollutant',
                                                  values='fraction').fillna(0).reset_index(drop=False)
            aux_CR = crossref.loc[:,['fullname',ipol]]
            aux_CR = pd.merge(aux_CR,aux_chem, left_on = ipol, 
                                          right_on = ['profile'], how='left')    
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
    
    County_Emissions_ChemSpec['total_emis'] = County_Emissions_ChemSpec.loc[:,nvhc].sum(axis=1)
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
    return County_Emissions_ChemSpec


# =============================================================================
# Function to generate houly emissions and Gridding
# =============================================================================
County_Emissions_ChemSpec = chemical_speciation(input_dir, County_Emissions,
                                                AD_SK, Chemical_Spec_Table)
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


aux_hourly = County_Emissions_ChemSpec.loc[:,['region_cd','pollutant']+list(AD_SK.fullname.vhc_name)]
hourly_emissions = County_Emissions_ChemSpec.loc[:,['region_cd','pollutant']]
for idt in diurnalPro_out.DateTime:
    auxDT = diurnalPro_out.loc[diurnalPro_out.DateTime == idt, nvhc].reset_index(drop=True)
    auxEmis = aux_hourly.copy()
    for icol in auxDT.columns:
        auxEmis[icol] = aux_hourly[icol] * auxDT[icol].values
    auxEmis[idt] = auxEmis.loc[:,nvhc].sum(axis=1)
    hourly_emissions.insert(loc=hourly_emissions.shape[1], column=idt, value=auxEmis.loc[:,nvhc].sum(axis=1))

grd_all = pd.merge(roads_RGS.surrogate.loc[:,['grid_id','region_cd','weight_factor']],
               hourly_emissions, left_on='region_cd', right_on='region_cd', how='left').fillna(0)

for icol in run_period.DateTime:
#    grd_all.loc[:,list(diurnalPro_out.DateTime.unique())] = \
#    grd_all.loc[:,list(diurnalPro_out.DateTime.unique())].apply(lambda x: np.asarray(x) * grd_all.weight_factor.values)
    grd_all.loc[:,icol] = grd_all.loc[:,icol] * grd_all.weight_factor.values
    
grd_all = grd_all.groupby(['grid_id','pollutant']).sum().drop(columns=['weight_factor','region_cd']).reset_index(drop=False)
grd_all = pd.merge(grd_all,roads_RGS.grid.loc[:,['geometry','grid_id','col','row']], left_on='grid_id', right_on='grid_id', how='left')
grd_all = gpd.GeoDataFrame(grd_all, crs=roads_RGS.grid.crs, geometry='geometry')

#
#grd_plot = pd.merge(County_Emissions_ChemSpec,
#                    roads_RGS.surrogate.loc[:,['grid_id','region_cd','weight_factor']],
#                    left_on='region_cd', right_on='region_cd', how='left').fillna(0)
#
##grd_plot = pd.merge(grd_all,roads_RGS.grid.loc[:,['geometry','grid_id','col','row']], left_on='grid_id', right_on='grid_id', how='left')
#
#grd_plot = gpd.GeoDataFrame(grd_plot, crs=roads_RGS.grid.crs, geometry='geometry')
#
#grd_plot['total_emis'] = grd_plot['total_emis'].mul(grd_plot['weight_factor'])
#
#aaa = grd_plot.loc[grd_plot.pollutant =='CO']
#aaa.plot(column='total_emis', cmap=GrayJet)


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



start_time_plot = time.time()

if plot_figures == 'yes':
    
    polls = ['CO','NOx','PM10','PM2.5','VOC']
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
        
        ax.set_title('Total Emissions of Seoul-{0} [ton/yr]'.format(ipol), fontsize = 22)
        sm = plt.cm.ScalarMappable(cmap=GrayJet, 
                                   norm=plt.Normalize(vmin=ccc.total_emis.min(),
                                                      vmax=ccc.total_emis.max() * 0.50))
        sm._A = []
        cbar = fig.colorbar(sm)
        #    ax.axis('off')
        name = 'Seoul_link_emis_SK_Total_emis_Fig_ax_{0}.png'.format(ipol)
        fig.savefig(output_dir+'/'+name, bbox_inches = 'tight',dpi=300)
        plt.close()
    #
    #nvhc = list(np.sort(AD_SK.fullname.vhc_name))
    #igeocd = 11110101
    #itime = run_period.DateTime[0]
    
    polls = ['CO','NOx','PM10','PM2.5','VOC']
    for ipol in polls:
        ddd = ddd = County_Emissions.county_total_WGT.loc[County_Emissions.county_total_WGT.pollutant == ipol]
        fig, ax = plt.subplots(1, figsize = (13,8))
        ddd.plot(column='total_emis', cmap=GrayJet, linewidth=0.8, ax=ax,
                           vmax = ddd.total_emis.max()*0.40)# ,edgecolor='0.8')
        
        ax.set_title('Total Emissions of Seoul-{0} [ton/yr]'.format(ipol), fontsize = 22)
        sm = plt.cm.ScalarMappable(cmap=GrayJet, 
                                   norm=plt.Normalize(vmin=ddd.total_emis.min(),
                                                      vmax=ddd.total_emis.max() * 0.40))
        sm._A = []
        cbar = fig.colorbar(sm)
        #    ax.axis('off')
        name = 'Seoul_county_emis_SK_Total_emis_{0}.png'.format(ipol)
        fig.savefig(output_dir+'/'+name, bbox_inches = 'tight',dpi=300)
        plt.close()


  
    polls = ['CO','NOx','PM10','PM2.5','VOC']
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
        fig.savefig(output_dir+'/'+name, bbox_inches = 'tight',dpi=300)
        plt.close()

    for ipol in polls:
        plot_pol = grd_all[grd_all.pollutant == ipol]
        plot_pol['total_emis'] = plot_pol.loc[:,run_period.DateTime.loc[0:24]].sum(axis=1)
        vmin = 0
        vmax = 0.4 * plot_pol['total_emis'].max().max()

        fig, ax = plt.subplots(1, figsize = (13,8))
        plot_pol.plot(column='total_emis', cmap=GrayJet, linewidth=0.8, ax=ax,
                           vmax = vmax)# ,edgecolor='0.8')
        
        ax.set_title('Seoul Grid Emissions - {0}'.format(ipol), fontsize = 22)
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

        name = 'Seoul_Grid_emis_SK_Total_emis_{0}.png'.format(ipol)
        fig.savefig(output_dir+'/'+name, bbox_inches = 'tight',dpi=200)
        plt.close()


        
        

if plot_24 == 'yes':
	polls = ['CO','NOx','PM10','PM2.5','VOC']
	for ipol in polls:
		plot_pol = grd_all[grd_all.pollutant == ipol]
		vmin = 0
		vmax = 0.9 * plot_pol.loc[:,run_period.DateTime.loc[0:24]].max().max()
		for itime in run_period.DateTime.loc[0:24]:
			fig, ax = plt.subplots(1, figsize = (13,8))
			plot_pol.plot(column=plot_pol[itime], cmap=GrayJet, linewidth=0.8, ax=ax,
							   vmax = vmax)# ,edgecolor='0.8')
			
			ax.set_title('Grid Emissions of Seoul - {0}'.format(ipol), fontsize = 22)
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
				name = 'Hourly_Seoul_Grid_emis_SK_Total_emis_{0}_hour_0{1}.png'.format(ipol,itime.hour)
			else:
				name = 'Hourly_Seoul_Grid_emis_SK_Total_emis_{0}_hour_{1}.png'.format(ipol,itime.hour)
			fig.savefig(output_dir+'/'+name, bbox_inches = 'tight',dpi=150)
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

