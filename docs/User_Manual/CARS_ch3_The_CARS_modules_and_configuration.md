The CARS model has been developed on python 3 platform. Therefore, it applied the third-party python modules to process and prepare the output. The required python modules of CARS including the **geopandas**, **shapely.geometry**, and **csv** modules for reading the shapefile and table data. The **NumPy** and **pandas** modules are applied to operate the arrays and scientific calculations, the **pyproj** is for converting the coordinate systems, and **matplotlib** is for plotting figures. Further, the CARS model can also read and write the NetCDF file by the **NetCDF4** package.

The first process module in CARS is Loading_function_path; users can use this module to define and check the input files paths is correct. After setup the input files, CARS has six processes to calculate the emission rate and generate output. The following figure is the model design flow chart for those six modules and their functions. The six process modules include **Process activity data, Process emission factors, Process shape file, Calculate district emissions, Grid4AQM, and Plot figures**. The rectangles presented the data array, and rectangles with round edges are the functions in the modules of CARS.

Each module is independent process to get the input data and generate the output data, because the CARS model design philosophy is to calculate the result fast in each module. Therefore the CARS model can modify the input for each module process and configure them without run unnecessary modules again.

## CARS model process flow chart

![CARS flow chart](https://github.com/CMASCenter/CARS/blob/master/docs/User_Manual/media/Picture2.png)


## The CARS module configuration:

### The **Process_activity_data** module
| CARS variable | Variable description| Example (default)|
| :------------ |:-------------------|:-------|
| process_activity_data | Apply the **process_activity_data** | 'yes'|
| activity_file | The location and the name of activity data | 'cars_input_v1.csv'|


### The **Process_emissions_factor** module
| CARS variable | Variable description| Example (default)|
| :------------ |:-------------------|:-------|
| Process_emissions_factor | Apply the **process_emission_factor** ['yes' or 'no']| 'yes'|
| Emis_Factor_list | The filename lists of the emission factor tables| ['gasoline.csv','diesel.csv','cng.csv','lpg.csv']
| temp_max        | The average of monthly maximum temperatures  | 17.8|
| temp_mean       | The annual average temperatures              | 12.9|
| temp_min        | The average of monthly minimum temperatures  | 7.8 |

### The **Process_shape_file** module
| CARS variable | Variable description| Example (default)|
| :---------------|:-------------------|:-------|
| Process_shape   | Apply the **Process_shape_file** ['yes' or 'no']| 'yes'|
| link_shape      | The location and name of road shapefile | '/shapes'+'/Road_by_County_SK_AADT_UTM52N.shp'|
| link_shape_att  | The attribute names in road shapefile | ['LINK_ID'  , 'EMD_CD' , 'EMD_ENG_NM', 'EMD_ENG_NM','ROAD_RANK', 'spd', 'SHAPE_STLe', 'VKT']|
| county_shape    | The location and name of county shapefile | '/shapes'+'/COUNTY.shp'|
| county_shape_att| The attribute names in county shapefile | ['EMD_CD', 'EMD_ENG_NM', 'EMD_KOR_NM']|


### The **Calculate_district_emission** module
| CARS variable | Variable description| Example (default)|
| :------------ |:-------------------|:-------|
| calculate_county_emissions | Apply the **Calculate_district_emission** ['yes' or 'no']| 'yes'|
| Cold_Start_list    | The filename of the cold start vehicle tables| ['cold_start_vehicles.csv']*
| avg_SPD_Dist_file  | The filename of the average speed distribution tables| ['avgSpeedDistribution.csv']
| process_road_restriction | Apply the road restriction function ['yes' or 'no']| 'yes'|
| road_restriction   | The filename of the road restriction tables| ['road_restriction_by_vehicle.csv']*
| Deterioration_factor     | Apply the deterioration function ['yes' or 'no']|'yes'|
| Deterioration_list | The filename of the deterioration tables|* ['degradation_rate_Diesel.csv','degradation_rate_Gasoline.csv','degradation_rate_LPG.csv']
| control_strategy   | Apply the control_strategy ['yes' or 'no']      |'yes'|
| control_list       | The filename of the control strategy list table | ['control_factors.csv']*


### The **Grid4AQM** module
| CARS variable | Variable description| Example (default)|
| :------------ |:-------------------|:-------|
| grid4AQM      | Apply the **Grid4AQM** ['yes' or 'no']| 'no' |
| IOAPI         | Apply the IOAPI output format         | 'no' |
| grid_size     | The grid cell size (meters)           | 1000 |
| gridfile_name | The filename of the grid cell file (NetCDF) format  | met_dir+'/GRIDDESC'|
| Radius_Earth  | The radius length of earth            | 6370000.0|
| Datum         | The projection system                 | WGS84|
| outFridShape  | output the GRID shapefile             | 'no' |
| future_case   | estimate emissions for the future     | 'no' |
| temporal_profile_folder | The folder of temporal profiles     | input_dir+'/temporal_profile'|
| temporal_monthly_file   | The monthly temporal profile table  | 'monthly_profile.csv'|
| temporal_week_file      | The week day temporal profile table | 'week_profile.csv'|
| temporal_weekday_file   | The weekday hourly temporal profile table | 'weekday_profile.csv'|
| temporal_weekend_file   | The weekend hourly temporal profile table | 'weekend_profile.csv'|
| temporal_CrossRef_file  | The cross reference table for temporal profile | 'temporal_crossref.csv'|
| chemical_profile_folder | The folder of chemical profiles | input_dir+'/chemical_speciation'|
| chemical_profile        | The filename of chemical profile| 'CMAQ_CB6_criteria.csv'|
| speciation_CrossRef     | The cross reference table for chemical profile | 'chemical_profile_crossref.csv'|

### The **Plot figure** module
| CARS variable | Variable description| Example (default)|
| :------------ |:-------------------|:-------|
| plot_figures  | Apply the **plot figure** ['yes' or 'no'] | 'no'|
| plot_24       | Plot 24 hours animated plot               | 'no'|
| adj_scale     | Adjust the scale of plot                  | 0.4 |
| cmap          | Adjust the color bar                      |'jet'|
