# Installation and setup variables

## Install python 3

1. Download and install python 3 from [Anaconda](https://www.anaconda.com/products/individual)

## Install the third-party packages:
The third-party packages can be installed by GUI (Anaconda-Navigator), or manual install in the terminal with python 3 environment. The geopandas must be verion 0.6.1.

Here are the steps to install packages manually in terminal with python 3 environment (not python console):

1. Install geopandas:
```
conda install geopandas>=0.7.0
```
2. Install pandas:
```
conda install pandas>=1.1.2
```
3. Install numpy:
```
conda install numpy
```
4. Install matplotlib:
```
conda install matplotlib
```
5. Install netCDF:
```
conda install netCDF4
```
6. Install pyproj:
```
conda install pyproj>=6.2.1
```
## The list of variables CARS operating script:

1. Open and edit the file CARS_rev_07.py by text app or vim.
2. Setup the model variables and directory in CARS python script (line 22 to 233):

-	case_name: The case name is for distinguishing one case from another, and it is also used in the figures and logfile names.
-	home_dir: This is the full path for the home directory. It is where CARS has been executed (e.g. /path/user/CARS_testCase). The home directory is the head of all CARS folders, so all folders other should be inside it, like src_dir, input_dir, etc.
-	src_dir: This is the source directory and it is where CARS source code is located.
-	input_dir:  The input directory holds all the input folders, which are: activity_data; chemical_speciation; emissions_factor; metdata; shapes; and temporal_profile. These folders are default in CARS and the code will look for them. Their names were set to facilitate their identification and which files should be inside them.
-	inter_dir: Intermediate folders is to save the intermediate files, like dataframes and binary data processed by CARS. These files are constantly accessed by CARS, and their objective is to speed up the calculations and CARS processing.
- output_dir: The output directory is for holding the output files, reports, logs, and plots.
-	met_dir: The meteorology folder is the location of the grid description file (e.g. GRIDDESC), which is used for setting the CMAQ grid, define the surrogate and create the gridded emissions.
-	STDATE: The start date of the simulation (format YYYY-MM-DD).
-	ENDATE: The end date of the simulation (format YYYY-MM-DD).
-	STTIME: The start hour of the simulation. Usually, it is 00h.
-	RUNLEN: Run length is the number of hours the user needs to calculate. This variable and STTIME are mostly to apply the temporal variation in the gridded emissions.
-	temp_max; temp_mean; and temp_min: These three temperatures, maximum, mean and minimum are setted for the Cold Start and Evaporative emissions. These values can affect the total VOC emissions. Set it correctly.
-	activity_file: The activity data file name. The physical file must inside the home_directory/ input_dir/activity_data folder.
-	Emis_Factor_list: This is a list of emissions factor files. This is set as a list because often the emissions factors file are divided by fuel type. However the list can contain a single file with all emissions factor, the user just need to pass it as a between brackets (e.g. [‘my_emission_factor.csv’]. These physical files must be inside of the home_directory/ input_dir/emissions_factor  folder.
-	Cold_Start_list: Cold start list is just the vehicle which CARS should calculate cold start emissions (e.g gasoline vehicles). Thus, the users just need to insert the names of the vehicles in the cold start list file. It must be passed as a list to CARS (e.g. [‘cols_start_vehicles.csv’]).
-	avg_SPD_Dist_file: Average speed distribution file name. This file contains the eight-speed profile for the eight road types handled by CARS ( e.g. 101, 102, 103, 104, 105, 106, 107, 108). The file must be inside the home_directory/ input_dir/emissions_factor  folder.
-	link_shape: The link shapefile name. This is the shapefile that contains all roads inside the area of interest. The file must be inside the home_directory/ input_dir/shapes folder.
-	link_shape_att: Shapefiles have the attribute table which holds the data. CARS needs to accesses these data based on the attribute name, so, the user must pass a list with the attributes name which holds information like: link ID, County code, County name, County name in Korean, link type (e.g. 101, 102), speed of the link, link length and average AADDT. These information is really important for CARS calculation. Be sure these data are set correctly.
-	county_shape: Similar to link shape, the county shape must have all counties polygons for the area of interest and it also must be inside the home_directory/ input_dir/shapes folder.
-	county_shape_att: As the link attribute, these list-like names for the county shapefile must pass for CARS the attribute name which holds the information of: County code, County name and County name in Korean. This information is really important for CARS calculation. Be sure these data are set correctly.
-	temporal_profile_folder: As the name says, this is the temporal profile folder name and the physical folder must be inside the home_directory/ input_dir folder.
-	temporal_monthly_file; temporal_week_file; temporal_weekday_file; temporal_weekend_file: These variables are the monthly, week, weekday and weekend temporal profile file names. These files are to apply the temporal variation into the emissions and to generate the hourly gridded emissions for air quality modeling.  All these files must be in the temporal_profile_folder.
-	temporal_CrossRef: This is the temporal cross-reference file name and it is necessary to link the monthly, week, weekday and weekend temporal profiles with the vehicles and road type. This file must be in the temporal_profile_folder.
-	 chemical_profile_folder: This variable is to set the Chemical profile folder name. To folder holds the chemical profile file and the speciation cross-reference.
-	chemical_profile: The chemical profile file name. The data inside this file is crucial for the chemical speciation of all species processed by CARS. Make sure you are using the correct chemical profile for each vehicle type. This file must be in the chemical_profile_folder.
-	speciation_CrossRef: Similar to the temporal cross-reference, the speciation cross-reference file must have link between the vehicle with the chemical profile inside the chemical_profile, so CARS can get these information and split the pollutants into air quality modeling species. This file must be in the chemical_profile_folder.


| CARS variable | Example |
| :------------ |:-------------------:|
| case_name           | TEST1 |
| home_dir            | /path/for/CARS |
| src_dir             | /path/for/src |
| input_dir           | /path/for/input_files |
| inter_dir           | /path/for/intermediate_files |
| output_dir          | /path/for/output_dir|
| met_dir             | /path/for/neteorology |
| STDATE              | '2011-01-01'      |
| ENDATE              | '2011-01-02'      |
| STTIME              | 00      |
| RUNLEN              | 24(ond day), 744(one month)  |
| temp_max            | 17.8      |
| temp_mean           | 12.9      |
| temp_min            | 7.8      |
| activity_file       | 'activity_data.csv' |
| Emis_Factor_list    | ['cng.csv','diesel.csv','gasoline.csv'] |
| Cold_Start_list     | ['cold_start.csv'] |
| avg_SPD_Dist_file   | 'avgSpeedDistribution.csv' |
| link_shape          | 'road_shape_file.shp' |
| link_shape_att      | ['LINK_ID'  , 'EMD_CD' , 'EMD_ENG_NM', 'EMD_ENG_NM','ROAD_RANK', 'spd', 'SHAPE_STLe', 'VKT']|
| county_shape        | 'County.shp' |
| county_shape_att    | ['EMD_CD', 'EMD_ENG_NM', 'EMD_KOR_NM']  |
| temporal_profile_folder   | input_dir+'temproal_profile'      |
| temporal_monthly_file     | monthly_profile.csv      |
| temporal_week_file        | week_profile.csv      |
| temporal_weekday_file     | weekday_profile.csv      |
| temporal_weekend_file     | weekend_profile.csv      |
| temporal_CrossRef         | temporal_profile_CrossRef.csv      |
| chemical_profile_folder   | input_dir+'chemical_speciation'      |
| chemical_profile          | 'CB6.csv'      |
| speciation_CrossRef       | 'chem_profile_CrossRef.csv'      |

# Run CARS model
In the CARS model folder, make sure that src folder is at the same folder with run script (CARS_rev_08_script.py):
```
python CARS_rev_08_script.py
```
