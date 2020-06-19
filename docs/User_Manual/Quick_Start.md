# CARS Quick Start

## Set up the python environment for CARS:

## Install python 3

1. Download and install python 3 from [Anaconda](https://www.anaconda.com/products/individual)

### Install the third-party packages:
The third-party packages can be installed by GUI (Anaconda-Navigator), or manual install in the terminal with python 3 environment. The geopandas must be verion 0.6.1.

Here are the steps to install packages manually in terminal (not python console):

1. Install geopandas:
```
conda install geopandas=0.6.1
```
2. Install pandas:
```
conda install pandas
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

## Setup the parameter in CARS model:

1. Open and edit the file CARS_rev_07.py by text app or vim.
2. Setup the model variables and directory in CARS python script (line 22 to 37):
```
case_name = 'YOUR_CASE_NAME'
home_dir   = r'YOUR_CARS_HOME_DIRECTORY'
src_dir    = home_dir+'/src'
input_dir  = home_dir+'/input_country'
inter_dir  = home_dir+'/intermediate'#+case_name# for different input location
output_dir = home_dir+'/output_'+case_name
met_dir    = input_dir+'/metdata'
```
3. Setup the model process duration of output
```
STDATE = '2017-01-01'  # start date
ENDATE = '2017-01-31'  # end date
STTIME =  00            # start time 
RUNLEN =  744            # run length  # oneday is 24
```

If the test case is downloaded and the home_dir has been edited to right directory, the CARS python script can be processed for the test case. Other detail setup for each module will be explained in other chaper. 

## Run CARS model:
```
python CARS_rev_07.py
```
