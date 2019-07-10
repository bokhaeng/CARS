'''
INFORMATIONS HOW TO DOWNLOAD AND INSTALL PYTHON WITH ANACONDA TO RUN CARS
'''

# Go to ANACONDA website and download the installer of Python 3.7 for your system.

https://www.anaconda.com/distribution/#download-section

# After the download finished, install the Anaconda python package

The program will install the Anaconda prompt and Spyder.
Anaconda prompt is the a terminal that rely on MSdos commands.
This terminal is good to install Python Packages, launch Spyder and also run scripts.
Spyder is the GUI for editing and running python scripts. It really helps you on your coding, debuging and ploting. I recommend to use Spyder.


# Configuring your Anaconda to run CARS - Instaling packages

After everything is installed launch the Anaconda Prompt and type the follwing command to install the required libraries:

> conda install -c conda-forge -y pandas geopandas pyproj shapely datetime 

The option "-c conda-forge" was set because it is a updated repository with newer packages
and the "-y" is to say yes to install all the packages. If you do not set it, the prompt you ask you if 
you want to install the packages.

Press return to run the command and to install the packages.


# Running CARS test case

Download the CARS_test_case folder and untar it in your work diretory.
This work directory is important to run CARS, because it is the home directory in the script.

After extraction, thre will be in the folder:
data_pre_processor_CARS_draft_rev_03.py
input_seoul

CARS test case has all necessary input data under "input_seoul" folder.
You may want to take a look on it, but do not change the names or delete the files. This will make CARS stop to run.

You need to edit the "home_dir" variable inside the data_pre_processor_CARS_draft_rev_03.py
You can do it by opening it on text editor and changing it

home_dir   = r'C:\Users\user01\CARS_test_Case'

***** WARNING *****: Remember to set the path inside the '', as shown.

Or you can open Spyder, which has its own python editor, to edit it (Recommended option)

*****   DO NOT DOUBLE CLICK ON SPYDER ICON TO OPEN IT.  *****

First, you need to open Anaconda prompt.
Type the follwing command to open spyder:

spyder --hide-console 

This metohds is need because Spyder needs Enviroment Variables which are set when you open Anaconda Prompt.

After Spyder is ready (sometimes takes few seconds) you can open the data_pre_processor_CARS_draft_rev_03.py using the Open file option.

The data_pre_processor_CARS_draft_rev_03.py will be loaded in the screen and what you need to do is to change the home directory on the script:

home_dir   = r'C:\Users\user01\CARS_test_Case'

***** WARNING *****: Remember to set the path inside the '', as shown.

After this, just run the scripts.

Remember that CARS has two variables for plotting:

plot_figures = 'yes'   #'yes' or 'no'
plot_24 = 'no'        #'yes' or 'no'


plot_figures will plot the PIECHARTS, Total Emissions by County, Link and Grid

plot_24 will plot the hourly emissions. This option takes a bit longer to run because there are several images. If you set the Run length to 24,
There will 24 figures by pollutants, so be careful and pacient with this option.


If you do not want to Spyder, open data_pre_processor_CARS_draft_rev_03.py in some trext editor and change the home directory as mentioned.
Open the Anaconda prompt and navegate to the work directory.
type command to run:

ipython data_pre_processor_CARS_draft_rev_03.py


# After Run

After the run, there will be the intermediate folder and the output folder.
Check the output folder. You will see the plots, the shape file of the grid and county total .csv file.

Enjoy!
  



