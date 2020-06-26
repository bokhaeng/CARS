# Introduction of CARS (Comprehensive Automobile Emissions Research Simulator)

Konkuk University and the University of North Carolina at Chapel Hill (UNC-CH) have been working closely to develop the bottom-up emission inventory for on-road and non-road mobile sources using their local emissions factors by vehicle, by fuel, by vehicle operating speed, and by ambient temperature (limited only to diesel hot engine exhaust emissions). While Konkuk University is responsible for collecting the latest of Korea’s official emissions factors and developing the local mobile source-related activity data, such as vehicle-specific total Vehicle Kilometer Travelled (VKT), vehicle ages, fuel type, and so on, UNC-CH is responsible for developing the new Korea on-road/non-road mobile emissions model called Comprehensive Automobile Research Modeling System (CARS) as a part of Korea Air Quality Modeling System (KAQMS) development to enhance the quality and accuracy of model prediction by providing a better quality emissions from mobile sources.

## CARS model scheme

![CARS scheme](https://github.com/CMASCenter/CARS/blob/master/docs/User_Manual/media/Picture1.png)

The CARS model has been developed on python 3 platform. Therefore, it applied the third-party python modules to process and prepare the output. The required python modules of CARS including the **geopandas**, **shapely.geometry**, and **csv** modules for reading the shapefile and table data. The **NumPy** and **pandas** modules are applied to operate the arrays and scientific calculations, the **pyproj** is for converting the coordinate systems, and **matplotlib** is for plotting figures. Further, the CARS model can also read and write the NetCDF file by the **NetCDF4** package.

The first process module in CARS is Loading_function_path; users can use this module to define and check the input files paths is correct. After setup the input files, CARS has six processes to calculate the emission rate and generate output. The following figure is the model design flow chart for those six modules and their functions. The six process modules include **Process activity data, Process emission factors, Process shape file, Calculate district emissions, Grid4AQM, and Plot figures**. The rectangles presented the data array, and rectangles with round edges are the functions in the modules of CARS. 



## CARS model process flow chart

![CARS flow chart](https://github.com/CMASCenter/CARS/blob/master/docs/User_Manual/media/Picture2.png)
