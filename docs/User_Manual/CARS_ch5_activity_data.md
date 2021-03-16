# Vehicle Activity Data format for CARS input
Activity data is one of the mandatory required input data to run CARS, and the format is described in Table 1. Each line contains the following information to represent daily total vehicle kilometers traveled (VKT) from one vehicle. The total number of lines in this activity data input file represents the total number of vehicles registered. These VKT values are sorted by an 8-digit region code and manufactured date.

In the Process activity data module, the vehicle activity data is applied to generate vehicle VKT (VKT data format table). The Process emission factors module uses the [emission factor table]((https://github.com/CMASCenter/CARS/blob/master/docs/User_Manual/CARS_ch4_emission_factor.md)) to estimate the appropriate emission factor for every vehicle type with ambient temperature and road information data, AADT. The district Shapefile and the road shapefile are processed in the Process shape file module to generate the data arrays for the emission calculation. After the CARS input modules process the activity data, emission factors, and shapefiles, the Calculate district emissions module is applied to produce the emission rate (tons y-1) table for each vehicle type in each district or each link level.

## VKT data format
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:---------------------------------------|
| 1 | fuel             | character | fuel type name (e.g. gasoline, diesel) |
| 2 | vehicle          | character | vehicle type name (e.g. sedan, truck)  |
| 3 | engine           | character | engine type name (e.g compact, full)   |
| 4 | daily_vkt        | real      | daily total vehicle kilometer traveled |
| 5 | region_code      | integer(8)| regional county code (e.g. 11290127)   |
| 6 | manufacture date | integer(8)| vehicle manufacture year (e.g.20110000)|



## Fuel types supported in CARS
| Fuel type | Fuel description|
| :------------ |:---------------------------------------:|
| diesel     | Diesel |
| gasoline   | GAsoline |
| lpg        | Liquefied Petroleum Gas (LPG) |
| cng        | Compressed Natural Gas (CNG) |
| h-diesel   | Hybrid Diesel |
| h-gasoline | Hybrid Gasoline |
| h-lpg      | Hybrid LPG |
| h-cng      | Hybrid CNG |
| electric   | Electric |
