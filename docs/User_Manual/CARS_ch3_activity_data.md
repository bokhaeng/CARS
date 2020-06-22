# Vehicle Activity Data
Activity data is one of the mandatory required input data to run CARS, and the format is described in Table 1. Each line contains the following information to represent daily total vehicle kilometers traveled (VKT) from one vehicle. The total number of lines in this activity data input file represents the total number of vehicles registered. These VKT values are sorted by an 8-digit region code and manufactured date.

## VKT data format
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:---------------------------------------|
| 1 | fuel             | character | fuel type name (e.g. gasoline, diesel) |
| 2 | vehicle          | character | vehicle type name (e.g. sedan, truck)  |
| 3 | engien           | character | engine type name (e.g compact, full)   |
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
