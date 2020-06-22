# Vehicle Activity Data
Activity data is one of the mandatory required input data to run CARS, and the format is described in Table 1. Each line contains the following information to represent daily total vehicle kilometers traveled (VKT) from one vehicle. The total number of lines in this activity data input file represents the total number of vehicles registered. These VKT values are sorted by an 8-digit region code and manufactured date.


| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:---------------------------------------|
| 1 | fuel             | character | Fuel type Name (e.g. gasoline, diesel) |
| 2 | vehicle          | character | Vehicle type Name (e.g. sedan, truck)  |
| 3 | engien           | character | Engine type Name (e.g compact, full)   |
| 4 | daily_vkt        | real      | Daily Total Vehicle Kilometer Traveled |
| 5 | region_code      | integer(8)| Regional County Code (e.g. 11290127)   |
| 6 | manufacture date | integer(8)| Year of Emissions Factor (e.g.20110000)|
