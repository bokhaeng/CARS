# The optional table format for district emission calculation in CARS

## the cold start vehicle tables format.
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:-------------------------------|
| 1 | vehicle| character   | Vehicle type |     
| 2 | engine | character   | Type name. Needs to match the Activity data names|
| 3 | fuel   | character   | Fuel name. Needs to match the Activity data names|
| 5 | ROAD_RANK | integer   | Road type (e.g. 101,102,...108)            |

## the road restriction tables format.
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:-------------------------------|
| 1 | vehicle| character   | Vehicle type. |     
| 2 | engine | character   | Type name. Needs to match the Activity data names|
| 3 | fuel   | character   | Fuel name. Needs to match the Activity data names|
| 5 | ROAD_RANK | integer   | The restriction of Road type (e.g. 101,102,...108)|
| 5 | ROAD_RANK | integer   | The restriction of Road typeRoad type (e.g. 101,102,...108)            |
| 5 | ROAD_RANK | integer   | The restriction of Road typeRoad type (e.g. 101,102,...108)            |
| 5 | ROAD_RANK | integer   | The restriction of Road typeRoad type (e.g. 101,102,...108)            |
| 5 | ROAD_RANK | integer   | The restriction of Road typeRoad type (e.g. 101,102,...108)            |
| 5 | ROAD_RANK | integer   | The restriction of Road typeRoad type (e.g. 101,102,...108)            |
| 5 | ROAD_RANK | integer   | The restriction of Road typeRoad type (e.g. 101,102,...108)            |
| 5 | ROAD_RANK | integer   | The restriction of Road typeRoad type (e.g. 101,102,...108)            |

## the deterioration tables, the table is varied by fuel types (gasoline, diesel, and LPG)
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:-------------------------------|
| 1 | vehicle| character   | Vehicle type. |     
| 2 | engine | character   | Type name. Needs to match the Activity data names|
| 3 | profile (SCC)| integer | Profile number to identify the vehicle, fuel and emission behavior |
| 4 | fuel   | character   | Fuel name. Needs to match the Activity data names|
| 5 | pollutant| character   | Pollutant. Each pollutant (VOC, NOx and PM2.5) should have its own file|
| 6 | Manufacture date | integer | year (e.g., 20170000)|
| 7 | 1   | real | hourly fraction|
| 8 | 2   | real | hourly fraction|
| 9 | 3   | real | hourly fraction|
| ... | ... | ... | ...|
| 20| 14  | real | hourly fraction|
| 21| 15  | real | hourly fraction|
| 22| 16  | real | hourly fraction|

## the control strategy list table
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:-------------------------------|
| 1 | vehicle| character   | Vehicle type. |     
| 2 | engine | character   | Type name. Needs to match the Activity data names|
| 3 | fuel   | character   | Fuel name. Needs to match the Activity data names|
| 4 | Manufacture date | integer | year (e.g., 2017)|
| 5 | region_cd     | integer   | district code or state code|   
| 5 | control_factor| integer   | 0-100, (%) the emission fractions under control |   
