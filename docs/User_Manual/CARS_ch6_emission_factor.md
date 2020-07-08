# Emission Factor table format for CARS input
The detail information of the EFs table with full descriptions is shown in emission factor table. These EFs are sorted by SCC, and varied by the age of vehicle as well as vehicle-road specific average speed.  Ambient temperature condition is only applicable for diesel fuel vehicle engines. Thus, depending on the type of vehicle and year of manufacture from Table below, appropriate EFs will be used to compute emissions inventory.

## Emission Factor table format
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:--------------------------------------------|
| 1 | vehicle     | character | vehicle type name (e.g. sedan)             |
| 2 | engine      | character | engine type name (e.g. compact)            |
| 3 | scc         | character | source category code 7 characters (e.g. 7010201)|
| 4 | fuel        | character | fuel name (e.g. Gasoline)              |
| 5 | pollutant   | character | pollutant (e.g. VOC, NOx, PM2.5)       |
| 6 | years       | integer(8)| vehicle manufacture year (e.g.20110000)|
| 7 | temperatures| character | temperature bins (e.g. LE0, GT0)       |
| 8 | speed       | real      | vehicle speed bins (e.g. LE65, GT65)   |
| 9 | a           | real      | emission factor coefficient a  |
| 10| b           | real      | emission factor coefficient b  |
| 11| c           | real      | emission factor coefficient c  |
| 12| d           | real      | emission factor coefficient d  |
| 13| e           | real      | emission factor coefficient e  |
| 14| f           | real      | emission factor coefficient f  |



# Vehicle Speed
This allows CARS to represent more realistic speed variation happening on each road type by calculating emissions factor by speed bin and then applying the weight factor by each speed bin. It allows CARS to generate a weighted emissions factor for each road and vehicle.

## Average Speed Distribution table format
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:--------------------------------------------|
| 1 | road type   | integer   | Road type (e.g. 101,102,...108)            |
| 2 | Speed bin 1 | float     | Average speed profile fraction from 0 to 1.|
| 3 | Speed bin 2 | float     | Average speed profile fraction from 0 to 1.|
| 4 | Speed bin 3 | float     | Average speed profile fraction from 0 to 1.|
| 5 | Speed bin 4 | float     | Average speed profile fraction from 0 to 1.|
| 6 | Speed bin 5 | float     | Average speed profile fraction from 0 to 1.|
| 7 | Speed bin 6 | float     | Average speed profile fraction from 0 to 1.|
| 8 | Speed bin 7 | float     | Average speed profile fraction from 0 to 1.|
| 9 | Speed bin 8 | float     | Average speed profile fraction from 0 to 1.|
| 10| Speed bin 9 | float     | Average speed profile fraction from 0 to 1.|
| 11| Speed bin 10| float     | Average speed profile fraction from 0 to 1.|
| 12| Speed bin 11| float     | Average speed profile fraction from 0 to 1.|
| 13| Speed bin 12| float     | Average speed profile fraction from 0 to 1.|
| 14| Speed bin 13| float     | Average speed profile fraction from 0 to 1.|
| 15| Speed bin 14| float     | Average speed profile fraction from 0 to 1.|
| 16| Speed bin 15| float     | Average speed profile fraction from 0 to 1.|
| 17| Speed bin 16| float     | Average speed profile fraction from 0 to 1.|
