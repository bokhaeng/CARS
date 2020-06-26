# Emission Factor
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
| 14| f           | real      | emission factor coefficient e  |
