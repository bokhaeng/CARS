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
