# The shapefile format

## The requirement parameters in the road link shapefile.
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:--------------------------------------------|
| 1 | LINK_ID   | integer   | The ID of road link (e.g. 1680254500)      |            
| 2 | EMD_CD    | integer   | district code|
| 3 | EMD_ENG_NM| character | english district name|
| 4 | EMD_ENG_NM| character | alternative district name|
| 5 | ROAD_RANK | integer   | Road type (e.g. 101,102,...108)            |
| 6 | spd       | float     | road speed (km/hour)|
| 7 | SHAPE_STLe| integer   | Link length (meters)|
| 8 | VKT       | float     | VKT (km)|
| 9 | geometry  | float     | The shapefile geometry data|

## The requirement parameters in the county shapefile.
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:--------------------------------------|
| 1 | EMD_CD     | integer   | district code|            
| 2 | EMD_ENG_NM | integer   | english district name|
| 3 | EMD_NM     | integer   | original district name|
| 4 | geometry   | shapefile | The shapefile geometry data|
