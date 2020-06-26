# The model output process in CARS

## 1. Chemical mechanism and species for the chemical transport model

### Chemical profile table
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:-------------------------------|
| 1 | profile  | integer | Profile number. The user can add as many as desired.|            
| 2 | pollutant| character   | Pollutant. Each pollutant (VOC, NOx and PM2.5) should have its own file|
| 3 | fraction | real | Fraction|
| 4 | MW       | real | Molecular weight|


### Chemical profile cross reference table
| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:--------------------------------------|
| 1 | vehicle| character   | Vehicle name. Needs to match the Activity data names|            
| 2 | engine | character   | Type name. Needs to match the Activity data names|
| 3 | fuel   | character   | Fuel name. Needs to match the Activity data names|
| 4 | PM2.5  | integer | Profile number specified in Chemical profile file|
| 5 | VOC    | integer | Profile number specified in Chemical profile file|
| 6 | NOX    | integer | Profile number specified in Chemical profile file|
| 7 | SOX    | integer | Profile number specified in Chemical profile file|
| 8 | NH3    | integer | Profile number specified in Chemical profile file|
| 9 | CO     | integer | Profile number specified in Chemical profile file|

## 2. Control factor table format table

| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:---------------------------------------|
| 1 | vehicle       | character | vehicle type name (e.g. sedan, truck)  |
| 2 | engien        | character | engine type name (e.g. compact         |
| 3 | fuel          | character | Fuel type (e.g. diesel, gasoline, LPG) |
| 4 | Year          | integer   | 19XX-20XX |
| 5 | data          | character | All|
| 6 | region_cd     | integer   | 8 digits Region/county code (e.g. 11290127)|
| 7 | Control factor| integer   | the percentage (%)|

## 3. Temporal profile

### Monthly Profile table

| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:--------------------------------------|
| 1 | profile| integer | profile index. The user can add as many as desired (e.g. 1-20)|            
| 2 | jan    |real| monthly fraction of January|
| 3 | feb    |real| monthly fraction of February|
| 4 | mar    |real| monthly fraction of March|
| 5 | apr    |real| monthly fraction of April|
| 6 | may    |real| monthly fraction of May|
| 7 | jun    |real| monthly fraction of June|
| 8 | jul    |real| monthly fraction of July|
| 9 | aug    |real| monthly fraction of August|
| 10| sep    |real| monthly fraction of September|
| 11| oct    |real| monthly fraction of October|
| 12| nov    |real| monthly fraction of November|
| 13| dec    |real| monthly fraction of December|


### Week Profile table

| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:--------------------------------------|
| 1 |profile| integer | profile index. The user can add as many as desired (e.g. 1-20) |            
| 2 | mon | real| daily fraction of Monday |
| 3 | tue | real| daily fraction of Tuesday|
| 4 | wed | real| daily fraction of Wednesday|
| 5 | thu | real| daily fraction of Thursday|
| 6 | fri | real| daily fraction of Friday|
| 7 | sat | real| daily fraction of Saturday|
| 8 | sun | real| daily fraction of Sunday|


### Weekday hourly Profile table

| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:--------------------------------------|
| 1 | profile| integer | profile index. The user can add as many as desired (e.g. 1-20)|            
| 2 | 0   | real | hourly fraction|
| 3 | 1   | real | hourly fraction|
| 4 | 2   | real | hourly fraction|
| 5 | 3   | real | hourly fraction|
| ... | ... | ... | ...|
| 23| 22  | real | hourly fraction|
| 24| 23  | real | hourly fraction|
| 25| 24  | real | hourly fraction|

### Weekend hourly Profile table

| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:--------------------------------------|
| 1 | profile| integer | profile index. The user can add as many as desired (e.g. 1-20)|            
| 2 | 0   | real | hourly fraction|
| 3 | 1   | real | hourly fraction|
| 4 | 2   | real | hourly fraction|
| 5 | 3   | real | hourly fraction|
| ... | ... | ... | ...|
| 23| 22  | real | hourly fraction|
| 24| 23  | real | hourly fraction|
| 25| 24  | real | hourly fraction|

### Temporal profile cross reference table

| Column | Name | Type | Description|
| :-------- |:------------------:| :-----------|:--------------------------------------------|
| 1 | vehicle     | character | vehicle type name (e.g. sedan)         |
| 2 | engine      | character | engine type name (e.g. compact)        |
| 3 | fuel        | character | fuel name (e.g. Gasoline)              |
| 4 | road type   | integer   | Road type (e.g. 101,102,...108)        |
| 5 | Monthly| integer | Temporal profile index from Temporal profile files|  
| 6 | Weekly | integer | Temporal profile index from Temporal profile files|  
| 7 | Weekday| integer | Temporal profile index from Temporal profile files|  
| 8 | Weekend| integer | Temporal profile index from Temporal profile files|  
