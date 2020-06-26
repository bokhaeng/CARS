# Emission Rate Calculation

Mobile emission includes on-road and non-road sections. The non-road part is transport emissions but not happened on-network roads but off-network, such as aviation, railways, construction, and boats. While non-road emissions are important, we will focus on describing the complex onroad emissions calculation in this study. The following part explains the approach of the on-road emission processes where happens on on-network roads. The on-road emission Eon-road in CARS is defined in Eq. 2, including three emission processes (Ntziachristos and Samaras, 2000):

![CARS scheme](https://github.com/CMASCenter/CARS/blob/master/docs/User_Manual/media/Picture3.png)

- The hot exhaust emissions (E<sub>hot</sub>) represent the tailpipe emission from the vehicles when the internal combustion engine (ICE) combusts the fuel to generate energy under the normal operating temperature.

- The cold start emissions (E<sub>cold</sub>) represent the tailpipe emissions from the ICE when the vehicle is just started, and the ICE operational temperature is below than normal condition.

- The evaporative VOC emissions (E<sub>vap</sub>) represent the emissions from the fuel systems (fuel tanks, injection system, and fuel lines) of petrol vehicles.

## Hot exhaust emissions

The hot exhaust emission, which emitted from the tailpipe of the vehicle, is the waste gas from the combustion process in ICE. The processes of ICE combustion cycle are not ideal and cause incomplete combustion processes. Therefore, those incomplete fuel combustion processes emit hydrocarbons, carbon monoxide (COCO), and particulate matter (PM). The sulfur compound in the fuel are oxidized and become the sulfur oxides (SO<sub>x</sub>); the nitrogen oxides (NO<sub>x</sub>) are also produced during the combustion process because of the abundant nitrogen (N<sub>2</sub>) and oxygen (O<sub>2</sub>) in the atmosphere.

![Hot emissions](https://github.com/CMASCenter/CARS/blob/master/docs/User_Manual/media/Picture4.png)

The hot emission equation represents the calculation of district daily total hot exhaust emission rate (E<sub>hot; p,v,s</sub> : g d<sup>-1</sup>) of pollutant p. In Eq. 3, the daily VKTn,v,age (km d<sup>-1</sup>) of an individual vehicle is calculated in Eq. 1. EF<sub>hot; p,v,s</sub> (g km<sup>-1</sup>) is the hot exhaust emission factor of pollutant p when the vehicle type is v and vehicle average speed is s. The district total emission rate is the summation of all individual vehicle emission rate, which is VKTn,v multiplied by the emission factor (EF<sub>hot; p,v,s</sub>).

![Emission Factor](https://github.com/CMASCenter/CARS/blob/master/docs/User_Manual/media/Picture5.png)

The equation above this paragraph is applied to calculate the emission factors, EFhot;p,v,s (g km<sup>-1</sup>), which is a function of vehicle speed (s) with other empirical coefficients: a, b, c, d, f, k. The emission factor formula and those coefficients were developed by CAPSS (Lee et al., 2011a). Those coefficients are varied by vehicle type v, pollutants p, and vehicle speed s. The vehicle speed s affects the combustion efficiency of ICE, and these can change the emission rate and emission composition from the tailpipe.
## Cold emission

The cold start emission is an additional emission of the hot exhaust emission and happened when the vehicle is started under the cold status. The cold status of ICE is caused by the vehicle has been stopped for a long period, and ICE has become cool. This lower temperature of ICE is not an optimum condition to completely combust the fuel. Therefore, this status reduces combustion efficiency (CE) and increase hydrocarbon and CO emission (Jang et al., 2007). These types of emissions occur in all vehicle types. The CARS model, however, now only estimate the cold start emissions factor for the gasoline, diesel, and Liquefied Petroleum Gas (LPG) vehicles. Besides the vehicle type and engine type, the road type is another major factor to affect the cold start emission. It rarely occurs on highways, but more on parking lots next to these road types. Further, the CARS model assumed that vehicles with cold start emission behave like passenger cars.

The cold start emission, E<sub>cold</sub> (g day<sup>-1</sup>), which is similar to the hot exhaust emission equation, is calculated in the following equation. In this equation, the emission ratio of hot and cold exhaust emission (EF<sub>cold</sub>/EF<sub>hot</sub> -1) and the percentage of the travel distance in the cold engine.

![Cold emissions](https://github.com/CMASCenter/CARS/blob/master/docs/User_Manual/media/Picture6.png)

## Evaporative VOC emission

The fuel in the vehicle evaporated into the atmosphere is called the evaporative emission. This emission happened in the petrol fuel system only, such as tanks, injection systems, and fuel lines of vehicles. The evaporative emission occurred in diesel vehicles can be ignored because of the low vapor pressure of diesel. The primary sources of evaporative emission from a vehicle are breathing losses through tank vent and fuel permeation/leakage. The CARS model considered the emission inventory guidebook (EEA, 2019) and made three mechanisms to estimate the evaporative VOC emissions.
The evaporative VOC emission, E<sub>vap</sub> (g), is calculated in CARS and shown in Eq. 8. In this equation, n is the individual vehicle belong the vehicle type v, and p is pollutant species. Those three different emission mechanisms are: diurnal emissions from tank (e<sub>d</sub>), hot and warm soak emission with fuel inject type (S<sub>fi</sub>), and running loss emissions (R). These emission mechanisms are discussed in the following part.


![Evaporative emissions](https://github.com/CMASCenter/CARS/blob/master/docs/User_Manual/media/Picture7.png)

## The emission redistribution process
The emission rate (tons y<sup>-1</sup>) data built by Calculate district emissions module only presents the annual emission rate for several grouped species, which cannot be directly applied in the chemical mechanisms of the CTM. Therefore, the module **Grid4AQM** can transfer the emission rate for the CTM. In the first part of **Grid4AQM**, the “Read_chemical” function read the chemical species profile and convert the pollutants to the specific model pollutants for the CTM chemical mechanism. The “Read_temporal” and “Read_griddesc” functions collect the monthly, weekly, hourly profiles, and grid cell format to assign the emissions to every hour and every grid cell in “Gridded_emis” function. Thus, the annual emission rate (tons y<sup>-1</sup>) data are distributed to hourly NetCDF format data in the **Grid4AQM** module. The Plot Figures module is for generating the spatial and temporal diagram for presenting the emission.
