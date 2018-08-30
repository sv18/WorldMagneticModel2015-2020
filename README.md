# World Magnetic Model 2015-2020

This project is a C/C++ implementation of the world magnetic model. ready to add the code to your projects. Validated with matlab WMM as well. forked https://www.ngdc.noaa.gov/geomag/WMM/soft.shtml

The input parameters and valid entries are:

Latitude -90.00 to +90.00 degrees 
Longitude -180.00 to +180.00 degrees 
Altitude referenced to the WGS 84 ellipsoid OR the Mean Sea Level (MSL)
Date 2015.0 to 2020.0

The seven magnetic components computed are:

F - Total Intensity of the geomagnetic field 
H - Horizontal Intensity of the geomagnetic field 
X - North Component of the geomagnetic field 
Y - East Component of the geomagnetic field 
Z - Vertical Component of the geomagnetic field 
I (DIP) - Geomagnetic Inclination 
D (DEC) - Geomagnetic Declination (Magnetic Variation)

make: (Linux Ubuntu 16)

git clone https://github.com/amirhosseinalikhahm/WorldMagneticModel2015-2020

cd WorldMagneticModel2015-2020

make

./WMM
