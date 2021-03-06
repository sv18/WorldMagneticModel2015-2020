/***********************************************************/
/* All rights reserved by Amir Hossein Alikhah Mishamandani*/
/******* Created for ADCS system of SSS-P1 Satellite *******/
/*********************** 2018 - 2019 ***********************/
/******************* Main ADCS Platform ********************/
/***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <list>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include "WMM.h"

int main(){
	
	geomag_vector now;
	time_t my_time;
	struct tm * timeinfo; 
	time (&my_time);
	timeinfo = localtime (&my_time);
	WorldMagneticModel(&now, 15, 15, 500, timeinfo->tm_year+1900);// Lat, Lon, Alt, Time
	printf("\nInitializing the World Magnetic Model Test:");
	printf("\nThe latitude = 15");
	printf("\nThe longitude = 15");
	printf("\nThe Altitude = 500");
	printf("\nThe year = %d",timeinfo->tm_year+1900);
	printf("\nThe World Magnetic Model Test successfull ...\n");
	printf("\nX - North Component of the geomagnetic field:	%lf	nT", now.X);
	printf("\nY - East Component of the geomagnetic field:	%lf	nT", now.Y);
	printf("\nZ - Vertical Component of the geomagnetic field:	%lf	nT", now.Z);
	printf("\nF - Total Intensity of the geomagnetic field:	%lf	nT", now.F);
	printf("\nH - Horizontal Intensity of the geomagnetic field:	%lf	nT", now.H);
	printf("\nI (DIP) - Geomagnetic Inclination:	%lf	deg	%lf	min", now.incld, now.inclm);
	printf("\nD (DEC) - Geomagnetic Declination (Magnetic Variation):	%lf	deg	%lf	min\n", now.decld, now.declm);
}
