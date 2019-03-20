/*
 * EphemTime.cpp
 *
 *  Created on: March 2019
 *      Author: Georgi Dragandjikov
 *
 *  EphemTime class implementation.
 *
 *
 *  libSatellitesEphemeris.so represents some of the ephemeris calculations I implemented in my
 *  planetarium software Power Age Sky Simulator (released between 1998 and 2003, www.powerageskysimulator.com).
 *  I mainly used the incredible "Astronomical Algorithms" by Jean Meeus, Second Edition.
 */

#include "EphemTime.h"



EphemTime::EphemTime( int y, int mon, int d, int h, int min, int s )
{
	year = y;
	month = mon;
	day = d;
	hour = h;
	minute = min;
	second = s;
	CalculateJD();
	CalculateJDE();
}

EphemTime::EphemTime( double jd_ut )
{
	jd = jd_ut;
	CalculateDate();
	CalculateJDE();
}

void EphemTime::CalculateDate()
{
	double jdTmp = jd + 0.5;
	double Z = (int)jdTmp;
	double F = jdTmp-Z;

	double A = Z;
	if( Z >= 2291161 ){
		double alpha = (int)( (Z - 1867216.25) / 36524.25 );
		A = Z + 1 + alpha - ((int)(alpha/4));
	}

	double B = A + 1524;
	double C = (int)( (B - 122.1) / 365.25 );
	double D = (int)( 365.25 * C );
	double E = (int)( (B-D) / 30.6001 );

	double res_day = B - D - ((int)(30.6001*E)) + F;
	double res_month = E - 1;
	if( E >= 14 )
		res_month = E - 13;
	double res_year = C - 4716;
	if( res_month <= 2 )
		res_year = C - 4715;

	double res_hour = (res_day - ((int)res_day))*24;
	double res_minute = (res_hour - ((int)res_hour))*60;
	double res_second = (res_minute - ((int)res_minute))*60;

	year = (int)res_year;
	month = (int)res_month;
	day = (int)res_day;
	hour = (int)res_hour;
	minute = (int)res_minute;
	second = (int)res_second;
}

void EphemTime::CalculateJD()
{
	int Y = year;
	int M = month;
	double D = day;
	D += (((double)hour)/24);
	D += (((double)minute)/(24*60));
	D += (((double)second)/(24*60*60));

	if( M <= 2 ){
		Y -= 1;
		M += 12;
	}

	int A = (int)( Y/100 );
	int B = 2 - A + ((int)( A/4 ));
	if( year <= 1582 )
		B = 0;

	jd = ((int)(365.25*(Y+4716))) +
		 ((int)(30.6001*(M+1))) +
		 D + B - 1524.5;
}

void EphemTime::CalculateJDE()
{
	double correction = 0;
	double t = (year - 2000.0) / 100.0;
	if( year < 948 ){ // before year 948
		correction = 2177 + 497*t + 44.1*t*t;
	}
	else if( (year >= 948 && year < 1620) || year >= 2000 ){ // between years 948 and 1620 and after 2000
		correction = 102 + 102*t + 25.3*t*t;
	}
	else{ // between years 1620 and 2000
		if( year  == 1999 ){
			correction = DT_correction[DT_ENTRIES-1][1];
		}
		else{
			for( int i = 0; i < DT_ENTRIES; i+=2 ){
				if( year == DT_correction[i][0] ){
					correction = DT_correction[i][1];
					break;
				}
				else if( year < DT_correction[i][0] ){
					correction = (DT_correction[i][1] + DT_correction[i-1][1]) / 2;
					break;
				}
			}
		}
	}

	if( year >= 2000 && year <= 2100 ){
		correction += (0.37 * (year - 2100));
	}

	jde = jd + correction/(60*60*24);
}


