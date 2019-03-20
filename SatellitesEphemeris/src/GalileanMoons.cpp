/*
 * GalileanMoons.cpp
 *
 *  Created on: March 2019
 *      Author: Georgi Dragandjikov
 *
 *  Calculates the position of the four great satellites of Jupiter (found by Galileo Galilei), low accuracy.
 *
 *
 *  libSatellitesEphemeris.so represents some of the ephemeris calculations I implemented in my
 *  planetarium software Power Age Sky Simulator (released between 1998 and 2003, www.powerageskysimulator.com).
 *  I mainly used the incredible "Astronomical Algorithms" by Jean Meeus, Second Edition.
 */

#include <math.h>
#include <vector>

#include "EphemDefs.h"
#include "EphemTime.h"

using namespace std;


void CalculateJupiterSatellites( EphemTime ephemTime, vector<Satellite>& satellites )
{
	double d = ephemTime.JDE() - 2451545;

	double V = 172.74 + 0.00111588*d;

	double M = 357.529 + 0.9856003*d;

	double N = 20.020 + 0.0830853*d + 0.329*sin( V/RADIAN_IN_DEG );

	double J = 66.115 + 0.9025179*d;

	double A = 1.915*sin( M/RADIAN_IN_DEG ) + 0.020*sin( (M/RADIAN_IN_DEG)*2 );

	double B = 5.555*sin( N/RADIAN_IN_DEG ) + 0.168*sin( (N/RADIAN_IN_DEG)*2 );

	double K = J + A - B;

	double R = 1.00014 - 0.01671*cos( M/RADIAN_IN_DEG ) - 0.00014*cos( (M/RADIAN_IN_DEG)*2 );

	double r = 5.20872 - 0.25208*cos( N/RADIAN_IN_DEG ) - 0.00611*cos( (N/RADIAN_IN_DEG)*2 );

	double delta = sqrt( r*r + R*R - 2*r*R*cos( K/RADIAN_IN_DEG ) );

	double phi = (asin( ( R/delta ) * sin( K/RADIAN_IN_DEG ) )) * RADIAN_IN_DEG;

	double lambda = 34.35 + 0.083091*d + 0.329*sin( V/RADIAN_IN_DEG ) + B;

	double Ds = 3.12 * sin( (lambda+42.8)/RADIAN_IN_DEG );

	double De = Ds -
				2.22*sin( phi/RADIAN_IN_DEG )*cos( (lambda+22)/RADIAN_IN_DEG ) -
				1.30*( (r-delta)/delta )*sin( (lambda-100.5)/RADIAN_IN_DEG );

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////

	double u1 = 163.8069 + 203.4058646*( d - delta/173 ) + phi - B;
	double u2 = 358.4140 + 101.2916335*( d - delta/173 ) + phi - B;
	double u3 =   5.7176 +  50.2345180*( d - delta/173 ) + phi - B;
	double u4 = 224.8092 +  21.4879800*( d - delta/173 ) + phi - B;

	// corrections of u1, u2, u3 and u4
	double G = 331.18 + 50.310482*( d - delta/173 );
	double H =  87.45 + 21.569231*( d - delta/173 );
	double u1minusu2 = u1 - u2;
	double u2minusu3 = u2 - u3;
	u1 += ( 0.473 * sin( 2*(u1minusu2/RADIAN_IN_DEG) ) );
	u2 += ( 1.065 * sin( 2*(u2minusu3/RADIAN_IN_DEG) ) );
	u3 += ( 0.165 * sin( G/RADIAN_IN_DEG ) );
	u4 += ( 0.843 * sin( H/RADIAN_IN_DEG ) );

	double r1 =  5.9057 - 0.0244 * cos( 2*(u1minusu2/RADIAN_IN_DEG) );
	double r2 =  9.3966 - 0.0882 * cos( 2*(u2minusu3/RADIAN_IN_DEG) );
	double r3 = 14.9883 - 0.0216 * cos( G/RADIAN_IN_DEG );
	double r4 = 26.3627 - 0.1939 * cos( H/RADIAN_IN_DEG );


	double x = r1 * sin( u1/RADIAN_IN_DEG );
	double y = -( r1 * cos( u1/RADIAN_IN_DEG ) * sin( De/RADIAN_IN_DEG ) );
	satellites.push_back( Satellite( eIo, x, y, 0 ) );

	x = r2 * sin( u2/RADIAN_IN_DEG );
	y = -( r2 * cos( u2/RADIAN_IN_DEG ) * sin( De/RADIAN_IN_DEG ) );
	satellites.push_back( Satellite( eEuropa, x, y, 0 ) );

	x = r3 * sin( u3/RADIAN_IN_DEG );
	y = -( r3 * cos( u3/RADIAN_IN_DEG ) * sin( De/RADIAN_IN_DEG ) );
	satellites.push_back( Satellite( eGanymede, x, y, 0 ) );

	x = r4 * sin( u4/RADIAN_IN_DEG );
	y = -( r4 * cos( u4/RADIAN_IN_DEG ) * sin( De/RADIAN_IN_DEG ) );
	satellites.push_back( Satellite( eCallisto, x, y, 0 ) );
}



