/*
 * Saturn.cpp
 *
 *  Created on: March 2019
 *      Author: Georgi Dragandjikov
 *
 *	GetSaturnRingPos() - calculates the position of Saturn's ring as it is observed from the Earth
 *	GetSaturnSatellites() - calculates visible positions of 7 Saturn satellites
 *
 *
 *  libSatellitesEphemeris.so represents some of the ephemeris calculations I implemented in my
 *  planetarium software Power Age Sky Simulator (released between 1998 and 2003, www.powerageskysimulator.com).
 *  I mainly used the incredible "Astronomical Algorithms" by Jean Meeus, Second Edition.
 */

#include <math.h>
#include <vector>
#include "EphemDefs.h"
#include "PeriodicTerms.h"
#include "EphemTime.h"


using namespace std;


void CalculateEarthLBR( double jd,
						double& l, double& b, double& r )
{
	double tau_earth = (jd - 2451545.0) / 365250;

	////////////// Earth L /////////////////
	double eL0 = 0;
	int ii = 0;
	for( ii = 0; ii < TERM_EARTH_L0_LEN; ++ii ){
		eL0 += term_Earth_L0[ii][0] *
			   cos( term_Earth_L0[ii][1] + term_Earth_L0[ii][2]*tau_earth );
	}
	double eL1 = 0;
	for( ii = 0; ii < TERM_EARTH_L1_LEN; ++ii ){
		eL1 += term_Earth_L1[ii][0] *
			   cos( term_Earth_L1[ii][1] + term_Earth_L1[ii][2]*tau_earth );
	}
	double eL2 = 0;
	for( ii = 0; ii < TERM_EARTH_L2_LEN; ++ii ){
		eL2 += term_Earth_L2[ii][0] *
			   cos( term_Earth_L2[ii][1] + term_Earth_L2[ii][2]*tau_earth );
	}
	double eL3 = 0;
	for( ii = 0; ii < TERM_EARTH_L3_LEN; ++ii ){
		eL3 += term_Earth_L3[ii][0] *
			   cos( term_Earth_L3[ii][1] + term_Earth_L3[ii][2]*tau_earth );
	}
	double eL4 = 0;
	for( ii = 0; ii < TERM_EARTH_L4_LEN; ++ii ){
		eL4 += term_Earth_L4[ii][0] *
			   cos( term_Earth_L4[ii][1] + term_Earth_L4[ii][2]*tau_earth );
	}
	double eL5 = 0;
	for( ii = 0; ii < TERM_EARTH_L5_LEN; ++ii ){
		eL5 += term_Earth_L5[ii][0] *
			   cos( term_Earth_L5[ii][1] + term_Earth_L5[ii][2]*tau_earth );
	}
	l = eL0 +
		 eL1*tau_earth +
		 eL2*tau_earth*tau_earth +
		 eL3*tau_earth*tau_earth*tau_earth +
		 eL4*tau_earth*tau_earth*tau_earth*tau_earth +
		 eL5*tau_earth*tau_earth*tau_earth*tau_earth*tau_earth;
	////////////////////////////////////////

	////////////// Earth B /////////////////
	double eB0 = 0;
	for( ii = 0; ii < TERM_EARTH_B0_LEN; ++ii ){
		eB0 += term_Earth_B0[ii][0] *
			   cos( term_Earth_B0[ii][1] + term_Earth_B0[ii][2]*tau_earth );
	}
	double eB1 = 0;
		for( ii = 0; ii < TERM_EARTH_B1_LEN; ++ii ){
		eB1 += term_Earth_B1[ii][0] *
			   cos( term_Earth_B1[ii][1] + term_Earth_B1[ii][2]*tau_earth );
	}
	b = eB0 +
		 eB1*tau_earth;
	////////////////////////////////////////

	////////////// Earth L /////////////////
	double eR0 = 0;
	for( ii = 0; ii < TERM_EARTH_R0_LEN; ++ii ){
		eR0 += term_Earth_R0[ii][0] *
			   cos( term_Earth_R0[ii][1] + term_Earth_R0[ii][2]*tau_earth );
	}
	double eR1 = 0;
	for( ii = 0; ii < TERM_EARTH_R1_LEN; ++ii ){
		eR1 += term_Earth_R1[ii][0] *
			   cos( term_Earth_R1[ii][1] + term_Earth_R1[ii][2]*tau_earth );
	}
	double eR2 = 0;
	for( ii = 0; ii < TERM_EARTH_R2_LEN; ++ii ){
		eR2 += term_Earth_R2[ii][0] *
			   cos( term_Earth_R2[ii][1] + term_Earth_R2[ii][2]*tau_earth );
	}
	double eR3 = 0;
	for( ii = 0; ii < TERM_EARTH_R3_LEN; ++ii ){
		eR3 += term_Earth_R3[ii][0] *
			   cos( term_Earth_R3[ii][1] + term_Earth_R3[ii][2]*tau_earth );
	}
	double eR4 = 0;
	for( ii = 0; ii < TERM_EARTH_R4_LEN; ++ii ){
		eR4 += term_Earth_R4[ii][0] *
			   cos( term_Earth_R4[ii][1] + term_Earth_R4[ii][2]*tau_earth );
	}
	r = eR0 +
		eR1*tau_earth +
		eR2*tau_earth*tau_earth +
		eR3*tau_earth*tau_earth*tau_earth +
		eR4*tau_earth*tau_earth*tau_earth*tau_earth;
	///////////////////////////////////////////
}


void CalculateSaturnLBR( double jd,
						 double& l, double& b, double& r )
{
	double tau_sat = (jd - 2451545.0) / 365250;

	int ii = 0;
	////////////// Saturn L /////////////////
	double sL0 = 0;
	for( ii = 0; ii < TERM_SAT_L0_LEN; ++ii ){
		sL0 += term_Sat_L0[ii][0] *
			   cos( term_Sat_L0[ii][1] + term_Sat_L0[ii][2]*tau_sat );
	}
	double sL1 = 0;
	for( ii = 0; ii < TERM_SAT_L1_LEN; ++ii ){
		sL1 += term_Sat_L1[ii][0] *
			   cos( term_Sat_L1[ii][1] + term_Sat_L1[ii][2]*tau_sat );
	}
	double sL2 = 0;
	for( ii = 0; ii < TERM_SAT_L2_LEN; ++ii ){
		sL2 += term_Sat_L2[ii][0] *
			   cos( term_Sat_L2[ii][1] + term_Sat_L2[ii][2]*tau_sat );
	}
	double sL3 = 0;
	for( ii = 0; ii < TERM_SAT_L3_LEN; ++ii ){
		sL3 += term_Sat_L3[ii][0] *
			   cos( term_Sat_L3[ii][1] + term_Sat_L3[ii][2]*tau_sat );
	}
	double sL4 = 0;
	for( ii = 0; ii < TERM_SAT_L4_LEN; ++ii ){
		sL4 += term_Sat_L4[ii][0] *
			   cos( term_Sat_L4[ii][1] + term_Sat_L4[ii][2]*tau_sat );
	}
	double sL5 = 0;
	for( ii = 0; ii < TERM_SAT_L5_LEN; ++ii ){
		sL5 += term_Sat_L5[ii][0] *
			   cos( term_Sat_L5[ii][1] + term_Sat_L5[ii][2]*tau_sat );
	}
	l = sL0 +
		sL1*tau_sat +
		sL2*tau_sat*tau_sat +
		sL3*tau_sat*tau_sat*tau_sat +
		sL4*tau_sat*tau_sat*tau_sat*tau_sat +
		sL5*tau_sat*tau_sat*tau_sat*tau_sat*tau_sat;
	/////////////////////////////////////////////////////////

	////////////// Saturn B /////////////////
	double sB0 = 0;
	for( ii = 0; ii < TERM_SAT_B0_LEN; ++ii ){
		sB0 += term_Sat_B0[ii][0] *
			   cos( term_Sat_B0[ii][1] + term_Sat_B0[ii][2]*tau_sat );
	}
	double sB1 = 0;
	for( ii = 0; ii < TERM_SAT_B1_LEN; ++ii ){
		sB1 += term_Sat_B1[ii][0] *
			   cos( term_Sat_B1[ii][1] + term_Sat_B1[ii][2]*tau_sat );
	}
	double sB2 = 0;
	for( ii = 0; ii < TERM_SAT_B2_LEN; ++ii ){
		sB2 += term_Sat_B2[ii][0] *
			   cos( term_Sat_B2[ii][1] + term_Sat_B2[ii][2]*tau_sat );
	}
	double sB3 = 0;
	for( ii = 0; ii < TERM_SAT_B3_LEN; ++ii ){
		sB3 += term_Sat_B3[ii][0] *
			   cos( term_Sat_B3[ii][1] + term_Sat_B3[ii][2]*tau_sat );
	}
	double sB4 = 0;
	for( ii = 0; ii < TERM_SAT_B4_LEN; ++ii ){
		sB4 += term_Sat_B4[ii][0] *
			   cos( term_Sat_B4[ii][1] + term_Sat_B4[ii][2]*tau_sat );
	}
	double sB5 = 0;
	for( ii = 0; ii < TERM_SAT_B5_LEN; ++ii ){
		sB5 += term_Sat_B5[ii][0] *
			   cos( term_Sat_B5[ii][1] + term_Sat_B5[ii][2]*tau_sat );
	}
	b = sB0 +
		sB1*tau_sat +
		sB2*tau_sat*tau_sat +
		sB3*tau_sat*tau_sat*tau_sat +
		sB4*tau_sat*tau_sat*tau_sat*tau_sat +
		sB5*tau_sat*tau_sat*tau_sat*tau_sat*tau_sat;
	/////////////////////////////////////////

	////////////// Saturn R /////////////////
	double sR0 = 0;
	for( ii = 0; ii < TERM_SAT_R0_LEN; ++ii ){
		sR0 += term_Sat_R0[ii][0] *
			   cos( term_Sat_R0[ii][1] + term_Sat_R0[ii][2]*tau_sat );
	}
	double sR1 = 0;
	for( ii = 0; ii < TERM_SAT_R1_LEN; ++ii ){
		sR1 += term_Sat_R1[ii][0] *
			   cos( term_Sat_R1[ii][1] + term_Sat_R1[ii][2]*tau_sat );
	}
	double sR2 = 0;
	for( ii = 0; ii < TERM_SAT_R2_LEN; ++ii ){
		sR2 += term_Sat_R2[ii][0] *
			   cos( term_Sat_R2[ii][1] + term_Sat_R2[ii][2]*tau_sat );
	}
	double sR3 = 0;
	for( ii = 0; ii < TERM_SAT_R3_LEN; ++ii ){
		sR3 += term_Sat_R3[ii][0] *
			   cos( term_Sat_R3[ii][1] + term_Sat_R3[ii][2]*tau_sat );
	}
	double sR4 = 0;
	for( ii = 0; ii < TERM_SAT_R4_LEN; ++ii ){
		sR4 += term_Sat_R4[ii][0] *
			   cos( term_Sat_R4[ii][1] + term_Sat_R4[ii][2]*tau_sat );
	}
	double sR5 = 0;
	for( ii = 0; ii < TERM_SAT_R5_LEN; ++ii ){
		sR5 += term_Sat_R5[ii][0] *
			   cos( term_Sat_R5[ii][1] + term_Sat_R5[ii][2]*tau_sat );
	}
	r = sR0 +
		sR1*tau_sat +
		sR2*tau_sat*tau_sat +
		sR3*tau_sat*tau_sat*tau_sat +
		sR4*tau_sat*tau_sat*tau_sat*tau_sat +
		sR5*tau_sat*tau_sat*tau_sat*tau_sat*tau_sat;
	////////////////////////////////////////
}


void CalculateSaturnDistance( double l0, double b0, double R,
							  double l, double b, double r,
							  double& x, double& y, double& z, double& delta )
{
	x = r*cos( b )*cos( l ) - R*cos( l0 );
	y = r*cos( b )*sin( l ) - R*sin( l0 );
	z = r*sin( b )          - R*sin( b0 );

	delta = sqrt( x*x + y*y + z*z );
}


void CalculateSaturnRingPosition( EphemTime ephemTime,
					   	   	      double& earthLat, double& sunLat,
								  double& posAngle,
								  double& majorAxe, double& minorAxe )
{
	// 1
	double T = (ephemTime.JDE() - 2451545.0) / 36525;
	double i = 28.075216/RADIAN_IN_DEG -
			   (0.012998/RADIAN_IN_DEG)*T +
			   (0.000004/RADIAN_IN_DEG)*T*T;
	double omi = 169.508470/RADIAN_IN_DEG +
				 (1.394681/RADIAN_IN_DEG)*T +
				 (0.000412/RADIAN_IN_DEG)*T*T;


	// 2; heliocentric longitude, latitude and radius vector of the Earth
	double l0 = 0, b0 = 0, R = 0;
	CalculateEarthLBR( ephemTime.JDE(), l0, b0, R );


	// 3; heliocentric longitude, latitude and radius vector of Saturn

	// Since the exact distance to Saturn is still not calculated,
	// we take 9 AU as an average for the time correction.
	// So, the correction of the time will be 9*AU_LIGHT_TIME
	double l = 0, b = 0, r = 0;
	CalculateSaturnLBR( ephemTime.JDE() - 9*AU_LIGHT_TIME, l, b, r );


	// 4; distance to Saturn
	double x = 0, y = 0, z = 0, delta = 0;
	CalculateSaturnDistance( l0, b0, R, l, b, r, x, y, z, delta );
	// few iterations for better precision;
	for( int ii = 0; ii < 3; ++ii ){
		CalculateSaturnLBR( ephemTime.JDE() - delta*AU_LIGHT_TIME, l, b, r );
		CalculateSaturnDistance( l0, b0, R, l, b, r, x, y, z, delta );
	}


	// 5; geocentric longitude and latitude
	double lambda = atan2( y, x );
	double tmp = sqrt( x*x + y*y );
	if( tmp == 0 )
		tmp = 0.00000001;
	double beta = atan( z / tmp );


	// 6; major and minor axes of the ring
	tmp = sin( i )*cos( beta )*sin(lambda-omi) -
		  cos( i )*sin( beta );
	if( tmp > 1 )
		tmp = 1;
	if( tmp < (-1) )
		tmp = (-1);
	double B = asin( tmp );
	double a_r = 375.35/delta; // major axis in "
	tmp = B;
	if( B < 0 )
		tmp *= (-1);
	double b_r = a_r * sin( tmp ); // minor axis in "


	// 7; longitude of the ascending node of Saturn's orbit;
	// correcting l and b for the Sun's aberration as seen from Saturn;
	double N = 113.6655/RADIAN_IN_DEG + (0.8771/RADIAN_IN_DEG)*T;
	double l_pr = l - (0.01759/RADIAN_IN_DEG)/r;
	double b_pr = b - (0.000764/RADIAN_IN_DEG)*(cos( l-N )/r);


	// 8
	tmp = sin( i )*cos( b_pr )*sin( l_pr-omi ) - cos( i )*sin( b_pr );
	if( tmp > 1 )
		tmp = 1;
	if( tmp < (-1) )
		tmp = (-1);
	double B_pr = asin( tmp );


	// 9


	// 10; nutation in longitude and nutation in obliquity
	double D_nut = 297.85036 +
				   445267.111480*T -
				   0.0019142*T*T +
				   (T*T*T)/189474;
	D_nut /= RADIAN_IN_DEG;
	double M_nut = 357.52772 +
				   35999.050340*T -
				   0.0001603*T*T -
				   (T*T*T)/300000;
	M_nut /= RADIAN_IN_DEG;
	double Mpr_nut = 134.96298 +
					 477198.867398*T +
					 0.0086972*T*T +
					 (T*T*T)/56250;
	Mpr_nut /= RADIAN_IN_DEG;
	double F_nut = 93.27191 +
				   483202.017538*T -
				   0.0036825*T*T +
				   (T*T*T)/327270;
	F_nut /= RADIAN_IN_DEG;
	double OMI_nut = 125.04452 -
					 1934.136261*T +
					 0.0020708*T*T +
					 (T*T*T)/450000;
	OMI_nut /= RADIAN_IN_DEG;

	double delta_phi = 0;
	double delta_eps = 0;

	delta_phi += ( (-171996-174.2*T) * sin( OMI_nut ) );
	delta_phi += ( (-13187-1.6*T) * sin( -2*D_nut + 2*F_nut + 2*OMI_nut ) );
	delta_phi += ( (-2274-0.2*T) * sin( 2*F_nut + 2*OMI_nut ) );
	delta_phi /= 10000;
	delta_phi /= 3600;
	delta_phi /= RADIAN_IN_DEG;

	delta_eps += ( (92025+8.9*T) * cos( OMI_nut ) );
	delta_eps += ( (5736-3.1*T) * cos( -2*D_nut + 2*F_nut + 2*OMI_nut ) );
	delta_eps += ( (977-0.5*T) * cos( 2*F_nut + 2*OMI_nut ) );
	delta_eps /= 10000;
	delta_eps /= 3600;
	delta_eps /= RADIAN_IN_DEG;

	double eps_mean = 84381.448 - 46.8150*T - 0.00059*T*T + 0.001813*T*T*T;
	eps_mean /= 3600;
	eps_mean /= RADIAN_IN_DEG;


	// 11; ecliptical longitude and latitude of the northern pole of the ring plane
	double lambda0 = omi - (90/RADIAN_IN_DEG);
	double beta0 = (90/RADIAN_IN_DEG) - i;


	// 12; correction of the geocentric longitude and latitude for the aberration of Saturn
	double corr_lambda = 0.005693 * ( cos( l0-lambda ) / cos( beta ) );
	corr_lambda /= RADIAN_IN_DEG;
	lambda += corr_lambda;
	double corr_beta = 0.005693 * ( sin( l0-lambda ) * sin( beta ) );
	corr_beta /= RADIAN_IN_DEG;
	beta += corr_beta;


	// 13
	lambda0 += delta_phi;
	lambda += delta_phi;


	// 14; transforming from ecliptical to equatorial coordinates
	double alpha0 = atan( (sin(lambda0)*cos(eps_mean) - tan(beta0)*sin(eps_mean)) /
					( cos(lambda0) ) );
	double delt0 = asin( sin(beta0)*cos(eps_mean) + cos(beta0)*sin(eps_mean)*sin(lambda0) );

	double alpha = atan2( (sin(lambda)*cos(eps_mean) - tan(beta)*sin(eps_mean)) ,
				   ( cos(lambda) ) );
	double delt = asin( sin(beta)*cos(eps_mean) + cos(beta)*sin(eps_mean)*sin(lambda) );


	// 15; the position angle
	double tmp1 = cos(delt0) * sin(alpha0-alpha);
	double tmp2 = sin(delt0)*cos(delt) - cos(delt0)*sin(delt)*cos(alpha0-alpha);
	double P = atan2( tmp1, tmp2 );


	earthLat = B*RADIAN_IN_DEG;		// degrees
	sunLat = B_pr*RADIAN_IN_DEG;	// degrees
	posAngle = P*RADIAN_IN_DEG;		// degrees
	majorAxe = a_r;					// arcseconds
	minorAxe = b_r;					// arcseconds

}


void SatelliteFinalCorrection( ESatellite satellite, double r, double delta, double Z,
							   double& X, double& Y )
{
	double K = 0;
	switch( satellite ){
		case eMimas:
			K = 20947;
			break;
		case eEnceladus:
			K = 23715;
			break;
		case eTethys:
			K = 26382;
			break;
		case eDione:
			K = 29876;
			break;
		case eRhea:
			K = 35313;
			break;
		case eTitan:
			K = 53800;
			break;
		case eHyperion:
			K = 59222;
			break;
		default:
			break;
	}
	double mod_Z = Z;
	if( mod_Z < 0 )
		mod_Z *= (-1);
	double tmp1 = X / r;
	double X_corr = (mod_Z/K) * sqrt( 1 - tmp1*tmp1 );
	X += X_corr;

	double W = delta / (delta + Z/2475);
	X *= W;
	Y *= W;
}


void CalculateSaturnSatellites( EphemTime ephemTime, vector<Satellite>& satellites )
{
	// heliocentric longitude, latitude and radius vector of the Earth
	double l0 = 0, b0 = 0, R = 0;
	CalculateEarthLBR( ephemTime.JDE(), l0, b0, R );


	// heliocentric longitude, latitude and radius vector of Saturn
	// Since the exact distance to Saturn is still not calculated,
	// we take 9 AU as an average for the time correction.
	// So, the correction of the time will be 9*AU_LIGHT_TIME
	double l = 0, b = 0, r = 0;
	CalculateSaturnLBR( ephemTime.JDE() - 9*AU_LIGHT_TIME, l, b, r );


	// distance to Saturn
	double x = 0, y = 0, z = 0, delta = 0;
	CalculateSaturnDistance( l0, b0, R, l, b, r, x, y, z, delta );
	// few iterations for better precision;
	for( int ii = 0; ii < 3; ++ii ){
		CalculateSaturnLBR( ephemTime.JDE() - delta*AU_LIGHT_TIME, l, b, r );
		CalculateSaturnDistance( l0, b0, R, l, b, r, x, y, z, delta );
	}


	// Saturn's geocentric longitude and latitude
	double lambda = atan2( y, x );
	double tmp = sqrt( x*x + y*y );
	if( tmp == 0 )
		tmp = 0.00000001;
	double beta = atan( z / tmp );


	// transforming Saturn's geocentric longitude and latitude to equinox B1950.0
	double sat_t = (2433282.4235 - ephemTime.JDE()) / 36525;
	double sat_T = 0; // J2000
	double eta = (47.0029 - 0.06603*sat_T + 0.000598*sat_T*sat_T)*sat_t +
					 ((-0.03302) + 0.000598*sat_T)*sat_t*sat_t +
					 0.000060*sat_t*sat_t*sat_t;
	eta /= 3600;
	eta /= RADIAN_IN_DEG;

	double pi = /*174.876384 degrees*/629554.9824 + 3289.4789*sat_T + 0.60622*sat_T*sat_T -
							(869.8089 + 0.50491*sat_T)*sat_t + 0.03536*sat_t*sat_t;
	pi /= 3600;
	pi /= RADIAN_IN_DEG;

	double rho = (5029.0966 + 2.22226*sat_T - 0.000042*sat_T*sat_T)*sat_t +
					 (1.11113 - 0.000042*sat_T)*sat_t*sat_t -
					 0.000006*sat_t*sat_t*sat_t;
	rho /= 3600;
	rho /= RADIAN_IN_DEG;

	double Apr = cos(eta)*cos(beta)*sin( pi-lambda ) -
						 sin(eta)*sin(beta);
	double Bpr = cos(beta)*cos( pi-lambda );
	double Cpr = cos(eta)*sin(beta) +
						 sin(eta)*cos(beta)*sin( pi-lambda );

	double lambda1950 = rho + pi - atan2( Apr, Bpr );
	double beta1950 = asin( Cpr );



	double JDE = ephemTime.JDE() - AU_LIGHT_TIME*delta;
	double t1 = JDE - 2411093.0;
	double t2 = t1 / 365.25;
	double t3 = (JDE - 2433282.423)/365.25 + 1950.0;
	double t4 = JDE - 2411368.0;
	double t5 = t4 / 365.25;
	double t6 = JDE - 2415020.0;
	double t7 = t6/36525;
	double t8 = t6/365.25;
	double t9 = (JDE - 2442000.5)/365.25;
//	double t10 = JDE - 2409786.5;
//	double t11 = t10/36525;


	double W0 = 5.095*( t3 - 1866.39 );
	double W1 = 74.4 + 32.39*t2;
	double W2 = 134.3 + 92.62*t2;
	double W3 = 42.0 - 0.5118*t5;
	double W4 = 276.59 + 0.5118*t5;
	double W5 = 267.2635 + 1222.1136*t7;
	double W6 = 175.4762 + 1221.5515*t7;
	double W7 = 2.4891 + 0.002435*t7;
	double W8 = 113.35 - 0.2597*t7;
	W0 /= RADIAN_IN_DEG;
	W1 /= RADIAN_IN_DEG;
	W2 /= RADIAN_IN_DEG;
	W3 /= RADIAN_IN_DEG;
	W4 /= RADIAN_IN_DEG;
	W5 /= RADIAN_IN_DEG;
	W6 /= RADIAN_IN_DEG;
	W7 /= RADIAN_IN_DEG;
	W8 /= RADIAN_IN_DEG;


	double s1 = sin( 28.0817/RADIAN_IN_DEG );
	double c1 = cos( 28.0817/RADIAN_IN_DEG );
	double s2 = sin( 168.8112/RADIAN_IN_DEG );
	double c2 = cos( 168.8112/RADIAN_IN_DEG );
	double e1 = 0.05589 - 0.000346*t7;


	// the ninth fictitious satellite //////////
	double fict_X = 0;
	double fict_Y = 0;
	double fict_Z = 1;

	double fict_A1 = fict_X;
	double fict_B1 = c1*fict_Y - s1*fict_Z;
	double fict_C1 = s1*fict_Y + c1*fict_Z;

	double fict_A2 = c2*fict_A1 - s2*fict_B1;
	double fict_B2 = s2*fict_A1 + c2*fict_B1;
	double fict_C2 = fict_C1;

	double fict_A3 = fict_A2*sin(lambda1950) - fict_B2*cos(lambda1950);
	double fict_B3 = fict_A2*cos(lambda1950) + fict_B2*sin(lambda1950);
	double fict_C3 = fict_C2 = fict_C1;

	double fict_A4 = fict_A3;
//	double fict_B4 = fict_B3*cos(beta1950) + fict_C3*sin(beta1950);
	double fict_C4 = fict_C3*cos(beta1950) - fict_B3*sin(beta1950);

	double fict_D = atan2( fict_A4, fict_C4 );
	////////////////////////////////////////////


	// Mimas (I) ///////////////////////////
	double mimas_L = 127.64 + 381.994497*t1 -
					 43.57*sin(W0) -
					 0.720*sin(3*W0) -
					 0.02144*sin(5*W0);
	mimas_L /= RADIAN_IN_DEG;

	double mimas_p = 106.1 + 365.549*t2;
	mimas_p /= RADIAN_IN_DEG;

	double mimas_M = mimas_L - mimas_p;

	double mimas_C = 2.18287*sin(mimas_M) +
					 0.025988*sin(2*mimas_M) +
					 0.00043*sin(3*mimas_M);
	mimas_C /= RADIAN_IN_DEG;

	double mimas_lambda = mimas_L + mimas_C;
	double mimas_r = 3.06879 / ( 1 + 0.01905*cos( mimas_M+mimas_C ) );
	double mimas_gama = 1.563/RADIAN_IN_DEG;
	double mimas_omi = 54.5 - 365.072*t2;
	mimas_omi /= RADIAN_IN_DEG;

	double mimas_u = mimas_lambda - mimas_omi;
	double mimas_w = mimas_omi - 168.8112/RADIAN_IN_DEG;

	double mimas_X = mimas_r*( cos(mimas_u)*cos(mimas_w) -
								   sin(mimas_u)*cos(mimas_gama)*sin(mimas_w) );
	double mimas_Y = mimas_r*( sin(mimas_u)*cos(mimas_w)*cos(mimas_gama) +
								   cos(mimas_u)*sin(mimas_w) );
	double mimas_Z = mimas_r*sin(mimas_u)*sin(mimas_gama);

	double mimas_A1 = mimas_X;
	double mimas_B1 = c1*mimas_Y - s1*mimas_Z;
	double mimas_C1 = s1*mimas_Y + c1*mimas_Z;

	double mimas_A2 = c2*mimas_A1 - s2*mimas_B1;
	double mimas_B2 = s2*mimas_A1 + c2*mimas_B1;
	double mimas_C2 = mimas_C1;

	double mimas_A3 = mimas_A2*sin(lambda1950) - mimas_B2*cos(lambda1950);
	double mimas_B3 = mimas_A2*cos(lambda1950) + mimas_B2*sin(lambda1950);
	double mimas_C3 = mimas_C2;

	double mimas_A4 = mimas_A3;
	double mimas_B4 = mimas_B3*cos(beta1950) + mimas_C3*sin(beta1950);
	double mimas_C4 = mimas_C3*cos(beta1950) - mimas_B3*sin(beta1950);

	double mimas_fin_X = mimas_A4*cos(fict_D) - mimas_C4*sin(fict_D);
	double mimas_fin_Y = mimas_A4*sin(fict_D) + mimas_C4*cos(fict_D);
	double mimas_fin_Z = mimas_B4;

	SatelliteFinalCorrection( eMimas, mimas_r, delta, mimas_fin_Z,
							  mimas_fin_X, mimas_fin_Y );

	satellites.push_back( Satellite( eMimas, mimas_fin_X, mimas_fin_Y, mimas_fin_Z ) );
	////////////////////////////////////////


	// Enceladus (II) //////////////////////
	double encel_L = 200.317 +
					 262.7319002*t1 +
					 0.25667*sin(W1) +
					 0.20883*sin(W2);
	encel_L /= RADIAN_IN_DEG;

	double encel_p = 309.107 + 123.44121*t2;
	encel_p /= RADIAN_IN_DEG;

	double encel_M = encel_L - encel_p;

	double encel_C = 0.55577*sin(encel_M) + 0.00168*sin(2*encel_M);
	encel_C /= RADIAN_IN_DEG;

	double encel_lambda = encel_L + encel_C;
	double encel_r = 3.94118 / ( 1 + 0.00485*cos( encel_M+encel_C ) );
	double encel_gama = 0.0262/RADIAN_IN_DEG;
	double encel_omi = 348 - 151.95*t2;
	encel_omi /= RADIAN_IN_DEG;

	double encel_u = encel_lambda - encel_omi;
	double encel_w = encel_omi - 168.8112/RADIAN_IN_DEG;

	double encel_X = encel_r*( cos(encel_u)*cos(encel_w) -
							   sin(encel_u)*cos(encel_gama)*sin(encel_w) );
	double encel_Y = encel_r*( sin(encel_u)*cos(encel_w)*cos(encel_gama) +
							   cos(encel_u)*sin(encel_w) );
	double encel_Z = encel_r*sin(encel_u)*sin(encel_gama);

	double encel_A1 = encel_X;
	double encel_B1 = c1*encel_Y - s1*encel_Z;
	double encel_C1 = s1*encel_Y + c1*encel_Z;

	double encel_A2 = c2*encel_A1 - s2*encel_B1;
	double encel_B2 = s2*encel_A1 + c2*encel_B1;
	double encel_C2 = encel_C1;

	double encel_A3 = encel_A2*sin(lambda1950) - encel_B2*cos(lambda1950);
	double encel_B3 = encel_A2*cos(lambda1950) + encel_B2*sin(lambda1950);
	double encel_C3 = encel_C2;

	double encel_A4 = encel_A3;
	double encel_B4 = encel_B3*cos(beta1950) + encel_C3*sin(beta1950);
	double encel_C4 = encel_C3*cos(beta1950) - encel_B3*sin(beta1950);

	double encel_fin_X = encel_A4*cos(fict_D) - encel_C4*sin(fict_D);
	double encel_fin_Y = encel_A4*sin(fict_D) + encel_C4*cos(fict_D);
	double encel_fin_Z = encel_B4;

	SatelliteFinalCorrection( eEnceladus, encel_r, delta, encel_fin_Z,
							  encel_fin_X, encel_fin_Y );

	satellites.push_back( Satellite( eEnceladus, encel_fin_X, encel_fin_Y, encel_fin_Z ) );
	////////////////////////////////////////


	// Tethys (III) ////////////////////////
	double tethys_lambda = 285.306 +
						   190.69791226*t1 +
						   2.063*sin(W0) +
						   0.03409*sin(3*W0) +
						   0.001015*sin(5*W0);
	tethys_lambda /= RADIAN_IN_DEG;
	double tethys_r = 4.880998;
	double tethys_gama = 1.0976/RADIAN_IN_DEG;
	double tethys_omi = 111.33 - 72.2441*t2;
	tethys_omi /= RADIAN_IN_DEG;

	double tethys_u = tethys_lambda - tethys_omi;
	double tethys_w = tethys_omi - 168.8112/RADIAN_IN_DEG;

	double tethys_X = tethys_r*( cos(tethys_u)*cos(tethys_w) -
								 sin(tethys_u)*cos(tethys_gama)*sin(tethys_w) );
	double tethys_Y = tethys_r*( sin(tethys_u)*cos(tethys_w)*cos(tethys_gama) +
								 cos(tethys_u)*sin(tethys_w) );
	double tethys_Z = tethys_r*sin(tethys_u)*sin(tethys_gama);

	double tethys_A1 = tethys_X;
	double tethys_B1 = c1*tethys_Y - s1*tethys_Z;
	double tethys_C1 = s1*tethys_Y + c1*tethys_Z;

	double tethys_A2 = c2*tethys_A1 - s2*tethys_B1;
	double tethys_B2 = s2*tethys_A1 + c2*tethys_B1;
	double tethys_C2 = tethys_C1;

	double tethys_A3 = tethys_A2*sin(lambda1950) - tethys_B2*cos(lambda1950);
	double tethys_B3 = tethys_A2*cos(lambda1950) + tethys_B2*sin(lambda1950);
	double tethys_C3 = tethys_C2;

	double tethys_A4 = tethys_A3;
	double tethys_B4 = tethys_B3*cos(beta1950) + tethys_C3*sin(beta1950);
	double tethys_C4 = tethys_C3*cos(beta1950) - tethys_B3*sin(beta1950);

	double tethys_fin_X = tethys_A4*cos(fict_D) - tethys_C4*sin(fict_D);
	double tethys_fin_Y = tethys_A4*sin(fict_D) + tethys_C4*cos(fict_D);
	double tethys_fin_Z = tethys_B4;

	SatelliteFinalCorrection( eTethys, tethys_r, delta, tethys_fin_Z,
							  tethys_fin_X, tethys_fin_Y );

	satellites.push_back( Satellite( eTethys, tethys_fin_X, tethys_fin_Y, tethys_fin_Z ) );
	////////////////////////////////////////


	// dione (IV) //////////////////////////
	double dione_L = 254.712 +
					 131.53493193*t1 -
					 0.0215*sin(W1) -
					 0.01733*sin(W2);
	dione_L /= RADIAN_IN_DEG;
	double dione_p = 174.8 + 30.820*t2;
	dione_p /= RADIAN_IN_DEG;
	double dione_M = dione_L - dione_p;
	double dione_C = 0.24717*sin(dione_M) + 0.00033*sin(2*dione_M);
	dione_C /= RADIAN_IN_DEG;

	double dione_lambda = dione_L + dione_C;
	double dione_r = 6.24871 / (1 + 0.002157*cos( dione_M+dione_C ) );
	double dione_gama = 0.0139/RADIAN_IN_DEG;
	double dione_omi = 232 - 30.27*t2;
	dione_omi /= RADIAN_IN_DEG;

	double dione_u = dione_lambda - dione_omi;
	double dione_w = dione_omi - 168.8112/RADIAN_IN_DEG;

	double dione_X = dione_r*( cos(dione_u)*cos(dione_w) -
							   sin(dione_u)*cos(dione_gama)*sin(dione_w) );
	double dione_Y = dione_r*( sin(dione_u)*cos(dione_w)*cos(dione_gama) +
							   cos(dione_u)*sin(dione_w) );
	double dione_Z = dione_r*sin(dione_u)*sin(dione_gama);

	double dione_A1 = dione_X;
	double dione_B1 = c1*dione_Y - s1*dione_Z;
	double dione_C1 = s1*dione_Y + c1*dione_Z;

	double dione_A2 = c2*dione_A1 - s2*dione_B1;
	double dione_B2 = s2*dione_A1 + c2*dione_B1;
	double dione_C2 = dione_C1;

	double dione_A3 = dione_A2*sin(lambda1950) - dione_B2*cos(lambda1950);
	double dione_B3 = dione_A2*cos(lambda1950) + dione_B2*sin(lambda1950);
	double dione_C3 = dione_C2;

	double dione_A4 = dione_A3;
	double dione_B4 = dione_B3*cos(beta1950) + dione_C3*sin(beta1950);
	double dione_C4 = dione_C3*cos(beta1950) - dione_B3*sin(beta1950);

	double dione_fin_X = dione_A4*cos(fict_D) - dione_C4*sin(fict_D);
	double dione_fin_Y = dione_A4*sin(fict_D) + dione_C4*cos(fict_D);
	double dione_fin_Z = dione_B4;

	SatelliteFinalCorrection( eDione, dione_r, delta, dione_fin_Z,
							  dione_fin_X, dione_fin_Y );

	satellites.push_back( Satellite( eDione, dione_fin_X, dione_fin_Y, dione_fin_Z ) );
	////////////////////////////////////////

	// Rhea (V) ////////////////////////////
	double rhea_p_pr = 342.7 + 10.057*t2;
	rhea_p_pr /= RADIAN_IN_DEG;
	double rhea_a1 = 0.000265*sin(rhea_p_pr) + 0.01*sin(W4);
	double rhea_a2 = 0.000265*cos(rhea_p_pr) + 0.01*cos(W4);
	double rhea_e = sqrt( rhea_a1*rhea_a1 + rhea_a2*rhea_a2 );

	double rhea_p = atan( rhea_a1/rhea_a2 );
	if( rhea_a2 < 0 )
		rhea_p += PI_VAL;

	double rhea_N = 345 - 10.057*t2;
	rhea_N /= RADIAN_IN_DEG;
	double rhea_lambda_pr = 359.244 + 79.69004720*t1 + 0.086754*sin(rhea_N);
	rhea_lambda_pr /= RADIAN_IN_DEG;
	double rhea_i = 28.0362 + 0.346898*cos(rhea_N) + 0.01930*cos(W3);
	rhea_i /= RADIAN_IN_DEG;
	double rhea_omi_prev = 168.8034 + 0.736936*sin(rhea_N) + 0.041*sin(W3);
	rhea_omi_prev /= RADIAN_IN_DEG;
	double rhea_a = 8.725924/RADIAN_IN_DEG;


	double _rhea_M = rhea_lambda_pr - rhea_p;
	double _rhea_C = (2*rhea_e - 0.25*rhea_e*rhea_e*rhea_e +
					 0.0520833333*rhea_e*rhea_e*rhea_e*rhea_e*rhea_e)*sin(_rhea_M) +
					 (1.25*rhea_e*rhea_e -
					 0.458333333*rhea_e*rhea_e*rhea_e*rhea_e)*sin(2*_rhea_M) +
					 (1.083333333*rhea_e*rhea_e*rhea_e -
					 0.671875*rhea_e*rhea_e*rhea_e*rhea_e*rhea_e)*sin(3*_rhea_M) +
					 1.072917*rhea_e*rhea_e*rhea_e*rhea_e*sin(4*_rhea_M) +
					 1.142708*rhea_e*rhea_e*rhea_e*rhea_e*rhea_e*sin(5*_rhea_M);
	double _rhea_r = (rhea_a*(1 - rhea_e*rhea_e)) / (1 + rhea_e*cos( _rhea_M+_rhea_C ));
	double _rhea_g = rhea_omi_prev - 168.8112/RADIAN_IN_DEG;
	double _rhea_a1 = sin(rhea_i)*sin(_rhea_g);
	double _rhea_a2 = c1*sin(rhea_i)*cos(_rhea_g) - s1*cos(rhea_i);
	double _rhea_gama = asin( sqrt(_rhea_a1*_rhea_a1 + _rhea_a2*_rhea_a2) );

	double _rhea_u = atan( _rhea_a1/_rhea_a2 );
	if( _rhea_a2 < 0 )
		_rhea_u += PI_VAL;

	double _rhea_w = 168.8112/RADIAN_IN_DEG + _rhea_u;
	double _rhea_h = c1*sin(rhea_i) - s1*cos(rhea_i)*cos(_rhea_g);

	double _rhea_phi = atan( (s1*sin(_rhea_g)) / _rhea_h );
	if( _rhea_h < 0 )
		_rhea_phi += PI_VAL;

	double _rhea_lambda = rhea_lambda_pr + _rhea_C + _rhea_u - _rhea_g - _rhea_phi;

	double rhea_lambda = _rhea_lambda;
	double rhea_gama = _rhea_gama;
	double rhea_omi = _rhea_w;
	double rhea_r = _rhea_r;
	rhea_r *= RADIAN_IN_DEG;

	double rhea_u = rhea_lambda - rhea_omi;
	double rhea_w = rhea_omi - 168.8112/RADIAN_IN_DEG;

	double rhea_X = rhea_r*( cos(rhea_u)*cos(rhea_w) -
							 sin(rhea_u)*cos(rhea_gama)*sin(rhea_w) );
	double rhea_Y = rhea_r*( sin(rhea_u)*cos(rhea_w)*cos(rhea_gama) +
							 cos(rhea_u)*sin(rhea_w) );
	double rhea_Z = rhea_r*sin(rhea_u)*sin(rhea_gama);

	double rhea_A1 = rhea_X;
	double rhea_B1 = c1*rhea_Y - s1*rhea_Z;
	double rhea_C1 = s1*rhea_Y + c1*rhea_Z;

	double rhea_A2 = c2*rhea_A1 - s2*rhea_B1;
	double rhea_B2 = s2*rhea_A1 + c2*rhea_B1;
	double rhea_C2 = rhea_C1;

	double rhea_A3 = rhea_A2*sin(lambda1950) - rhea_B2*cos(lambda1950);
	double rhea_B3 = rhea_A2*cos(lambda1950) + rhea_B2*sin(lambda1950);
	double rhea_C3 = rhea_C2;

	double rhea_A4 = rhea_A3;
	double rhea_B4 = rhea_B3*cos(beta1950) + rhea_C3*sin(beta1950);
	double rhea_C4 = rhea_C3*cos(beta1950) - rhea_B3*sin(beta1950);

	double rhea_fin_X = rhea_A4*cos(fict_D) - rhea_C4*sin(fict_D);
	double rhea_fin_Y = rhea_A4*sin(fict_D) + rhea_C4*cos(fict_D);
	double rhea_fin_Z = rhea_B4;

	SatelliteFinalCorrection( eRhea, rhea_r, delta, rhea_fin_Z,
							  rhea_fin_X, rhea_fin_Y );

	satellites.push_back( Satellite( eRhea, rhea_fin_X, rhea_fin_Y, rhea_fin_Z ) );
	////////////////////////////////////////


	// Titan (VI) ////////////////////////////
	double titan_L = 261.1582 + 22.57697855*t4 + 0.074025*sin(W3);
	titan_L /= RADIAN_IN_DEG;
	double titan_i_pr = 27.45141 + 0.295999*cos(W3);
	titan_i_pr /= RADIAN_IN_DEG;
	double titan_omi_pr = 168.66925 + 0.628808*sin(W3);
	titan_omi_pr /= RADIAN_IN_DEG;
	double titan_a1 = sin(W7)*sin( titan_omi_pr-W8 );
	double titan_a2 = cos(W7)*sin(titan_i_pr) -
					  sin(W7)*cos(titan_i_pr)*cos( titan_omi_pr-W8 );
	double titan_g0 = 102.8623/RADIAN_IN_DEG;

	double titan_phi = atan( titan_a1/titan_a2 );
	if( titan_a2 < 0 )
		titan_phi += PI_VAL;

	double titan_s = sqrt( titan_a1*titan_a1 + titan_a2*titan_a2 );
	double titan_g = W4 - titan_omi_pr - titan_phi;

	double titan_omega;
	for( int titan_ii = 0; titan_ii < 3; titan_ii++ ){
		titan_omega = W4 +
					  (0.37515/RADIAN_IN_DEG)*( sin(2*titan_g) -
					  sin(2*titan_g0) );
		titan_g = titan_omega - titan_omi_pr - titan_phi;
	}

	double titan_e_pr = 0.029092 + 0.00019048*( cos(2*titan_g) - cos(2*titan_g0) );
	double titan_q = 2*( W5 - titan_omega );
	double titan_b1 = sin(titan_i_pr)*sin( titan_omi_pr - W8 );
	double titan_b2 = cos(W7)*sin(titan_i_pr)*cos( titan_omi_pr-W8 ) -
					  sin(W7)*cos(titan_i_pr);

	double titan_tita = atan( titan_b1/titan_b2 );
	if( titan_b2 < 0 )
		titan_tita += PI_VAL;
	titan_tita += W8;

	double titan_e = titan_e_pr + 0.002778797*titan_e_pr*cos(titan_q);
	double titan_p = titan_omega + (0.159215/RADIAN_IN_DEG)*sin(titan_q);
	double titan_u_prev = 2*W5 - 2*titan_tita + titan_phi;
	double titan_h = 0.9375*titan_e_pr*titan_e_pr*sin(titan_q) +
					 0.1875*titan_s*titan_s*sin( 2*(W5-titan_tita) );
	double titan_lambda_pr = titan_L -
							 (0.254744/RADIAN_IN_DEG)*
							 (e1*sin(W6) + 0.75*e1*e1*sin(2*W6) + titan_h);
	double titan_i = titan_i_pr + (0.031843/RADIAN_IN_DEG)*titan_s*cos(titan_u_prev);
	double titan_omi_prev = titan_omi_pr +
							((0.031843/RADIAN_IN_DEG)*titan_s*sin(titan_u_prev))/sin(titan_i_pr);
	double titan_a = 20.216193/RADIAN_IN_DEG;


	double _titan_M = titan_lambda_pr - titan_p;
	double _titan_C = (2*titan_e - 0.25*titan_e*titan_e*titan_e +
					  0.0520833333*titan_e*titan_e*titan_e*titan_e*titan_e)*sin(_titan_M) +
					  (1.25*titan_e*titan_e -
					  0.458333333*titan_e*titan_e*titan_e*titan_e)*sin(2*_titan_M) +
					  (1.083333333*titan_e*titan_e*titan_e -
					  0.671875*titan_e*titan_e*titan_e*titan_e*titan_e)*sin(3*_titan_M) +
					  1.072917*titan_e*titan_e*titan_e*titan_e*sin(4*_titan_M) +
					  1.142708*titan_e*titan_e*titan_e*titan_e*titan_e*sin(5*_titan_M);
	double _titan_r = (titan_a*(1 - titan_e*titan_e)) / (1 + titan_e*cos( _titan_M+_titan_C ));
	double _titan_g = titan_omi_prev - 168.8112/RADIAN_IN_DEG;
	double _titan_a1 = sin(titan_i)*sin(_titan_g);
	double _titan_a2 = c1*sin(titan_i)*cos(_titan_g) - s1*cos(titan_i);
	double _titan_gama = asin( sqrt(_titan_a1*_titan_a1 + _titan_a2*_titan_a2) );

	double _titan_u = atan( _titan_a1/_titan_a2 );
	if( _titan_a2 < 0 )
		_titan_u += PI_VAL;

	double _titan_w = 168.8112/RADIAN_IN_DEG + _titan_u;
	double _titan_h = c1*sin(titan_i) - s1*cos(titan_i)*cos(_titan_g);

	double _titan_phi = atan( (s1*sin(_titan_g)) / _titan_h );
	if( _titan_h < 0 )
		_titan_phi += PI_VAL;

	double _titan_lambda = titan_lambda_pr + _titan_C + _titan_u - _titan_g - _titan_phi;

	double titan_lambda = _titan_lambda;
	double titan_gama = _titan_gama;
	double titan_omi = _titan_w;
	double titan_r = _titan_r;
	titan_r *= RADIAN_IN_DEG;


	double titan_u = titan_lambda - titan_omi;
	double titan_w = titan_omi - 168.8112/RADIAN_IN_DEG;

	double titan_X = titan_r*( cos(titan_u)*cos(titan_w) -
							 sin(titan_u)*cos(titan_gama)*sin(titan_w) );
	double titan_Y = titan_r*( sin(titan_u)*cos(titan_w)*cos(titan_gama) +
							 cos(titan_u)*sin(titan_w) );
	double titan_Z = titan_r*sin(titan_u)*sin(titan_gama);

	double titan_A1 = titan_X;
	double titan_B1 = c1*titan_Y - s1*titan_Z;
	double titan_C1 = s1*titan_Y + c1*titan_Z;

	double titan_A2 = c2*titan_A1 - s2*titan_B1;
	double titan_B2 = s2*titan_A1 + c2*titan_B1;
	double titan_C2 = titan_C1;

	double titan_A3 = titan_A2*sin(lambda1950) - titan_B2*cos(lambda1950);
	double titan_B3 = titan_A2*cos(lambda1950) + titan_B2*sin(lambda1950);
	double titan_C3 = titan_C2;

	double titan_A4 = titan_A3;
	double titan_B4 = titan_B3*cos(beta1950) + titan_C3*sin(beta1950);
	double titan_C4 = titan_C3*cos(beta1950) - titan_B3*sin(beta1950);

	double titan_fin_X = titan_A4*cos(fict_D) - titan_C4*sin(fict_D);
	double titan_fin_Y = titan_A4*sin(fict_D) + titan_C4*cos(fict_D);
	double titan_fin_Z = titan_B4;

	SatelliteFinalCorrection( eTitan, titan_r, delta, titan_fin_Z,
							  titan_fin_X, titan_fin_Y );

	satellites.push_back( Satellite( eTitan, titan_fin_X, titan_fin_Y, titan_fin_Z ) );
	//////////////////////////////////////////


	// Hyperion (VII) ////////////////////////////
	double hyperion_eta = 92.39 + 0.5621071*t6;
	hyperion_eta /= RADIAN_IN_DEG;
	double hyperion_sigma = 148.19 - 19.18*t8;
	hyperion_sigma /= RADIAN_IN_DEG;
	double hyperion_tita = 184.8 - 35.41*t9;
	hyperion_tita /= RADIAN_IN_DEG;
	double hyperion_tita_pr = hyperion_tita - 7.5/RADIAN_IN_DEG;
	double hyperion_as = 176 + 12.22*t8;
	hyperion_as /= RADIAN_IN_DEG;
	double hyperion_bs = 8 + 24.44*t8;
	hyperion_bs /= RADIAN_IN_DEG;
	double hyperion_cs = hyperion_bs + 5/RADIAN_IN_DEG;
	double hyperion_omega = 69.898 - 18.67088*t8;
	hyperion_omega /= RADIAN_IN_DEG;
	double hyperion_fi = 2*( hyperion_omega - W5 );
	double hyperion_ksi = 94.9 - 2.292*t8;
	hyperion_ksi /= RADIAN_IN_DEG;
	double hyperion_a = 24.50601 -
						0.08686*cos(hyperion_eta) -
						0.00166*cos(hyperion_sigma+hyperion_eta) +
						0.00175*cos(hyperion_sigma-hyperion_eta);
	double hyperion_e = 0.103458 -
						0.004099*cos(hyperion_eta) -
						0.000167*cos(hyperion_sigma+hyperion_eta) +
						0.000235*cos(hyperion_sigma-hyperion_eta) +
						0.02303*cos(hyperion_sigma) -
						0.00212*cos(2*hyperion_sigma) +
						0.000151*cos(3*hyperion_sigma) +
						0.00013*cos(hyperion_fi);
	double hyperion_p = hyperion_omega*RADIAN_IN_DEG +
						0.15648*sin(hyperion_ksi) -
						0.4457*sin(hyperion_eta) -
						0.2657*sin(hyperion_sigma+hyperion_eta) -
						0.3573*sin(hyperion_sigma-hyperion_eta) -
						12.872*sin(hyperion_sigma) +
						1.668*sin(2*hyperion_sigma) -
						0.2419*sin(3*hyperion_sigma) -
						0.07*sin(hyperion_fi);
	hyperion_p /= RADIAN_IN_DEG;
	double hyperion_lambda_pr = 177.047 +
								16.91993829*t6 +
								0.15648*sin(hyperion_ksi) +
								9.142*sin(hyperion_eta) +
								0.007*sin(2*hyperion_eta) -
								0.014*sin(3*hyperion_eta) +
								0.2275*sin(hyperion_sigma+hyperion_eta) +
								0.2112*sin(hyperion_sigma-hyperion_eta) -
								0.26*sin(hyperion_sigma) -
								0.0098*sin(2*hyperion_sigma) -
								0.013*sin(hyperion_as) +
								0.017*sin(hyperion_bs) -
								0.0303*sin(hyperion_fi);
	hyperion_lambda_pr /= RADIAN_IN_DEG;
	double hyperion_i = 27.3347 +
						0.643486*cos(hyperion_ksi) +
						0.315*cos(W3) +
						0.018*cos(hyperion_tita) -
						0.018*cos(hyperion_cs);
	hyperion_i /= RADIAN_IN_DEG;
	double hyperion_omi_prev = 168.6812 +
							   1.40136*cos(hyperion_ksi) +
							   0.68599*sin(W3) -
							   0.0392*sin(hyperion_cs) +
							   0.0366*sin(hyperion_tita_pr);
	hyperion_omi_prev /= RADIAN_IN_DEG;


	double _hyperion_M = hyperion_lambda_pr - hyperion_p;
	double _hyperion_C = (2*hyperion_e - 0.25*hyperion_e*hyperion_e*hyperion_e +
						 0.0520833333*hyperion_e*hyperion_e*hyperion_e*hyperion_e*hyperion_e)*sin(_hyperion_M) +
						 (1.25*hyperion_e*hyperion_e -
						 0.458333333*hyperion_e*hyperion_e*hyperion_e*hyperion_e)*sin(2*_hyperion_M) +
						 (1.083333333*hyperion_e*hyperion_e*hyperion_e -
						 0.671875*hyperion_e*hyperion_e*hyperion_e*hyperion_e*hyperion_e)*sin(3*_hyperion_M) +
						 1.072917*hyperion_e*hyperion_e*hyperion_e*hyperion_e*sin(4*_hyperion_M) +
						 1.142708*hyperion_e*hyperion_e*hyperion_e*hyperion_e*hyperion_e*sin(5*_hyperion_M);
	double _hyperion_r = (hyperion_a*(1 - hyperion_e*hyperion_e)) / (1 + hyperion_e*cos( _hyperion_M+_hyperion_C ));
	double _hyperion_g = hyperion_omi_prev - 168.8112/RADIAN_IN_DEG;
	double _hyperion_a1 = sin(hyperion_i)*sin(_hyperion_g);
	double _hyperion_a2 = c1*sin(hyperion_i)*cos(_hyperion_g) - s1*cos(hyperion_i);
	double _hyperion_gama = asin( sqrt(_hyperion_a1*_hyperion_a1 + _hyperion_a2*_hyperion_a2) );

	double _hyperion_u = atan( _hyperion_a1/_hyperion_a2 );
	if( _hyperion_a2 < 0 )
		_hyperion_u += PI_VAL;

	double _hyperion_w = 168.8112/RADIAN_IN_DEG + _hyperion_u;
	double _hyperion_h = c1*sin(hyperion_i) - s1*cos(hyperion_i)*cos(_hyperion_g);

	double _hyperion_phi = atan( (s1*sin(_hyperion_g)) / _hyperion_h );
	if( _hyperion_h < 0 )
		_hyperion_phi += PI_VAL;

	double _hyperion_lambda = hyperion_lambda_pr + _hyperion_C + _hyperion_u - _hyperion_g - _hyperion_phi;

	double hyperion_lambda = _hyperion_lambda;
	double hyperion_gama = _hyperion_gama;
	double hyperion_omi = _hyperion_w;
	double hyperion_r = _hyperion_r;
//	hyperion_r *= RADIAN_IN_DEG;


	double hyperion_u = hyperion_lambda - hyperion_omi;
	double hyperion_w = hyperion_omi - 168.8112/RADIAN_IN_DEG;

	double hyperion_X = hyperion_r*( cos(hyperion_u)*cos(hyperion_w) -
							 sin(hyperion_u)*cos(hyperion_gama)*sin(hyperion_w) );
	double hyperion_Y = hyperion_r*( sin(hyperion_u)*cos(hyperion_w)*cos(hyperion_gama) +
							 cos(hyperion_u)*sin(hyperion_w) );
	double hyperion_Z = hyperion_r*sin(hyperion_u)*sin(hyperion_gama);

	double hyperion_A1 = hyperion_X;
	double hyperion_B1 = c1*hyperion_Y - s1*hyperion_Z;
	double hyperion_C1 = s1*hyperion_Y + c1*hyperion_Z;

	double hyperion_A2 = c2*hyperion_A1 - s2*hyperion_B1;
	double hyperion_B2 = s2*hyperion_A1 + c2*hyperion_B1;
	double hyperion_C2 = hyperion_C1;

	double hyperion_A3 = hyperion_A2*sin(lambda1950) - hyperion_B2*cos(lambda1950);
	double hyperion_B3 = hyperion_A2*cos(lambda1950) + hyperion_B2*sin(lambda1950);
	double hyperion_C3 = hyperion_C2;

	double hyperion_A4 = hyperion_A3;
	double hyperion_B4 = hyperion_B3*cos(beta1950) + hyperion_C3*sin(beta1950);
	double hyperion_C4 = hyperion_C3*cos(beta1950) - hyperion_B3*sin(beta1950);

	double hyperion_fin_X = hyperion_A4*cos(fict_D) - hyperion_C4*sin(fict_D);
	double hyperion_fin_Y = hyperion_A4*sin(fict_D) + hyperion_C4*cos(fict_D);
	double hyperion_fin_Z = hyperion_B4;

	SatelliteFinalCorrection( eHyperion, hyperion_r, delta, hyperion_fin_Z,
							  hyperion_fin_X, hyperion_fin_Y );

	satellites.push_back( Satellite( eHyperion, hyperion_fin_X, hyperion_fin_Y, hyperion_fin_Z ) );
	//////////////////////////////////////////////

}
