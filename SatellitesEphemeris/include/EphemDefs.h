/*
 * EphemDefs.h
 *
 *  Created on: March 2019
 *      Author: Georgi Dragandjikov
 *
 *  Common definitions.
 *
 *
 *  libSatellitesEphemeris.so represents some of the ephemeris calculations I implemented in my
 *  planetarium software Power Age Sky Simulator (released between 1998 and 2003, www.powerageskysimulator.com).
 *  I mainly used the incredible "Astronomical Algorithms" by Jean Meeus, Second Edition.
 */

#ifndef INCLUDE_EPHEMDEFS_H_
#define INCLUDE_EPHEMDEFS_H_

constexpr double PI_VAL = 3.141592653589793;
constexpr double RADIAN_IN_DEG = 57.295779513;
constexpr double AU_VAL = 149597870.7;
constexpr double AU_LIGHT_TIME = 0.0057755183; // time for the light to pass 1AU (in days)

enum ESatellite{
	eIo,
	eEuropa,
	eGanymede,
	eCallisto,
	eMimas,
	eEnceladus,
	eTethys,
	eDione,
	eRhea,
	eTitan,
	eHyperion
};

struct Satellite{
	ESatellite	sat;
	double		x;
	double		y;
	double		z;
	Satellite( ESatellite _sat, double _x, double _y, double _z ) : sat(_sat), x(_x), y(_y), z(_z){}
};

#endif /* INCLUDE_EPHEMDEFS_H_ */
