/*
 * SatellitesEphemeris.h
 *
 *  Created on: March 2019
 *      Author: Georgi Dragandjikov
 *
 *  Interface to the libSatellitesEphemeris.so
 *
 *
 *  libSatellitesEphemeris.so represents some of the ephemeris calculations I implemented in my
 *  planetarium software Power Age Sky Simulator (released between 1998 and 2003, www.powerageskysimulator.com).
 *  I mainly used the incredible "Astronomical Algorithms" by Jean Meeus, Second Edition.
 */

#ifndef INCLUDE_SATELLITESEPHEMERIS_H_
#define INCLUDE_SATELLITESEPHEMERIS_H_

#include <vector>
#include "EphemDefs.h"
#include "EphemTime.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates visible positions of 7 Saturn satellites: Mimas, Enceladus, Tethys, Dione, Rhea, Titan and Hyperion.
// The unit of measurement is the equatorial radius of Saturn.
// X axis coinciding with the equator of the planet and is positive to west and negative to east.
// Y axis coinciding with the planet's rotation axis and is positive to north and negative to south.
// Z is negative if the satellite is between Saturn and the Earth, and positive if it is behind the planet.
//
// parameters:
//    - ephemTime; the given calculation time;
//    - satellites; the calculated results.
//
void CalculateSaturnSatellites( EphemTime ephemTime, std::vector<Satellite>& satellites );
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the visible Saturn's ring position. Namely, the following parameters are calculated over the
// given ephemTime time:
//   - earthLat; the latitude of the Earth as it could be seen from Saturn. When it is positive, we see
//     Saturn's north pole, and when it is negative - the Saturn's south pole. Measured in degrees;
//   - sunLat; the latitude of the Sun as it could be seen from Saturn. When it is positive, the Saturn's
//     north pole is illuminated, and when it is negative - the Saturn's south pole is illuminated; Measured in degrees;
//   - posAngle; this is the angle between the Saturn's north pole rotation axis and the northern celestial pole.
//     Measured in degrees;
//   - majorAxe, minorAxe; the major and the minor axes of the outer edge of the outer ring. Measured in arcseconds.
//
void CalculateSaturnRingPosition( EphemTime ephemTime,
								  double& earthLat, double& sunLat,
								  double& posAngle,
								  double& majorAxe, double& minorAxe );
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the positions of the four great Jupiter moons: Io, Europa, Ganymede and Callisto.
// The unit of measurement is the equatorial radius of Jupiter. Low accuracy method is used.
// X axis coinciding with the equator of the planet and is positive to west and negative to east.
// Y axis coinciding with the planet's rotation axis and is positive to north and negative to south.
// Z is not used (in this particular low accuracy method of calculation). Still, each satellite's
// position behind/before Jupiter could be easily considered by making two consecutive calculations and
// finding the move direction by X.
//
// parameters:
//    - ephemTime; the given calculation time;
//    - satellites; the calculated results.
//
void CalculateJupiterSatellites( EphemTime ephemTime, std::vector<Satellite>& satellites );
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif /* INCLUDE_SATELLITESEPHEMERIS_H_ */
