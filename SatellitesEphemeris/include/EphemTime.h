/*
 * EphemTime.h
 *
 *  Created on: March 2019
 *      Author: Georgi Dragandjikov
 *
 *  EphemTime class;
 *  Between years 1620 and 1998 a precise dynamic time correction is available. For the the rest of the
 *  periods approximation formulas are used.
 *
 *
 *  libSatellitesEphemeris.so represents some of the ephemeris calculations I implemented in my
 *  planetarium software Power Age Sky Simulator (released between 1998 and 2003, www.powerageskysimulator.com).
 *  I mainly used the incredible "Astronomical Algorithms" by Jean Meeus, Second Edition.
 */

#ifndef INCLUDE_EPHEMTIME_H_
#define INCLUDE_EPHEMTIME_H_


constexpr int DT_ENTRIES = 190;
constexpr int DT_correction[DT_ENTRIES][2] = {

//			year	correction (seconds)
		{	1620,	121	},
		{	1622,	112	},
		{	1624,	103	},
		{	1626,	95	},
		{	1628,	88	},

		{	1630,	82	},
		{	1632,	77	},
		{	1634,	72	},
		{	1636,	68	},
		{	1638,	63	},

		{	1640,	60	},
		{	1642,	56	},
		{	1644,	53	},
		{	1646,	51	},
		{	1648,	48	},

		{	1650,	46	},
		{	1652,	44	},
		{	1654,	42	},
		{	1656,	40	},
		{	1658,	38	},

		{	1660,	35	},
		{	1662,	33	},
		{	1664,	31	},
		{	1666,	29	},
		{	1668,	26	},

		{	1670,	24	},
		{	1672,	22	},
		{	1674,	20	},
		{	1676,	18	},
		{	1678,	16	},

		{	1680,	14	},
		{	1682,	12	},
		{	1684,	11	},
		{	1686,	10	},
		{	1688,	9	},

		{	1690,	8	},
		{	1692,	7	},
		{	1694,	7	},
		{	1696,	7	},
		{	1698,	7	},

		{	1700,	7	},
		{	1702,	7	},
		{	1704,	8	},
		{	1706,	8	},
		{	1708,	9	},

		{	1710,	9	},
		{	1712,	9	},
		{	1714,	9	},
		{	1716,	9	},
		{	1718,	10	},

		{	1720,	10	},
		{	1722,	10	},
		{	1724,	10	},
		{	1726,	10	},
		{	1728,	10	},

		{	1730,	10	},
		{	1732,	10	},
		{	1734,	11	},
		{	1736,	11	},
		{	1738,	11	},

		{	1740,	11	},
		{	1742,	11	},
		{	1744,	12	},
		{	1746,	12	},
		{	1748,	12	},

		{	1750,	12	},
		{	1752,	13	},
		{	1754,	13	},
		{	1756,	13	},
		{	1758,	14	},

		{	1760,	14	},
		{	1762,	14	},
		{	1764,	14	},
		{	1766,	15	},
		{	1768,	15	},

		{	1770,	15	},
		{	1772,	15	},
		{	1774,	15	},
		{	1776,	16	},
		{	1778,	16	},

		{	1780,	16	},
		{	1782,	16	},
		{	1784,	16	},
		{	1786,	16	},
		{	1788,	16	},

		{	1790,	16	},
		{	1792,	15	},
		{	1794,	15	},
		{	1796,	14	},
		{	1798,	13	},

		{	1800,	13	},
		{	1802,	13	},
		{	1804,	12	},
		{	1806,	12	},
		{	1808,	12	},

		{	1810,	12	},
		{	1812,	12	},
		{	1814,	12	},
		{	1816,	12	},
		{	1818,	12	},

		{	1820,	12	},
		{	1822,	11	},
		{	1824,	10	},
		{	1826,	9	},
		{	1828,	8	},

		{	1830,	7	},
		{	1832,	6	},
		{	1834,	6	},
		{	1836,	5	},
		{	1838,	5	},

		{	1840,	5	},
		{	1842,	6	},
		{	1844,	6	},
		{	1846,	6	},
		{	1848,	7	},

		{	1850,	7	},
		{	1852,	7	},
		{	1854,	7	},
		{	1856,	8	},
		{	1858,	8	},

		{	1860,	8	},
		{	1862,	7	},
		{	1864,	6	},
		{	1866,	5	},
		{	1868,	3	},

		{	1870,	1	},
		{	1872,	-1	},
		{	1874,	-3	},
		{	1876,	-4	},
		{	1878,	-5	},

		{	1880,	-6	},
		{	1882,	-5	},
		{	1884,	-6	},
		{	1886,	-6	},
		{	1888,	-6	},

		{	1890,	-6	},
		{	1892,	-6	},
		{	1894,	-7	},
		{	1896,	-6	},
		{	1898,	-5	},

		{	1900,	-3	},
		{	1902,	0	},
		{	1904,	3	},
		{	1906,	5	},
		{	1908,	8	},

		{	1910,	10	},
		{	1912,	13	},
		{	1914,	16	},
		{	1916,	18	},
		{	1918,	20	},

		{	1920,	21	},
		{	1922,	22	},
		{	1924,	24	},
		{	1926,	24	},
		{	1928,	24	},

		{	1930,	24	},
		{	1932,	24	},
		{	1934,	24	},
		{	1936,	24	},
		{	1938,	24	},

		{	1940,	24	},
		{	1942,	25	},
		{	1944,	26	},
		{	1946,	27	},
		{	1948,	28	},

		{	1950,	29	},
		{	1952,	30	},
		{	1954,	31	},
		{	1956,	31	},
		{	1958,	32	},

		{	1960,	33	},
		{	1962,	34	},
		{	1964,	35	},
		{	1966,	37	},
		{	1968,	38	},

		{	1970,	40	},
		{	1972,	42	},
		{	1974,	45	},
		{	1976,	47	},
		{	1978,	49	},

		{	1980,	51	},
		{	1982,	52	},
		{	1984,	54	},
		{	1986,	55	},
		{	1988,	56	},

		{	1990,	57	},
		{	1992,	58	},
		{	1994,	60	},
		{	1996,	62	},
		{	1998,	63	}

};


class EphemTime{

	int		year;
	int		month;
	int		day;
	int		hour;
	int		minute;
	int		second;
	double	jd;			// universal time (julian day)
	double	jde;		// dynamical time (julian day), julian ephemeris day

public:
	EphemTime( int y, int mon, int d, int h, int min, int s );
	EphemTime( double jd_ut );
	int Year() { return year; }
	int Month() { return month; }
	int Day() { return day; }
	int Hour() { return hour; }
	int Minute() { return minute; }
	int Second() { return second; }
	double JD() { return jd; }
	double JDE() { return jde; }

private:
	void CalculateDate();
	void CalculateJD();
	void CalculateJDE();

};


#endif /* INCLUDE_EPHEMTIME_H_ */
