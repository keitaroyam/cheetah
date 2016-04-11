/*
 *  cheetahmodules.h
 *  cheetah
 *
 *  Created by Anton Barty on 7/2/11.
 *  Copyright 2011 CFEL. All rights reserved.
 *
 */

/*
 *	Function prototypes
 */

// peakfinders.cpp
int peakfinder3(tPeakList*, float*, char*, long, long, long, long, float, float, long, long, long);
int peakfinder6(tPeakList*, float*, char*, long, long, long, long, float, float, long, long, long, float);
int peakfinder8(tPeakList*, float*, char*, float*, long, long, long, long, float, float, long, long, long);
int peakfinder8old(tPeakList*, float*, char*, float*, long, long, long, long, float, float, long, long, long);
int killNearbyPeaks(tPeakList*, float );

