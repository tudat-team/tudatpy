/*
 *  EigenRoutines.cpp
 *  Propagator
 *
 *  Created by Melman, J.C.P. (jcpmelman) on 8/7/09.
 *  Copyright 2009 Delft University of Technology. All rights reserved.
 *
 */

#include "linearAlgebra.h"

// Determine the cosine of the angle between two vectors
double cosAngle(const Vector3d& v0, const Vector3d& v1)
{
	// Determine the length of the vectors
	double l0 = v0.norm( );
	double l1 = v1.norm( );

	// Check for a division by zero, which is obviously not allowed
	if (l0 > 0.0 && l1 > 0.0)
	{
		// Normalize both vectors
		Vector3d v0n = v0 / l0;
		Vector3d v1n = v1 / l1;

		// Get the cosine of the angle by dotting the normalized vectors
		double dotNorm = v0n.dot(v1n);
		// Explicitly define the extreme cases, which can give problems with the
		// acos function.
		if (dotNorm >= 1.0)
			return 1.0;
		else if (dotNorm <= -1.0)
			return -1.0;
		// Determine the actual angle
		else
			return dotNorm;
	}
	else
	{
		return 1.0e10;
	}
}

// Determine the angle between two vectors
double angle(const Vector3d& v0, const Vector3d& v1)
{
	// Determine the cosine of the angle by using another routine
	double dotNorm = cosAngle( v0, v1 );

	// Check for a division by zero, which is obviously not allowed
	if (fabs(dotNorm) <= 1.0)
	{
		return acos(dotNorm);

	}
	else
	{
		return 1.0e10;
	}
}

// Determine the average of the components of a vector
double average(const VectorXd& v)
{
	return v.sum( ) / v.rows( );
}

// Determine the standard deviation of the components of a vector
double standDev(const VectorXd& v)
{
	double variance = 0.;
	double ave = average( v );
	for( int i = 0; i < v.rows( ); i++ )
	{
		variance += pow( ( v(i) - ave ), 2. );
	}
	variance /= v.rows( ) - 1;
	return sqrt( variance );
}
