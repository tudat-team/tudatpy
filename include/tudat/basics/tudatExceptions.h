/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EXCEPTIONS_H
#define TUDAT_EXCEPTIONS_H

#include <iostream>
#include <string>

namespace tudat
{

namespace exceptions
{

class TudatError : public std::runtime_error
{
private:
public:
    TudatError( const std::string& errorMessage ): std::runtime_error( errorMessage ) { }
    ~TudatError( ) { }
};

class InterpolationOutOfBoundsError : public TudatError
{
private:
public:
    InterpolationOutOfBoundsError( const std::string& errorMessage ): exceptions::TudatError( errorMessage ) { }
    ~InterpolationOutOfBoundsError( ) { }
};

class LagrangeInterpolationOutOfBoundsError : public InterpolationOutOfBoundsError
{
private:
public:
    LagrangeInterpolationOutOfBoundsError( const std::string& errorMessage ): exceptions::InterpolationOutOfBoundsError( errorMessage ) { }
    ~LagrangeInterpolationOutOfBoundsError( ) { }
};

}  // namespace exceptions

}  // namespace tudat

#endif  // TUDAT_EXCEPTIONS_H
