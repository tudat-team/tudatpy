/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/External/JsonInterface/Mathematics/interpolation.h"

#include "aerodynamics.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicSettings )
{
    if ( ! aerodynamicSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Body::Aerodynamics;

    const AerodynamicCoefficientTypes aerodynamicCoefficientType =
            aerodynamicSettings->getAerodynamicCoefficientType( );
    jsonObject[ K::coefficientsType ] = aerodynamicCoefficientType;
    jsonObject[ K::referenceArea ] = aerodynamicSettings->getReferenceArea( );

    switch ( aerodynamicCoefficientType ) {
    case constant_aerodynamic_coefficients:
    {
        boost::shared_ptr< ConstantAerodynamicCoefficientSettings > constantAerodynamicSettings =
                boost::dynamic_pointer_cast< ConstantAerodynamicCoefficientSettings >( aerodynamicSettings );
        enforceNonNullPointer( constantAerodynamicSettings );
        jsonObject[ K::referenceArea ] = constantAerodynamicSettings->getReferenceArea( );
        jsonObject[ K::forceCoefficients ] = constantAerodynamicSettings->getConstantForceCoefficient( );

        // Moments
        if ( ! isNaN( constantAerodynamicSettings->getReferenceLength( ) ) )
        {
            jsonObject[ K::referenceLength ] = constantAerodynamicSettings->getReferenceLength( );
            jsonObject[ K::lateralReferenceLength ] = constantAerodynamicSettings->getLateralReferenceLength( );
            jsonObject[ K::momentReferencePoint ] = constantAerodynamicSettings->getMomentReferencePoint( );
            jsonObject[ K::momentCoefficients ] = constantAerodynamicSettings->getConstantMomentCoefficient( );
        }

        jsonObject[ K::areCoefficientsInAerodynamicFrame ] =
                constantAerodynamicSettings->getAreCoefficientsInAerodynamicFrame( );
        jsonObject[ K::areCoefficientsInNegativeAxisDirection ] =
                constantAerodynamicSettings->getAreCoefficientsInNegativeAxisDirection( );
        return;
    }
    case tabulated_coefficients:
    {
        boost::shared_ptr< TabulatedAerodynamicCoefficientSettings< 1 > > tabulatedAerodynamicSettings =
                boost::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< 1 > >( aerodynamicSettings );
        if ( tabulatedAerodynamicSettings )  // uni-dimensional
        {
            jsonObject[ K::independentVariables ] =
                    getMapKeys< std::map >( tabulatedAerodynamicSettings->getForceCoefficients( ) );
            jsonObject[ K::forceCoefficients ] =
                    getMapValues< std::map >( tabulatedAerodynamicSettings->getForceCoefficients( ) );
            jsonObject[ K::referenceArea ] = tabulatedAerodynamicSettings->getReferenceArea( );

            // Moments
            if ( ! isNaN( tabulatedAerodynamicSettings->getReferenceLength( ) ) )
            {
                jsonObject[ K::momentCoefficients ] =
                        getMapValues< std::map >( tabulatedAerodynamicSettings->getMomentCoefficients( ) );
                jsonObject[ K::referenceLength ] = tabulatedAerodynamicSettings->getReferenceLength( );
                jsonObject[ K::lateralReferenceLength ] = tabulatedAerodynamicSettings->getLateralReferenceLength( );
                jsonObject[ K::momentReferencePoint ] = tabulatedAerodynamicSettings->getMomentReferencePoint( );
            }

            jsonObject[ K::independentVariableName ] =
                    tabulatedAerodynamicSettings->getIndependentVariableNames( ).front( );
            jsonObject[ K::interpolator ] = tabulatedAerodynamicSettings->getInterpolationSettings( );

            jsonObject[ K::areCoefficientsInAerodynamicFrame ] =
                    tabulatedAerodynamicSettings->getAreCoefficientsInAerodynamicFrame( );
            jsonObject[ K::areCoefficientsInNegativeAxisDirection ] =
                    tabulatedAerodynamicSettings->getAreCoefficientsInNegativeAxisDirection( );
        }
        else
        {
            std::cerr << "Multi-dimensional tabulated aerodynamic coefficients not yet supported by json_interface."
                      << std::endl;
        }
        return;
    }
    default:
        handleUnimplementedEnumValue( aerodynamicCoefficientType, aerodynamicCoefficientTypes,
                                      unsupportedAerodynamicCoefficientTypes );
    }
}

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
void from_json( const json& jsonObject, boost::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicSettings )
{
    using namespace aerodynamics;
    using namespace interpolators;
    using namespace json_interface;
    using K = Keys::Body::Aerodynamics;

    // Aerodynamic coefficient type (constant by default)
    const AerodynamicCoefficientTypes coefficientsType =
            getValue( jsonObject, K::coefficientsType, constant_aerodynamic_coefficients );

    // Reference area (use fallback area if reference area not provided, final value cannont be NaN)
    double fallbackReferenceArea = getNumeric< double >(
                jsonObject, SpecialKeys::up / K::referenceArea, TUDAT_NAN, true );

    switch ( coefficientsType )
    {
    case constant_aerodynamic_coefficients:
    {
        // Read forceCoefficients. If not defined, use [ dragCoefficient, 0, 0 ].
        Eigen::Vector3d forceCoefficients = Eigen::Vector3d::Zero( );
        if ( defined( jsonObject, K::forceCoefficients ) )
        {
            forceCoefficients = getValue< Eigen::Vector3d >( jsonObject, K::forceCoefficients );
        }
        else
        {
            forceCoefficients( 0 ) = getValue< double >( jsonObject, K::dragCoefficient );
        }

        if ( defined( jsonObject, K::momentCoefficients ) )  // moments
        {
            ConstantAerodynamicCoefficientSettings defaults( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN,
                                                             Eigen::Vector3d( ), Eigen::Vector3d( ) );
            aerodynamicSettings = boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                        getNumeric< double >( jsonObject, K::referenceLength ),
                        getNumeric( jsonObject, K::referenceArea, fallbackReferenceArea ),
                        getNumeric< double >( jsonObject, K::lateralReferenceLength ),
                        getValue< Eigen::Vector3d >( jsonObject, K::momentReferencePoint ),
                        forceCoefficients,
                        getValue< Eigen::Vector3d >( jsonObject, K::momentCoefficients ),
                        getValue( jsonObject, K::areCoefficientsInAerodynamicFrame,
                                  defaults.getAreCoefficientsInAerodynamicFrame( ) ),
                        getValue( jsonObject, K::areCoefficientsInNegativeAxisDirection,
                                  defaults.getAreCoefficientsInNegativeAxisDirection( ) ) );
        }
        else  // no moments
        {
            ConstantAerodynamicCoefficientSettings defaults( TUDAT_NAN, Eigen::Vector3d( ) );
            aerodynamicSettings = boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                        getNumeric( jsonObject, K::referenceArea, fallbackReferenceArea ),
                        forceCoefficients,
                        getValue( jsonObject, K::areCoefficientsInAerodynamicFrame,
                                  defaults.getAreCoefficientsInAerodynamicFrame( ) ),
                        getValue( jsonObject, K::areCoefficientsInNegativeAxisDirection,
                                  defaults.getAreCoefficientsInNegativeAxisDirection( ) ) );
        }
        return;
    }
    case tabulated_coefficients:
    {
        const unsigned int dimensions = getValue< unsigned int >( jsonObject, K::numberOfDimensions, 1 );
        if ( dimensions == 1 )
        {
            TabulatedAerodynamicCoefficientSettings< 1 > defaults(
            { }, { }, { }, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, Eigen::Vector3d( ), mach_number_dependent, NULL );
            aerodynamicSettings = boost::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >(
                        getValue< std::vector< double > >( jsonObject, K::independentVariables ),
                        getValue< std::vector< Eigen::Vector3d > >( jsonObject, K::forceCoefficients ),
                        getValue< std::vector< Eigen::Vector3d > >( jsonObject, K::momentCoefficients ),
                        getNumeric< double >( jsonObject, K::referenceLength ),
                        getNumeric( jsonObject, K::referenceArea, fallbackReferenceArea ),
                        getNumeric< double >( jsonObject, K::lateralReferenceLength ),
                        getValue< Eigen::Vector3d >( jsonObject, K::momentReferencePoint ),
                        getValue< AerodynamicCoefficientsIndependentVariables >(
                            jsonObject, K::independentVariableName ),
                        getValue< boost::shared_ptr< InterpolatorSettings > >( jsonObject, K::interpolator ),
                        getValue( jsonObject, K::areCoefficientsInAerodynamicFrame,
                                  defaults.getAreCoefficientsInAerodynamicFrame( ) ),
                        getValue( jsonObject, K::areCoefficientsInNegativeAxisDirection,
                                  defaults.getAreCoefficientsInNegativeAxisDirection( ) ) );
        }
        else
        {
            /*
            if ( defined( jsonObject, K::momentCoefficients ) )  // moments
            {
                TabulatedAerodynamicCoefficientSettings< dimensions > defaults(
                { }, { }, { }, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, Eigen::Vector3d( ), { } );
            }
            else  // no moments
            {
                TabulatedAerodynamicCoefficientSettings< dimensions > defaults( { }, { }, TUDAT_NAN, { } );
            }
            */
            std::cerr << "Multi-dimensional tabulated aerodynamic coefficients not yet supported by json_interface."
                      << std::endl;
            throw;
        }
        return;
    }
    default:
        handleUnimplementedEnumValue( coefficientsType, aerodynamicCoefficientTypes,
                                      unsupportedAerodynamicCoefficientTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
