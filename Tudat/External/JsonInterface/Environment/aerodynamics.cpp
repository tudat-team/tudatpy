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
        boost::shared_ptr< FromFileAerodynamicCoefficientSettings > fromFileAerodynamicCoefficientSettings =
                boost::dynamic_pointer_cast< FromFileAerodynamicCoefficientSettings >( aerodynamicSettings );
        enforceNonNullPointer( fromFileAerodynamicCoefficientSettings );
        const bool hasMoments = ! isNaN( fromFileAerodynamicCoefficientSettings->getReferenceLength( ) );

        boost::shared_ptr< TabulatedAerodynamicCoefficientSettings< 1 > > tabulated1AerodynamicSettings =
                boost::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< 1 > >( aerodynamicSettings );
        if ( tabulated1AerodynamicSettings )  // uni-dimensional
        {
            jsonObject[ K::independentVariables ] =
                    getMapKeys< std::map >( tabulated1AerodynamicSettings->getForceCoefficients( ) );
            jsonObject[ K::forceCoefficients ] =
                    getMapValues< std::map >( tabulated1AerodynamicSettings->getForceCoefficients( ) );

            if ( hasMoments )
            {
                jsonObject[ K::momentCoefficients ] =
                        getMapValues< std::map >( tabulated1AerodynamicSettings->getMomentCoefficients( ) );
            }

            jsonObject[ K::interpolator ] = tabulated1AerodynamicSettings->getInterpolationSettings( );
        }
        else  // N-dimensional
        {
            const std::map< int, std::string > forceCoefficientsFilesMap =
                    fromFileAerodynamicCoefficientSettings->getForceCoefficientsFiles( );
            jsonObject[ K::forceCoefficients ] = { "", "", "" };
            for ( auto entry : forceCoefficientsFilesMap )
            {
                jsonObject[ K::forceCoefficients ][ entry.first ] = path( entry.second );
            }

            if ( hasMoments )
            {
                const std::map< int, std::string > momentCoefficientsFilesMap =
                        fromFileAerodynamicCoefficientSettings->getMomentCoefficientsFiles( );
                jsonObject[ K::momentCoefficients ] = { "", "", "" };
                for ( auto entry : momentCoefficientsFilesMap )
                {
                    jsonObject[ K::momentCoefficients ][ entry.first ] = path( entry.second );
                }
            }
        }

        if ( hasMoments )
        {
            jsonObject[ K::referenceLength ] = fromFileAerodynamicCoefficientSettings->getReferenceLength( );
            jsonObject[ K::lateralReferenceLength ] =
                    fromFileAerodynamicCoefficientSettings->getLateralReferenceLength( );
            jsonObject[ K::momentReferencePoint ] = fromFileAerodynamicCoefficientSettings->getMomentReferencePoint( );
        }

        jsonObject[ K::independentVariableNames ] =
                fromFileAerodynamicCoefficientSettings->getIndependentVariableNames( );
        jsonObject[ K::areCoefficientsInAerodynamicFrame ] =
                fromFileAerodynamicCoefficientSettings->getAreCoefficientsInAerodynamicFrame( );
        jsonObject[ K::areCoefficientsInNegativeAxisDirection ] =
                fromFileAerodynamicCoefficientSettings->getAreCoefficientsInNegativeAxisDirection( );

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
        const json jsonForceCoefficients = getValue< json >( jsonObject, K::forceCoefficients );
        const json jsonMomentCoefficients = getValue( jsonObject, K::momentCoefficients, json( ) );
        const bool useMoments = ! jsonMomentCoefficients.is_null( );

        bool readFromFiles = true;
        std::map< int, std::string > forceCoefficientsFiles;
        std::vector< Eigen::Vector3d > forceCoefficients;
        std::map< int, std::string > momentCoefficientsFiles;
        std::vector< Eigen::Vector3d > momentCoefficients;
        try
        {
            const std::vector< path > forceCoefficientsFilesVector =
                    getAs< std::vector< path > >( jsonForceCoefficients );
            for ( unsigned int i = 0; i < forceCoefficientsFilesVector.size( ); ++i )
            {
                if ( ! jsonObject.at( K::forceCoefficients ).at( i ).get< std::string >( ).empty( ) )
                {
                    forceCoefficientsFiles[ i ] = forceCoefficientsFilesVector.at( i ).string( );
                }
            }

            if ( useMoments )
            {
                const std::vector< path > momentCoefficientsFilesVector =
                        getAs< std::vector< path > >( jsonMomentCoefficients );
                for ( unsigned int i = 0; i < momentCoefficientsFilesVector.size( ); ++i )
                {
                    if ( ! jsonObject.at( K::momentCoefficients ).at( i ).get< std::string >( ).empty( ) )
                    {
                        momentCoefficientsFiles[ i ] = momentCoefficientsFilesVector.at( i ).string( );
                    }
                }
            }
        }
        catch ( ... )
        {
            readFromFiles = false;
            forceCoefficients = getAs< std::vector< Eigen::Vector3d > >( jsonForceCoefficients );
            if ( useMoments )
            {
                momentCoefficients = getAs< std::vector< Eigen::Vector3d > >( jsonMomentCoefficients );
            }
        }

        const double referenceArea = getNumeric( jsonObject, K::referenceArea, fallbackReferenceArea );
        const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames =
                getValue< std::vector< AerodynamicCoefficientsIndependentVariables > >(
                    jsonObject, K::independentVariableNames );
        const bool areCoefficientsInAerodynamicFrame =
                getValue( jsonObject, K::areCoefficientsInAerodynamicFrame, true );
        const bool areCoefficientsInNegativeAxisDirection =
                getValue( jsonObject, K::areCoefficientsInNegativeAxisDirection, true );

        if ( readFromFiles )
        {
            if ( useMoments )
            {
                aerodynamicSettings = readTabulatedAerodynamicCoefficientsFromFiles(
                            forceCoefficientsFiles,
                            momentCoefficientsFiles,
                            getNumeric< double >( jsonObject, K::referenceLength ),
                            referenceArea,
                            getNumeric< double >( jsonObject, K::lateralReferenceLength ),
                            getValue< Eigen::Vector3d >( jsonObject, K::momentReferencePoint ),
                            independentVariableNames,
                            areCoefficientsInAerodynamicFrame,
                            areCoefficientsInNegativeAxisDirection );
            }
            else
            {
                aerodynamicSettings = readTabulatedAerodynamicCoefficientsFromFiles(
                            forceCoefficientsFiles,
                            referenceArea,
                            independentVariableNames,
                            areCoefficientsInAerodynamicFrame,
                            areCoefficientsInNegativeAxisDirection );
            }
        }
        else
        {
            if ( useMoments )
            {
                aerodynamicSettings = boost::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >(
                            getValue< std::vector< double > >( jsonObject, K::independentVariables ),
                            forceCoefficients,
                            momentCoefficients,
                            getNumeric< double >( jsonObject, K::referenceLength ),
                            referenceArea,
                            getNumeric< double >( jsonObject, K::lateralReferenceLength ),
                            getValue< Eigen::Vector3d >( jsonObject, K::momentReferencePoint ),
                            independentVariableNames.front( ),
                            getValue< boost::shared_ptr< InterpolatorSettings > >( jsonObject, K::interpolator ),
                            areCoefficientsInAerodynamicFrame,
                            areCoefficientsInNegativeAxisDirection );
            }
            else
            {
                aerodynamicSettings = boost::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >(
                            getValue< std::vector< double > >( jsonObject, K::independentVariables ),
                            forceCoefficients,
                            referenceArea,
                            independentVariableNames.front( ),
                            getValue< boost::shared_ptr< InterpolatorSettings > >( jsonObject, K::interpolator ),
                            areCoefficientsInAerodynamicFrame,
                            areCoefficientsInNegativeAxisDirection );
            }
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
