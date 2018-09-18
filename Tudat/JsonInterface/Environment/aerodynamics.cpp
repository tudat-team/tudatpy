/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/JsonInterface/Mathematics/interpolation.h"

#include "Tudat/JsonInterface/Environment/aerodynamics.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicSettings )
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
        std::shared_ptr< ConstantAerodynamicCoefficientSettings > constantAerodynamicSettings =
                std::dynamic_pointer_cast< ConstantAerodynamicCoefficientSettings >( aerodynamicSettings );
        assertNonnullptrPointer( constantAerodynamicSettings );
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
        std::shared_ptr< TabulatedAerodynamicCoefficientSettingsBase > fromFileAerodynamicCoefficientSettings =
                std::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettingsBase >( aerodynamicSettings );
        assertNonnullptrPointer( fromFileAerodynamicCoefficientSettings );
        const bool hasMoments = ! isNaN( fromFileAerodynamicCoefficientSettings->getReferenceLength( ) );

        std::shared_ptr< TabulatedAerodynamicCoefficientSettings< 1 > > tabulated1AerodynamicSettings =
                std::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< 1 > >( aerodynamicSettings );
        if ( tabulated1AerodynamicSettings )  // uni-dimensional
        {
            jsonObject[ K::independentVariableValues ] =
                    getMapKeys< std::map >( tabulated1AerodynamicSettings->getForceCoefficients( ) );
            jsonObject[ K::forceCoefficients ] =
                    getMapValues< std::map >( tabulated1AerodynamicSettings->getForceCoefficients( ) );

            if ( hasMoments )
            {
                jsonObject[ K::momentCoefficients ] =
                        getMapValues< std::map >( tabulated1AerodynamicSettings->getMomentCoefficients( ) );
            }

            jsonObject[ K::interpolator ] = tabulated1AerodynamicSettings->getInterpolatorSettings( );
        }
        else  // N-dimensional
        {
            const std::map< int, std::string > forceCoefficientsFilesMap =
                    fromFileAerodynamicCoefficientSettings->getForceCoefficientsFiles( );
            jsonObject[ K::forceCoefficients ] = { "", "", "" };
            for ( auto entry : forceCoefficientsFilesMap )
            {
                jsonObject[ K::forceCoefficients ][ entry.first ] = boost::filesystem::path( entry.second );
            }

            if ( hasMoments )
            {
                const std::map< int, std::string > momentCoefficientsFilesMap =
                        fromFileAerodynamicCoefficientSettings->getMomentCoefficientsFiles( );
                jsonObject[ K::momentCoefficients ] = { "", "", "" };
                for ( auto entry : momentCoefficientsFilesMap )
                {
                    jsonObject[ K::momentCoefficients ][ entry.first ] = boost::filesystem::path( entry.second );
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
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicSettings )
{
    using namespace aerodynamics;
    using namespace interpolators;
    using namespace json_interface;
    using K = Keys::Body::Aerodynamics;

    // Aerodynamic coefficient type (constant by default)
    const AerodynamicCoefficientTypes coefficientsType =
            getValue( jsonObject, K::coefficientsType, constant_aerodynamic_coefficients );

    // Reference area (either from the current object or from the current object's parent, i.e. the body)
    const double referenceArea =
            getValue< double >( jsonObject, { K::referenceArea, SpecialKeys::up / Keys::Body::referenceArea } );

    switch ( coefficientsType )
    {
    case constant_aerodynamic_coefficients:
    {
        // Read forceCoefficients. If not defined, use [ dragCoefficient, 0, 0 ].
        Eigen::Vector3d forceCoefficients = Eigen::Vector3d::Zero( );
        if ( isDefined( jsonObject, K::forceCoefficients ) )
        {
            forceCoefficients = getValue< Eigen::Vector3d >( jsonObject, K::forceCoefficients );
        }
        else
        {
            forceCoefficients( 0 ) = getValue< double >( jsonObject, K::dragCoefficient );
        }

        if ( isDefined( jsonObject, K::momentCoefficients ) )  // moments
        {
            ConstantAerodynamicCoefficientSettings defaults( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN,
                                                             Eigen::Vector3d( ), Eigen::Vector3d( ) );
            aerodynamicSettings = std::make_shared< ConstantAerodynamicCoefficientSettings >(
                        getValue< double >( jsonObject, K::referenceLength ),
                        referenceArea,
                        getValue< double >( jsonObject, K::lateralReferenceLength ),
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
            aerodynamicSettings = std::make_shared< ConstantAerodynamicCoefficientSettings >(
                        referenceArea,
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
        const nlohmann::json jsonForceCoefficients = getValue< nlohmann::json >( jsonObject, K::forceCoefficients );
        const nlohmann::json jsonMomentCoefficients = getValue( jsonObject, K::momentCoefficients, nlohmann::json( ) );
        const bool useMoments = ! jsonMomentCoefficients.is_null( );

        bool readFromFiles = true;
        std::map< int, std::string > forceCoefficientsFiles;
        std::vector< Eigen::Vector3d > forceCoefficients;
        std::map< int, std::string > momentCoefficientsFiles;
        std::vector< Eigen::Vector3d > momentCoefficients;
        try
        {
            const std::vector< boost::filesystem::path > forceCoefficientsFilesVector =
                    getAs< std::vector< boost::filesystem::path > >( jsonForceCoefficients );
            for ( unsigned int i = 0; i < forceCoefficientsFilesVector.size( ); ++i )
            {
                if ( ! jsonObject.at( K::forceCoefficients ).at( i ).get< std::string >( ).empty( ) )
                {
                    forceCoefficientsFiles[ i ] = forceCoefficientsFilesVector.at( i ).string( );
                }
            }

            if ( useMoments )
            {
                const std::vector< boost::filesystem::path > momentCoefficientsFilesVector =
                        getAs< std::vector< boost::filesystem::path > >( jsonMomentCoefficients );
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
                            getValue< double >( jsonObject, K::referenceLength ),
                            referenceArea,
                            getValue< double >( jsonObject, K::lateralReferenceLength ),
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
                aerodynamicSettings = std::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >(
                            getValue< std::vector< double > >( jsonObject, K::independentVariableValues ),
                            forceCoefficients,
                            momentCoefficients,
                            getValue< double >( jsonObject, K::referenceLength ),
                            referenceArea,
                            getValue< double >( jsonObject, K::lateralReferenceLength ),
                            getValue< Eigen::Vector3d >( jsonObject, K::momentReferencePoint ),
                            independentVariableNames.front( ),
                            areCoefficientsInAerodynamicFrame,
                            areCoefficientsInNegativeAxisDirection,
                            getValue< std::shared_ptr< InterpolatorSettings > >( jsonObject, K::interpolator ) );
            }
            else
            {
                aerodynamicSettings = std::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >(
                            getValue< std::vector< double > >( jsonObject, K::independentVariableValues ),
                            forceCoefficients,
                            referenceArea,
                            independentVariableNames.front( ),
                            areCoefficientsInAerodynamicFrame,
                            areCoefficientsInNegativeAxisDirection,
                            getValue< std::shared_ptr< InterpolatorSettings > >( jsonObject, K::interpolator ) );
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
