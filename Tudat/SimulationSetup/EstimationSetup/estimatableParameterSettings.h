/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ESTIMATABLEPARAMETERSETTINGS_H
#define TUDAT_ESTIMATABLEPARAMETERSETTINGS_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Class for providing settings for parameters to estimate
/*!
 *  Class for providing settings for parameters to estimate. This class is a functional base class for parameters that
 *  require no information in addition to their type.
 *  This class can be used for the easy setup of parameter objects (see createEstimatableParameters.h), but users may also
 *  chose to do so manually. (Derived) Class members are all public, for ease of access and modification.
 */
class EstimatableParameterSettings
{
public:

    //! Constructor.
    /*!
     *  Constructor, takes parameter type and body of which it is a property.
     *  \param associatedBody Body of which parameter is a property.
     *  \param parameterType Type of parameter.
     *  \param pointOnBodyId Reference point on body associated with parameter (empty by default).
     */
    EstimatableParameterSettings( const std::string associatedBody ,
                                  const EstimatebleParametersEnum parameterType,
                                  const std::string pointOnBodyId = "" ):
        parameterType_( std::make_pair( parameterType, std::make_pair( associatedBody, pointOnBodyId ) ) ){ }

    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    virtual ~EstimatableParameterSettings( ){ }


    //! Identifier for parameter.
    /*!
     *  Identifier for parameter, contains type of parameter and body of which parameter is a property.
     */
    EstimatebleParameterIdentifier parameterType_;

};

//! Class for providing settings spherical harmonic coefficient(s) parameter
class SphericalHarmonicEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor
    /*!
     * Constructor, used to set user-defined list of degrees and orders.
     * \param blockIndices List of degrees and orders that are to estimated (first and second of each entry are
     * degree and order).
     * \param associatedBody Body for which coefficients are to be estimated.
     * \param parameterType Type of parameter that is to be estimated (must be spherical_harmonics_cosine_coefficient_block
     * of spherical_harmonics_sine_coefficient_block).
     */
    SphericalHarmonicEstimatableParameterSettings( const std::vector< std::pair< int, int > > blockIndices,
                                                   const std::string associatedBody,
                                                   const EstimatebleParametersEnum parameterType ):
        EstimatableParameterSettings( associatedBody, parameterType ), blockIndices_( blockIndices )
    {
        if( ( parameterType != spherical_harmonics_cosine_coefficient_block ) &&
                ( parameterType != spherical_harmonics_sine_coefficient_block ) )
        {
            throw std::runtime_error ("Error when making spherical harmonic parameter settings, input parameter type is inconsistent." );
        }
    }

    //! Constructor
    /*!
     * Constructor, used to set a full block of degrees and orders that are to be estimated
     * \param minimumDegree Minimum degree of field that is to be estimated.
     * \param minimumOrder Minimum order of field that is to be estimated.
     * \param maximumDegree Maximum degree of field that is to be estimated.
     * \param maximumOrder Maximum order of field that is to be estimated.
     * \param associatedBody Body for which coefficients are to be estimated.
     * \param parameterType Type of parameter that is to be estimated (must be spherical_harmonics_cosine_coefficient_block
     * of spherical_harmonics_sine_coefficient_block).
     */
    SphericalHarmonicEstimatableParameterSettings( const int minimumDegree,
                                                   const int minimumOrder,
                                                   const int maximumDegree,
                                                   const int maximumOrder,
                                                   const std::string associatedBody,
                                                   const EstimatebleParametersEnum parameterType ):
        EstimatableParameterSettings( associatedBody, parameterType )
    {
        if( ( parameterType != spherical_harmonics_cosine_coefficient_block ) &&
                ( parameterType != spherical_harmonics_sine_coefficient_block ) )
        {
            throw std::runtime_error( "Error when making spherical harmonic parameter settings, input parameter type is inconsistent." );
        }

        for( int i = minimumDegree; i <= maximumDegree; i++ )
        {
            for( int j = minimumOrder; ( ( j <= i ) && ( j <= maximumOrder ) ); j++ )
            {
                blockIndices_.push_back( std::make_pair( i, j ) );
            }
        }
    }

    //! List of degrees and orders that are to estimated (first and second of each entry are degree and order.
    std::vector< std::pair< int, int > > blockIndices_;
};

//! Class to define settings for estimating an initial translational state.
template< typename InitialStateParameterType >
class InitialTranslationalStateEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor, sets initial value of translational state.
    /*!
     * Constructor, sets initial value of translational state.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialStateValue Current value of initial state (w.r.t. centralBody)
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, 6, 1 > initialStateValue,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, initial_body_state ), initialTime_( TUDAT_NAN ),
        initialStateValue_( initialStateValue ),
        centralBody_( centralBody ), frameOrientation_( frameOrientation ){ }

    //! Constructor, without initial value of translational state.
    /*!
     * Constructor, without initial value of translational state. Current initial state is retrieved from environment
     * (ephemeris objects) during creation of parameter object.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialTime Time at which initial state is defined.
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const double initialTime,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, initial_body_state ), initialTime_( initialTime ),
        centralBody_( centralBody ), frameOrientation_( frameOrientation ){ }

    //! Time at which initial state is defined (NaN for user-defined initial state value).
    double initialTime_;

    //! Current value of initial state (w.r.t. centralBody), set manually by used.
    Eigen::Matrix< InitialStateParameterType, 6, 1 > initialStateValue_;

    //! Body w.r.t. which the initial state is to be estimated.
    std::string centralBody_;

    //! Orientation of the frame in which the state is defined.
    std::string frameOrientation_;

};

template< typename InitialStateParameterType >
class ArcWiseInitialTranslationalStateEstimatableParameterSettings: public EstimatableParameterSettings
{
public:
    ArcWiseInitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStateValue,
            const std::vector< std::pair< double, double > >& arcStartAndEndTimes,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, arc_wise_initial_body_state ), initialStateValue_( initialStateValue ),
        arcStartAndEndTimes_( arcStartAndEndTimes ), centralBody_( centralBody ), frameOrientation_( frameOrientation ),
        isStateSet_( 1 ){ }

    ArcWiseInitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const std::vector< std::pair< double, double > >& arcStartAndEndTimes,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, arc_wise_initial_body_state ),
        arcStartAndEndTimes_( arcStartAndEndTimes ), centralBody_( centralBody ), frameOrientation_( frameOrientation ),
        isStateSet_( 0 ){ }


    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStateValue_;

    std::vector< std::pair< double, double > > arcStartAndEndTimes_;

    std::string centralBody_;

    std::string frameOrientation_;

    bool isStateSet_;

};

} // namespace estimatable_parameters

} // namespace tudat
#endif // TUDAT_ESTIMATABLEPARAMETERSETTINGS_H
