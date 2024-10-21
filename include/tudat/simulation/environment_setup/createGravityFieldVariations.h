/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEGRAVITYFIELDVARIATIONS_H
#define TUDAT_CREATEGRAVITYFIELDVARIATIONS_H

#include <string>
#include <vector>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/math/interpolators/createInterpolator.h"

namespace tudat
{

namespace simulation_setup
{

//! Object to define settings to be used for interpolating time-variations computed by an
//! environment model.
/*!
 * Object to define settings to be used for interpolating time-variations computed by an environment
 * model. For instance, gravity field variations may be computed a priori (in certain cases) and
 * interpolated during propagation, which can significantly reduce the computation time.class
 * ModelInterpolationSettings
 */
class ModelInterpolationSettings
{
public:

    //! Constructor.
    /*! Constructor
     * \param initialTime Start time for interpolator.
     * \param finalTime End time for interpolator.
     * \param timeStep Time step with which to evaluate model, and provide input to interpolator
     * \param interpolatorSettings Settings to use to crate the interpolator (i.e. type and any
     * required associated information).
     */
    ModelInterpolationSettings(
            const double initialTime = 0.0,
            const double finalTime = 0.0,
            const double timeStep = 0.0,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ) ):
        interpolatorSettings_( interpolatorSettings ), initialTime_( initialTime ),
        finalTime_( finalTime ), timeStep_( timeStep ){ }

    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings_;

    //! Start time for interpolator.
    double initialTime_;

    //! End time for interpolator.
    double finalTime_;

    //! Time step with which to evaluate model, and provide input to interpolator
    double timeStep_;
};

//! Base class for defining settings for gravity field variations.
/*!
 * Base class (non-functiona;) for defining settings for gravity field variations. Each type fo
 * gravity field variations requires its own dedicated derived class, in which the properties to be
 * used for the specific model are to be defined.
 */
class GravityFieldVariationSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param bodyDeformationType Type of gravity field variation to be used.
     * \param interpolatorSettings Settings that are to be used to create an interpolator for the
     * gravity field variations immediately upon creation (to be used during propagation). Default
     * is nullptr, in which no interpolation is used, and the model is evaluated during propagation.
     */
    GravityFieldVariationSettings( const gravitation::BodyDeformationTypes bodyDeformationType,
                                   const std::shared_ptr< ModelInterpolationSettings > interpolatorSettings = nullptr ):
        bodyDeformationType_( bodyDeformationType ),
        interpolatorSettings_( interpolatorSettings ){ }

    //! Virtual destructor.
    virtual ~GravityFieldVariationSettings( ){ }

    //! Function to retrieve type of gravity field variation to be used.
    /*!
     * \brief Function to retrieve type of gravity field variation to be used.
     * \return Type of gravity field variation to be used.
     */
    gravitation::BodyDeformationTypes getBodyDeformationType( ){ return bodyDeformationType_;  }

    //! Function to retrieve settings that are to be used to create an interpolator for the gravity
    //! field variations
    /*!
     * \brief Function to retrieve settings that are to be used to create an interpolator for the
     * gravity field variations
     * \return Settings that are to be used to create an interpolator for the gravity field
     * variations (if not nullptr)
     */
    std::shared_ptr< ModelInterpolationSettings > getInterpolatorSettings( ){ return interpolatorSettings_; }

protected:

    //! Type of gravity field variation to be used.
    gravitation::BodyDeformationTypes bodyDeformationType_;

    //! Settings that are to be used to create an interpolator for the gravity field variations
    /*!
     * Settings that are to be used to create an interpolator for the gravity field variations
     * immediately upon creation. If nullptr, no interpolation is used, and the model is evaluated
     * during propagation.
     */
    std::shared_ptr< ModelInterpolationSettings > interpolatorSettings_;

};

//! Class to define settings for basic tidal gravity field variations, i.e. according to Eq. (6.6)
//! of IERS 2010 conventions.
class BasicSolidBodyGravityFieldVariationSettings: public GravityFieldVariationSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param deformingBodies List of bodies causing tidal deformation
     * \param loveNumbers List of Love number for the deformed body. First vector level denotes
     *  degree (index 0 = degree 2), second vector level denotes order and must be of maximum size
     *  (loveNumbers.size( ) + 2, i.e. maximum degree >= maximum order)
     * \param interpolatorSettings Settings that are to be used to create an interpolator for the
     * gravity field variations immediately upon creation (to be used during propagation). Default
     * is nullptr, in which no interpolation is used, and the model is evaluated during propagation.
     */
    BasicSolidBodyGravityFieldVariationSettings(
            const std::vector< std::string > deformingBodies,
            const std::map< int, std::vector< std::complex< double > > > loveNumbers,
            const std::shared_ptr< ModelInterpolationSettings > interpolatorSettings = nullptr ):
        GravityFieldVariationSettings( gravitation::basic_solid_body, interpolatorSettings ),
        deformingBodies_( deformingBodies ), loveNumbers_( loveNumbers ){ }

    virtual ~BasicSolidBodyGravityFieldVariationSettings( ){ }

    //! Function to retrieve list of bodies causing tidal deformation
    /*!
     * \brief Function to retrieve list of bodies causing tidal deformation
     * \return List of bodies causing tidal deformation
     */
    std::vector< std::string > getDeformingBodies( ){ return deformingBodies_;}

    //! Function to retrieve list of Love number for the deformed body.
    /*!
     * \brief Function to retrieve list of Love number for the deformed body.
     * \return List of Love number for the deformed body.
     */
    std::map< int, std::vector< std::complex< double > > > getLoveNumbers( ){ return loveNumbers_; }


    //! Function to reset list of bodies causing tidal deformation
    /*!
     * \brief Function to reset list of bodies causing tidal deformation
     * \param deformingBodies New list of bodies causing tidal deformation
     */
    void resetDeformingBodies( const std::vector< std::string >& deformingBodies ){
        deformingBodies_ = deformingBodies; }

protected:

    //! List of bodies causing tidal deformation
    std::vector< std::string > deformingBodies_;

    //! List of Love number for the deformed body.
    std::map< int, std::vector< std::complex< double > > > loveNumbers_;

};


class ModeCoupledSolidBodyGravityFieldVariationSettings: public GravityFieldVariationSettings
{
public:


    ModeCoupledSolidBodyGravityFieldVariationSettings(
        const std::vector< std::string > deformingBodies,
        const std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > loveNumbers ):
        GravityFieldVariationSettings( gravitation::mode_coupled_solid_body ),
        deformingBodies_( deformingBodies ), loveNumbers_( loveNumbers ){ }

    virtual ~ModeCoupledSolidBodyGravityFieldVariationSettings( ){ }

    //! Function to retrieve list of bodies causing tidal deformation
    /*!
     * \brief Function to retrieve list of bodies causing tidal deformation
     * \return List of bodies causing tidal deformation
     */
    std::vector< std::string > getDeformingBodies( ){ return deformingBodies_;}

    //! Function to retrieve list of Love number for the deformed body.
    /*!
     * \brief Function to retrieve list of Love number for the deformed body.
     * \return List of Love number for the deformed body.
     */
    std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > getLoveNumbers( ){ return loveNumbers_; }


    //! Function to reset list of bodies causing tidal deformation
    /*!
     * \brief Function to reset list of bodies causing tidal deformation
     * \param deformingBodies New list of bodies causing tidal deformation
     */
    void resetDeformingBodies( const std::vector< std::string >& deformingBodies ){
        deformingBodies_ = deformingBodies; }

protected:

    //! List of bodies causing tidal deformation
    std::vector< std::string > deformingBodies_;

    //! List of Love number for the deformed body.
    std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > loveNumbers_;

};

//! Class to define settings for tabulated gravity field variations.
class TabulatedGravityFieldVariationSettings: public GravityFieldVariationSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param cosineCoefficientCorrections Map of corrections to cosine coefficients (key is time)
     * \param sineCoefficientCorrections Map of corrections to sine coefficients (key is time)
     * \param minimumDegree Minimum degree (i.e. degree represented by first row of correction
     * matrix values in (co)sineCoefficientCorrections) of spherical harmonic corrections.
     * \param minimumOrder Minimum order (i.e. degree represented by first column of correction
     * matrix values in (co)sineCoefficientCorrections of spherical harmonic corrections.
     * \param interpolatorSettings Settings to use for the interpolator that is to be used.
     */
    TabulatedGravityFieldVariationSettings(
            const std::map< double, Eigen::MatrixXd > cosineCoefficientCorrections,
            const std::map< double, Eigen::MatrixXd > sineCoefficientCorrections,
            const int minimumDegree, const int minimumOrder,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings ):
        GravityFieldVariationSettings(
            gravitation::tabulated_variation, std::make_shared< ModelInterpolationSettings >(
                TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, interpolatorSettings ) ),
        cosineCoefficientCorrections_( cosineCoefficientCorrections ),
        sineCoefficientCorrections_( sineCoefficientCorrections ),
        minimumDegree_( minimumDegree ), minimumOrder_( minimumOrder ){ }

    //! Function to retrieve map of corrections to cosine coefficients (key is time)
    /*!
     * \brief Function to retrieve map of corrections to cosine coefficients (key is time)
     * \return Map of corrections to cosine coefficients (key is time)
     */
    std::map< double, Eigen::MatrixXd > getCosineCoefficientCorrections( )
    {
        return cosineCoefficientCorrections_;
    }

    //! Function to retrieve map of corrections to sine coefficients (key is time)
    /*!
     * \brief Function to retrieve map of corrections to sine coefficients (key is time)
     * \return Map of corrections to sine coefficients (key is time)
     */
    std::map< double, Eigen::MatrixXd > getSineCoefficientCorrections( )
    {
        return sineCoefficientCorrections_;
    }

    //! Function to retrieve minimum degree of spherical harmonic corrections.
    /*!
     * \brief Function to retrieve minimum degree of spherical harmonic corrections.
     * \return Minimum degree of spherical harmonic corrections.
     */
    int getMinimumDegree( )
    {
        return minimumDegree_;
    }

    //! Function to retrieve minimum order of spherical harmonic corrections.
    /*!
     * \brief Function to retrieve minimum order of spherical harmonic corrections.
     * \return Minimum order of spherical harmonic corrections.
     */
    int getMinimumOrder( )
    {
        return minimumOrder_;
    }
private:

    //! Map of corrections to cosine coefficients (key is time)
    std::map< double, Eigen::MatrixXd > cosineCoefficientCorrections_;

    //! Map of corrections to sine coefficients (key is time)
    std::map< double, Eigen::MatrixXd > sineCoefficientCorrections_;

    //! Minimum degree of spherical harmonic corrections.
    int minimumDegree_;

    //! Minimum order of spherical harmonic corrections.
    int minimumOrder_;

};

class PeriodicGravityFieldVariationsSettings: public GravityFieldVariationSettings
{

public:
    PeriodicGravityFieldVariationsSettings(
            const std::vector< Eigen::MatrixXd >& cosineShAmplitudesCosineTime,
            const std::vector< Eigen::MatrixXd >& cosineShAmplitudesSineTime,
            const std::vector< Eigen::MatrixXd >& sineShAmplitudesCosineTime,
            const std::vector< Eigen::MatrixXd >& sineShAmplitudesSineTime,
            const std::vector< double >& frequencies,
            const double referenceEpoch,
            const int minimumDegree = 2,
            const int minimumOrder = 0 ):
        GravityFieldVariationSettings( gravitation::periodic_variation ),
        cosineShAmplitudesCosineTime_( cosineShAmplitudesCosineTime ),
        cosineShAmplitudesSineTime_( cosineShAmplitudesSineTime ),
        sineShAmplitudesCosineTime_( sineShAmplitudesCosineTime ),
        sineShAmplitudesSineTime_( sineShAmplitudesSineTime ),
        frequencies_( frequencies ),
        referenceEpoch_( referenceEpoch ),
        minimumDegree_( minimumDegree ),
        minimumOrder_( minimumOrder ){ }

    std::vector< Eigen::MatrixXd > getCosineShAmplitudesCosineTime( )
    {
        return cosineShAmplitudesCosineTime_;
    }

    std::vector< Eigen::MatrixXd > getCosineShAmplitudesSineTime( )
    {
        return cosineShAmplitudesSineTime_;
    }

    std::vector< Eigen::MatrixXd > getSineShAmplitudesCosineTime( )
    {
        return sineShAmplitudesCosineTime_;
    }

    std::vector< Eigen::MatrixXd > getSineShAmplitudesSineTime( )
    {
        return sineShAmplitudesSineTime_;
    }

    std::vector< double > getFrequencies( )
    {
        return frequencies_;
    }

    double getReferenceEpoch( )
    {
        return referenceEpoch_;
    }

    int getMinimumDegree( )
    {
        return minimumDegree_;
    }

    int getMinimumOrder( )
    {
        return minimumOrder_;
    }

private:

    std::vector< Eigen::MatrixXd > cosineShAmplitudesCosineTime_;

    std::vector< Eigen::MatrixXd > cosineShAmplitudesSineTime_;

    std::vector< Eigen::MatrixXd > sineShAmplitudesCosineTime_;

    std::vector< Eigen::MatrixXd > sineShAmplitudesSineTime_;

    std::vector< double > frequencies_;

    double referenceEpoch_;

    int minimumDegree_;

    int minimumOrder_;

};

class PolynomialGravityFieldVariationsSettings: public GravityFieldVariationSettings
{

public:
    PolynomialGravityFieldVariationsSettings(
        const std::map<int, Eigen::MatrixXd> &cosineAmplitudes,
        const std::map<int, Eigen::MatrixXd> &sineAmplitudes,
        const double referenceEpoch,
        const int minimumDegree,
        const int minimumOrder ):
        GravityFieldVariationSettings( gravitation::polynomial_variation ),
        cosineAmplitudes_( cosineAmplitudes ),
        sineAmplitudes_( sineAmplitudes ),
        referenceEpoch_( referenceEpoch ),
        minimumDegree_( minimumDegree ),
        minimumOrder_( minimumOrder ){ }

    std::map< int, Eigen::MatrixXd > getCosineAmplitudes( )
    {
        return cosineAmplitudes_;
    }

    std::map< int, Eigen::MatrixXd > getSineAmplitudes( )
    {
        return sineAmplitudes_;
    }


    double getReferenceEpoch( )
    {
        return referenceEpoch_;
    }

    int getMinimumDegree( )
    {
        return minimumDegree_;
    }

    int getMinimumOrder( )
    {
        return minimumOrder_;
    }

private:
    std::map< int, Eigen::MatrixXd > cosineAmplitudes_;

    std::map< int, Eigen::MatrixXd > sineAmplitudes_;

    double referenceEpoch_;

    int minimumDegree_;

    int minimumOrder_;
};
//! Function to create constant complex Love number list for a range of degrees and orders.
/*!
 * Function to create constant complex Love number list for a range of degrees and orders, maximum degree and order
 * are given as input, minimum degree and order are 2 and 0, respectively.
 * \param constantLoveNumber Love number to be set at each degree and order
 * \param maximumDegree Maximum degree for Love numbers.
 * \param maximumOrder Maximum order for Love numbers.
 * \return List of Love numbers with requested settings.
 */
std::map< int, std::vector< std::complex< double > > > getFullLoveNumbersVector(
        const std::complex< double > constantLoveNumber, const int maximumDegree, const int maximumOrder );

//! Function to create constant real Love number list for a range of degrees and orders
/*!
* Function to create constant real Love number list for a range of degrees and orders, maximum degree and order
* are given as input, minimum degree and order are 2 and 0, respectively.
* \param constantLoveNumber Love number to be set at each degree and order
* \param maximumDegree Maximum degree for Love numbers.
* \param maximumOrder Maximum order for Love numbers.
* \return List of Love numbers with requested settings.
*/
std::map< int, std::vector< std::complex< double > > > getFullLoveNumbersVector(
        const double constantLoveNumber, const int maximumDegree, const int maximumOrder );


std::vector< std::complex< double > > getLoveNumberPerDegree(
        const std::complex< double > loveNumber,
        const int degree );

std::vector< std::complex< double > > getLoveNumberPerDegree(
        const double loveNumber,
        const int degree );

inline std::shared_ptr< GravityFieldVariationSettings > fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
        const std::string deformingBody,
        const double loveNumber,
        const int degree )
{
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    loveNumbers[ degree ] = getLoveNumberPerDegree( loveNumber, degree );

    return std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >( { deformingBody } ), loveNumbers, nullptr );
}

inline std::shared_ptr< GravityFieldVariationSettings > fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
        const std::string deformingBody,
        const std::complex< double > loveNumber,
        const int degree   )
{
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    loveNumbers[ degree ] = getLoveNumberPerDegree( loveNumber, degree );

    return std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >( { deformingBody } ), loveNumbers, nullptr );
}

inline std::shared_ptr< GravityFieldVariationSettings > fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
        const std::string deformingBody,
        const std::map< int, double > loveNumberPerDegree )
{
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    for( auto loveNumberIt : loveNumberPerDegree )
    {
        loveNumbers[ loveNumberIt.first ] = getLoveNumberPerDegree( loveNumberIt.second, loveNumberIt.first );
    }

    return std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >( { deformingBody } ), loveNumbers, nullptr );
}

inline std::shared_ptr< GravityFieldVariationSettings > fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
        const std::string deformingBody,
        const std::map< int, std::complex< double > > loveNumberPerDegree  )
{
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    for( auto loveNumberIt : loveNumberPerDegree )
    {
        loveNumbers[ loveNumberIt.first ] = getLoveNumberPerDegree( loveNumberIt.second, loveNumberIt.first );
    }

    return std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >( { deformingBody } ), loveNumbers, nullptr );
}


inline std::shared_ptr< GravityFieldVariationSettings > orderVariableSingleDegreeLoveNumberGravityFieldVariationSettings(
        const std::string deformingBody,
        const std::vector< double > loveNumber,
        const int degree,
        const std::shared_ptr< ModelInterpolationSettings > interpolatorSettings = nullptr )
{
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    for( unsigned int i = 0; i < loveNumber.size( ); i++ )
    {
        loveNumbers[ degree ].push_back( std::complex< double >( loveNumber.at( i ), 0 ) );
    }
    return std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >( { deformingBody } ), loveNumbers, interpolatorSettings );
}

inline std::shared_ptr< GravityFieldVariationSettings > orderVariableSingleDegreeLoveNumberGravityFieldVariationSettings(
        const std::string deformingBody,
        const std::vector< std::complex< double > > loveNumber,
        const int degree,
        const std::shared_ptr< ModelInterpolationSettings > interpolatorSettings = nullptr )
{
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    loveNumbers[ degree ] = loveNumber;
    return std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >( { deformingBody } ), loveNumbers, interpolatorSettings );
}

inline std::shared_ptr< GravityFieldVariationSettings > degreeOrderVariableLoveNumberGravityFieldVariationSettings(
        const std::string deformingBody,
        const std::map< int, std::vector< double > > loveNumber,
        const std::shared_ptr< ModelInterpolationSettings > interpolatorSettings = nullptr )
{
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    for( auto loveNumberIt : loveNumber )
    {
        for( unsigned int i = 0; i < loveNumberIt.second.size( ); i++ )
        {
            loveNumbers[ loveNumberIt.first ].push_back( std::complex< double >( loveNumberIt.second.at( i ), 0 ) );
        }
    }
    return std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >( { deformingBody } ), loveNumbers, interpolatorSettings );
}

inline std::shared_ptr< GravityFieldVariationSettings > degreeOrderVariableLoveNumberGravityFieldVariationSettings(
        const std::string deformingBody,
        const std::map< int, std::vector< std::complex< double > > > loveNumber,
        const std::shared_ptr< ModelInterpolationSettings > interpolatorSettings = nullptr )
{
    return std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >( { deformingBody } ), loveNumber, interpolatorSettings );
}

inline std::shared_ptr< GravityFieldVariationSettings > tabulatedGravityFieldVariationSettings(
        const std::map< double, Eigen::MatrixXd > cosineCoefficientCorrections,
        const std::map< double, Eigen::MatrixXd > sineCoefficientCorrections,
        const int minimumDegree, const int minimumOrder,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    return std::make_shared< TabulatedGravityFieldVariationSettings >(
                cosineCoefficientCorrections, sineCoefficientCorrections, minimumDegree, minimumOrder,
                interpolatorSettings );
}

inline std::shared_ptr< GravityFieldVariationSettings > periodicGravityFieldVariationsSettings(
        const std::vector< Eigen::MatrixXd >& cosineShAmplitudesCosineTime,
        const std::vector< Eigen::MatrixXd >& cosineShAmplitudesSineTime,
        const std::vector< Eigen::MatrixXd >& sineShAmplitudesCosineTime,
        const std::vector< Eigen::MatrixXd >& sineShAmplitudesSineTime,
        const std::vector< double >& frequencies,
        const double referenceEpoch,
        const int minimumDegree = 2,
        const int minimumOrder = 0 )
{
    return std::make_shared< PeriodicGravityFieldVariationsSettings >(
        cosineShAmplitudesCosineTime, cosineShAmplitudesSineTime, sineShAmplitudesCosineTime, sineShAmplitudesSineTime, frequencies, referenceEpoch, minimumDegree, minimumOrder );
}

inline std::shared_ptr< GravityFieldVariationSettings > periodicGravityFieldVariationsSettingsSingleFrequency(
        const Eigen::MatrixXd& cosineShAmplitudesCosineTime,
        const Eigen::MatrixXd& cosineShAmplitudesSineTime,
        const Eigen::MatrixXd& sineShAmplitudesCosineTime,
        const Eigen::MatrixXd& sineShAmplitudesSineTime,
        const double frequency,
        const double referenceEpoch,
        const int minimumDegree = 2,
        const int minimumOrder = 0 )
{
    return std::make_shared< PeriodicGravityFieldVariationsSettings >(
                std::vector< Eigen::MatrixXd >( { cosineShAmplitudesCosineTime } ),
                std::vector< Eigen::MatrixXd >( { cosineShAmplitudesSineTime } ),
                std::vector< Eigen::MatrixXd >( { sineShAmplitudesCosineTime } ),
                std::vector< Eigen::MatrixXd >( { sineShAmplitudesSineTime } ),
                std::vector< double >( { frequency } ),
                referenceEpoch, minimumDegree, minimumOrder );

}

inline std::shared_ptr< GravityFieldVariationSettings > polynomialGravityFieldVariationsSettings(
    const std::map< int, Eigen::MatrixXd >& cosineAmplitudes,
    const std::map< int, Eigen::MatrixXd >& sineAmplitudes,
    const double referenceEpoch,
    const int minimumDegree = 2,
    const int minimumOrder = 0 )
{
    return std::make_shared< PolynomialGravityFieldVariationsSettings >(
        cosineAmplitudes, sineAmplitudes, referenceEpoch, minimumDegree, minimumOrder );
}

inline std::shared_ptr< GravityFieldVariationSettings > polynomialGravityFieldVariationsSettingsSinglePower(
    const Eigen::MatrixXd& cosineAmplitude,
    const Eigen::MatrixXd& sineAmplitude,
    const int polynomialPower,
    const double referenceEpoch,
    const int minimumDegree = 2,
    const int minimumOrder = 0 )
{
    std::map< int, Eigen::MatrixXd > cosineAmplitudes;
    cosineAmplitudes[ polynomialPower ] = cosineAmplitude;
    std::map< int, Eigen::MatrixXd > sineAmplitudes;
    sineAmplitudes[ polynomialPower ] = sineAmplitude;

    return std::make_shared< PolynomialGravityFieldVariationsSettings >(
        cosineAmplitudes, sineAmplitudes, referenceEpoch, minimumDegree, minimumOrder );
}

inline std::shared_ptr< GravityFieldVariationSettings > modeCoupledSolidBodyGravityFieldVariationSettings(
    const std::vector< std::string > deformingBodies,
    const std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > loveNumbers )
{
    return std::make_shared< ModeCoupledSolidBodyGravityFieldVariationSettings >( deformingBodies, loveNumbers );
}


//! Function to create a set of gravity field variations, stored in the associated interface class
/*!
 * Function to create a set of gravity field variations, stored in the associated interface class of
 * type GravityFieldVariationsSet
 * \param body Body for which gravity field variations are createad
 * \param bodies List of body objects in simulations.
 * \param gravityFieldVariationSettings List of settings for gravity field variations
 * \return Interface class containing list of GravityFieldVariations.
 */
std::shared_ptr< gravitation::GravityFieldVariationsSet > createGravityFieldModelVariationsSet(
        const std::string& body,
        const SystemOfBodies& bodies,
        const std::vector< std::shared_ptr< GravityFieldVariationSettings > >&
        gravityFieldVariationSettings );

//! Function to create a single gravity field variation object.
/*!
 * Function to create a single gravity field variation object.
 * \param gravityFieldVariationSettings Settings for the gravity field variation object that is to
 * be created. Depending on the settings and type of the gravity field variations, spherical
 * harmonic corrections are calculated a priori and handled by an interpolator during propagation,
 * or they are directly calculated from the current state during numerical propagation.
 * \param body Body for which gravity field variations are createad
 * \param bodies List of body objects in simulations.
 * \return Single gravity field variation object.
 */
std::shared_ptr< gravitation::GravityFieldVariations > createGravityFieldVariationsModel(
        const std::shared_ptr< GravityFieldVariationSettings > gravityFieldVariationSettings,
        const std::string body,
        const SystemOfBodies& bodies );

} // namespace simulation_setup
} // namespace tudat
#endif // TUDAT_CREATEGRAVITYFIELDVARIATIONS_H
