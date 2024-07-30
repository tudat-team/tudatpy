/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_BASICSOLIDBODYTIDEGRAVITYFIELDVARIATIONS_H
#define TUDAT_BASICSOLIDBODYTIDEGRAVITYFIELDVARIATIONS_H

#include <functional>
#include <boost/multi_array.hpp>


#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vector>
#include <complex>
#include <map>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"



namespace tudat
{

namespace gravitation
{


//! Function to calculate solid body tide gravity field variations due to single body at single degree and order from
//! precomputed quantaties.
/*!
 *  Function to calculate solid body tide gravity field variations due to single body at single
 *  degree and order, (frequency-independent part), after (Petit et al. 2010, eq. 6.6), using
 *  complex tidal love numbers. Using this function requires the computation of a number of precomputed quantities
 *  (as is done by the BasicSolidBodyTideGravityFieldVariations class).
 *  \param loveNumber Complex Love number for given degree and order.
 *  \param massRatio Ratio of masses of body causing deformation to body being deformed.
 *  \param radiusRatioPowerN Ratio of equatorial radius of body being deformed over distance
 *  between center of body being deformed to center of body causing deformation to the power (n+1).
 *  \param amplitude Amplitude of the tide.
 *  \param tideArgument Argument of the tide.
 *  \param degree Degree of current coefficients
 *  \param order Order of current coefficients
 *  \return Combined variation in cosine (Delta C_{n,m}) and sine (Delta S_{n,m}) coefficients as:
 *  Delta C_{n,m} - i * Delta S_{n,m}
 */
std::complex< double > calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
        const std::complex< double > loveNumber, const double massRatio,
        const double radiusRatioPowerN, const double amplitude,
        const std::complex< double > tideArgument, const int degree, const int order );

//! Function to calculate solid body tide gravity field variations due to single body at single degree and order directly
//! from perturbing body's Cartesian state.
/*!
 *  Function to calculate solid body tide gravity field variations due to single body at single
 *  degree and order, (frequency-independent part), after (Petit et al. 2010, eq. 6.6), using
 *  complex tidal love numbers. This function directly uses the perturbing body's Cartesian state, and should be used for
 *  'external' computations. The BasicSolidBodyTideGravityFieldVariations uses the alternative overloaded function.
 *  \param loveNumber Complex Love number for given degree and order.
 *  \param massRatio Ratio of masses of body causing deformation to body being deformed.
 *  \param referenceRadius Equatorial radius of body being deformed.
 *  \param relativeBodyFixedPosition Cartesian position of body causing deformation in a frame centered on and fixed to
 *  the body that is being deformed.
 *  \param degree Degree of current coefficients
 *  \param order Order of current coefficients
 *  \return Combined variation in cosine (Delta C_{n,m}) and sine (Delta S_{n,m}) coefficients as:
 *  Delta C_{n,m} - i * Delta S_{n,m}
 */
std::complex< double > calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
        const std::complex< double > loveNumber, const double massRatio,
        const double referenceRadius, const Eigen::Vector3d& relativeBodyFixedPosition, const int degree, const int order );

//! Function to calculate solid body tide gravity field variations due to single body at a set of degrees and orders
//! from perturbing body's Cartesian state.
/*!
 *  Function to calculate solid body tide gravity field variations due to single body at a set of degrees and orders
 *  (frequency-independent part), after (Petit et al. 2010, eq. 6.6), using complex tidal love numbers.
 *  This function directly uses the perturbing body's Cartesian state, and should be used for
 *  'external' computations. The BasicSolidBodyTideGravityFieldVariations uses the alternative overloaded function.
 *  This function defines the maximum degree and order, and always starts at degree 2 and order 0
 *  \param loveNumbers Complex Love numbers for each degree and order (index of first vector is degree-2; index of second
 *  vector is order.
 *  \param massRatio Ratio of masses of body causing deformation to body being deformed.
 *  \param referenceRadius Equatorial radius of body being deformed.
 *  \param relativeBodyFixedPosition Cartesian position of body causing deformation in a frame centered on and fixed to
 *  the body that is being deformed.
 *  \param maximumDegree Maximum degree of current coefficient corrections.
 *  \param maximumOrder Maximum order of current coefficient corrections.
 *  \return Combined variation in cosine (Delta C_{n,m}) and sine (Delta S_{n,m}) matrices as first and second entry of
 *  pair, respectively. Both matrices start at degree and order 0.
 */
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
        const std::map< int, std::vector< std::complex< double > > > loveNumbers, const double massRatio,
        const double referenceRadius, const Eigen::Vector3d& relativeBodyFixedPosition,
        const int maximumDegree, const int maximumOrder );



class SolidBodyTideGravityFieldVariations: public GravityFieldVariations
{
public:
    //! Constructor
    /*!
     *  Contructor, sets Love numbers and properties of deformed and tide-raising bodies.
     *  \param deformedBodyStateFunction Function returning state of body being deformed.
     *  \param deformedBodyOrientationFunction Function providing rotation from inertial to body
     *  being deformed-fixed frame
     *  \param deformingBodyStateFunctions List of state functions of body causing deformations.
     *  \param deformedBodyReferenceRadius Reference radius (typically equatorial) of body being
     *  deformed's spherical harmonic gravity field.
     *  \param deformedBodyMass Function returning mass (or gravitational parameter) of body being
     *  deformed.
     *  \param deformingBodyMasses List of functions returning masses (or gravitational parameters)
     *  of bodies causing deformation.
     *  \param deformingBodies List of names of bodies causing deformation
     */
    SolidBodyTideGravityFieldVariations(
        const std::function< Eigen::Vector6d( const double ) >
        deformedBodyStateFunction,
        const std::function< Eigen::Quaterniond( const double ) >
        deformedBodyOrientationFunction,
        const std::vector< std::function< Eigen::Vector6d( const double ) > >
        deformingBodyStateFunctions,
        const double deformedBodyReferenceRadius,
        const std::function< double( ) > deformedBodyMass,
        const std::vector< std::function< double( ) > > deformingBodyMasses,
        const std::vector< std::string > deformingBodies,
        const int maximumDegree, const int maximumOrder ):
    // LOVE NUMBERS TODO: FIX MIN/MAX STUFF
        GravityFieldVariations( 2, 0, maximumDegree, maximumOrder ),
        deformedBodyStateFunction_( deformedBodyStateFunction ),
        deformedBodyOrientationFunction_( deformedBodyOrientationFunction ),
        deformingBodyStateFunctions_( deformingBodyStateFunctions ),
        deformedBodyReferenceRadius_( deformedBodyReferenceRadius ),
        deformedBodyMass_( deformedBodyMass ),
        deformingBodyMasses_( deformingBodyMasses ),
        deformingBodies_( deformingBodies )
    {
        currentCosineCorrections_ = Eigen::MatrixXd::Zero(
            maximumDegree_ - minimumDegree_ + 1, maximumOrder_ - minimumOrder_ + 1 );
        currentSineCorrections_ = Eigen::MatrixXd::Zero(
            maximumDegree_ - minimumDegree_ + 1, maximumOrder_ - minimumOrder_ + 1 );
    }

    //! Destructor
    /*!
     *  Destructor
     */
    virtual ~SolidBodyTideGravityFieldVariations( ){ }

    //! Function for calculating basic spherical harmonic coefficient corrections.
    /*!
     *  Function for calculating basic spherical harmonic coefficient corrections.
     *  \param time Time at which variations are to be calculated.
     *  \return Pair of matrices containing variations in (cosine,sine) coefficients.
     */
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateBasicSphericalHarmonicsCorrections(
        const double time );

    //! Derived function for calculating spherical harmonic coefficient corrections.
    /*!
     *  Derived function for calculating spherical harmonic coefficient corrections.
     *  \param time Time at which variations are to be calculated.
     *  \return Pair of matrices containing variations in (cosine,sine) coefficients.
     */
    virtual std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections(
        const double time )
    {
        return calculateBasicSphericalHarmonicsCorrections( time );
    }

    //! Function to return reference radius the spherical harmonic gravity field of deformed body.
    /*!
     *  Function to return reference radius (typically equatorial)  the spherical harmonic gravity
     * field of deformed body.
     *  \return Reference radius body being deformed's spherical harmonic gravity field.
     */
    double getDeformedBodyReferenceRadius( )
    {
        return deformedBodyReferenceRadius_;
    }

    //! Function to return the mass function of the deformed body.
    /*!
     *  Function to return the mass function of the deformed body.
     *  \return Mass function of the deformed body.
     */
    std::function< double( ) > getDeformedBodyMassFunction( )
    {
        return deformedBodyMass_;
    }

    //! Function to return list of the mass functions of the bodies causing the deformation.
    /*!
     *  Function to return list of the mass functions of the bodies causing the deformation.
     *  \return List of the mass functions of the bodies causing the deformation.
     */
    std::vector< std::function< double( ) > > getDeformingBodyMasses( )
    {
        return deformingBodyMasses_;
    }

    //! Function to return list of the names of the bodies causing the deformation
    /*!
     *  Function to return list of the names of the bodies causing the deformation
     *  \return List of the names of the bodies causing the deformation
     */
    std::vector< std::string > getDeformingBodies( )
    {
        return deformingBodies_;
    }

    //! Function to return the state function of the deformed body.
    /*!
     *  Function to return the state function of the deformed body.
     *  \return State function of the deformed body.
     */
    std::function< Eigen::Vector6d( const double ) > getDeformedBodyStateFunction( )
    {
        return deformedBodyStateFunction_;
    }

    //! Function to return the orientation function (rotation to body-fixed frame) of the deformed
    //! body.
    /*!
     *  Function to return the orientation function (rotation to body-fixed frame) of the deformed
     *  body.
     *  \return Orientation function (rotation to body-fixed frame) of the deformed body.
     */
    std::function< Eigen::Quaterniond( const double ) > getDeformedBodyOrientationFunction( )
    {
        return deformedBodyOrientationFunction_;
    }

    //! Function to return list of the state functions of the bodies causing the deformation.
    /*!
     *  Function to return list of the state functions of the bodies causing the deformation.
     *  \return List of the state functions of the bodies causing the deformation.
     */
    std::vector< std::function< Eigen::Vector6d( const double ) > >
    getDeformingBodyStateFunctions( )
    {
        return deformingBodyStateFunctions_;
    }

    //! Get a string with the concatenation of all the bodies causing the deformation.
    /*!
     *  Get a string with the concatenation of all the bodies causing the deformation.
     *  \return A string with the concatenation of all the bodies causing the deformation.
     */
    std::string getConcatenatedDeformingBodies( )
    {
        // Initialize string
        std::string id = "";

        // Iterate over all bodies and add name to concatenates string
        for( unsigned int i = 0; i < deformingBodies_.size( ); i++ )
        {
            id += deformingBodies_[ i ];
            if( i != deformingBodies_.size( ) - 1 )
            {
                id += "_";
            }
        }

        return id;
    }

    //! Function to return the current corrections to the cosine coefficients.
    /*!
     * Function to return the current corrections to the cosine coefficients (as calculated by
     * current instance of this class).
     * \return Current corrections to the cosine coefficients.
     */
    Eigen::MatrixXd getCurrentCosineCorrections( )
    {
        return currentCosineCorrections_;
    }

    //! Function to return the current corrections to the sine coefficients.
    /*!
     * Function to return the current corrections to the sine coefficients (as calculated by
     * current instance of this class).
     * \return Current corrections to the sine coefficients.
     */
    Eigen::MatrixXd getCurrentSineCorrections( )
    {
        return currentSineCorrections_;
    }

protected:

    //! List of functions to call for calculating spherical harmonic corrections.
    /*!
     *  List of functions to call for calculating spherical harmonic corrections. Each function
     *  modifies MatrixXd arguments (cosine, sine) as they are passed by reference, and adds
     *  the required (tidal) correction.
     */
    std::vector< std::function< void( Eigen::MatrixXd&, Eigen::MatrixXd& ) > >
        correctionFunctions;

    //! Calculates basic solid body gravity field corrections due to single body.
    /*!
     *  Calculates basic solid body gravity field corrections for all degrees and orders set.
     *  The arguments are modified as they are passed by reference, through which the corrections
     *  are returned.
     *  Class variables denoting properties of currently considered body must have been set before
     *  this function is called.
     *  \param cTermCorrections Corrections to cosine terms
     *  (passed by reference; correction added to input value).
     *  \param sTermCorrections Corrections to sine terms.
     *
     *  (passed by reference; correction added to input value).
     */
    virtual void addBasicSolidBodyTideCorrections(
        Eigen::MatrixXd& cTermCorrections, Eigen::MatrixXd& sTermCorrections ) = 0;

    //! Sets current properties (mass state) of body causing tidal deformation.
    /*!
     *  Sets current properties (mass state) of body causing tidal deformation.
     * \param bodyIndex Index of body causing deformation for which data is to be retrieved.
     * Deformed body is also updated to bodyIndex = 0.
     * \param evaluationTime Time at which properties are to be evaluated.
     */
    virtual void setBodyGeometryParameters(
        const int bodyIndex, const double evaluationTime);

    //! Calculate tidal amplitude and argument at current degree and order.
    /*!
     * Calculate tidal amplitude and argument at current degree and order.
     * \param degree Degree of tide.
     * \param order Order of tide.
     */
    virtual void updateTidalAmplitudeAndArgument(
        const int degree, const int order )
    {
        tideAmplitude = basic_mathematics::computeLegendrePolynomialExplicit(
            degree, order, sineOfLatitude );
        tideArgument = static_cast< double >( order ) * iLongitude;
    }


    //! Function returning state of body being deformed.
    /*!
     *  Function returning state of body being deformed.
     */
    std::function< Eigen::Vector6d( const double ) > deformedBodyStateFunction_;

    //! Function providing rotation from inertial to body being deformed-fixed frame
    /*!
     *  Function providing rotation from inertial to body being deformed-fixed frame.
     */
    std::function< Eigen::Quaterniond( const double ) > deformedBodyOrientationFunction_;

    //! List of state functions of body causing deformations.
    /*!
     *  List of state functions of body causing deformations.
     */
    std::vector< std::function< Eigen::Vector6d( const double ) > >
        deformingBodyStateFunctions_;

    //! Reference radius (typically equatorial) of body being deformed's spherical harmonic
    //! gravity field.
    /*!
     *  Reference radius (typically equatorial) of body being deformed's spherical harmonic gravity
     *  field.
     */
    double deformedBodyReferenceRadius_;

    //! Function returning mass of body being deformed.
    /*!
     *  Function returning mass of body being deformed.
     */
    std::function< double( ) > deformedBodyMass_;

    //! List of functions returning masses of bodies causing deformation.
    /*!
     *  List of functions returning masses of bodies causing deformation.
     */
    std::vector< std::function< double( ) > > deformingBodyMasses_;

    //! List of names of bodies causing deformation.
    /*!
     *  List of names of bodies causing deformation.
     */
    std::vector< std::string > deformingBodies_;


    //! Ratio of masses in current calculation step
    /*!
     *  Ratio of masses of body causing deformation to body being deformed in current calculation
     *  step.
     */
    double massRatio;

    //! Ratio of radii in current calculation step.
    /*!
     *  Ratio of equatorial radius of body being deformed over distance between center of body
     *  being deformed to center of body causing deformation in current calculation step
     */
    double radiusRatio;

    //! Ratio of radii in current calculation step to the power (degree+1).
    /*!
     *  Ratio of equatorial radius of body being deformed over distance between center of body
     *  being deformed to center of body causing deformation, to the power (degree+1), in current
     *  calculation step
     */
    double radiusRatioPower;

    //! Sine of latitude of currently considered body in current calculation step
    /*!
     *  Sine of latitude of body causing deformation in frame fixed to body being deformed in
     *  current calculation step.
     */
    double sineOfLatitude;

    //! Longitude of currently considered body times i in current calculation step
    /*!
     *  i (sqrt(-1)) times longitude of body causing deformation in frame fixed to body being
     *  deformed in current calculation step.
     */
    std::complex< double > iLongitude;

    //! Current argument of the tide.
    /*!
     *  Current argument of the tide, to be used as input to
     *  calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude, modulates the tide as
     *  exp( -tideArgument).
     */
    std::complex< double > tideArgument;

    //! Current amplitude of the tide.
    /*!
     *  Current amplitude of the tide, to be used as input to
     *  calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude.
     */
    double tideAmplitude;

    //! Current position of body being deformed.
    Eigen::Vector3d deformedBodyPosition;

    //! Current rotation to frame fixed to deformed body.
    Eigen::Quaterniond toDeformedBodyFrameRotation;

    //! Tidal corrections to cosine coefficients at current calculation step.
    Eigen::MatrixXd currentCosineCorrections_;

    //! Tidal corrections to sine coefficients at current calculation step.
    Eigen::MatrixXd currentSineCorrections_;

};


//! Class to calculate first-order solid body tide gravity field variations on a single body raised
//! by any number of bodies up to any degree and order.
class BasicSolidBodyTideGravityFieldVariations: public SolidBodyTideGravityFieldVariations
{
public:
    //! Constructor
    /*!
     *  Contructor, sets Love numbers and properties of deformed and tide-raising bodies.
     *  \param deformedBodyStateFunction Function returning state of body being deformed.
     *  \param deformedBodyOrientationFunction Function providing rotation from inertial to body
     *  being deformed-fixed frame
     *  \param deformingBodyStateFunctions List of state functions of body causing deformations.
     *  \param deformedBodyReferenceRadius Reference radius (typically equatorial) of body being
     *  deformed's spherical harmonic gravity field.
     *  \param deformedBodyMass Function returning mass (or gravitational parameter) of body being
     *  deformed.
     *  \param deformingBodyMasses List of functions returning masses (or gravitational parameters)
     *  of bodies causing deformation.
     *  \param loveNumbers List of love numbers for each degree and order. First vector level
     *  denotes degree (index 0 = degree 2), second vector level denotes order and must be of maximum size
     *  (loveNumbers.size( ) + 2, i.e. maximum degree >= maximum order)
     *  \param deformingBodies List of names of bodies causing deformation
     */
    BasicSolidBodyTideGravityFieldVariations(
            const std::function< Eigen::Vector6d( const double ) >
            deformedBodyStateFunction,
            const std::function< Eigen::Quaterniond( const double ) >
            deformedBodyOrientationFunction,
            const std::vector< std::function< Eigen::Vector6d( const double ) > >
            deformingBodyStateFunctions,
            const double deformedBodyReferenceRadius,
            const std::function< double( ) > deformedBodyMass,
            const std::vector< std::function< double( ) > > deformingBodyMasses,
            const std::map< int, std::vector< std::complex< double > > > loveNumbers,
            const std::vector< std::string > deformingBodies ):
        // LOVE NUMBERS TODO: FIX MIN/MAX STUFF
        SolidBodyTideGravityFieldVariations(
            deformedBodyStateFunction, deformedBodyOrientationFunction, deformingBodyStateFunctions, deformedBodyReferenceRadius,
            deformedBodyMass, deformingBodyMasses, deformingBodies, loveNumbers.rbegin( )->first, loveNumbers.rbegin( )->first ),
        loveNumbers_( loveNumbers )
    {
        // Set basic deformation functon as function to be evaluated when requesting variations.
        correctionFunctions.push_back(
            std::bind(
                &BasicSolidBodyTideGravityFieldVariations::addBasicSolidBodyTideCorrections,
                this, std::placeholders::_1, std::placeholders::_2 ));
    }

    //! Destructor
    /*!
     *  Destructor
     */
    virtual ~BasicSolidBodyTideGravityFieldVariations( ){ }


    //! Function to retrieve the love numbers at given degree.
    /*!
     *  Function to retrieve the love numbers at given degree. Returns a vector containing (complex)
     *  love numbers at all orders in current degree.
     *  \param degree Degree from which love numbers are to be retrieved.
     *  \return Vector of love numbers (i^{th} entry representing i^{th} order in requested degree)
     *  containing love numbers at current degree.
     */
    std::vector< std::complex< double > > getLoveNumbersOfDegree( const int degree )
    {
        return loveNumbers_[ degree ];
    }

    //! Function to return all love numbers.
    /*!
     *  Function to return all love numbers, i.e. at all degrees and orders.
     *  \return Complete set of available love numbers.
     */
    const std::map< int, std::vector< std::complex< double > > > getLoveNumbers( )
    {
        return loveNumbers_;
    }

    //! Function to reset the love numbers at given degree.
    /*!
     *  Function to reset the love numbers at given degree. Input requires a vector containing
     *  (complex) love numbers at all orders in current degree.
     *  \param degree Degree from which love numbers are to be retrieved.
     *  \param loveNumbers Vector of love numbers (i^{th} entry representing i^{th} order in requested degree)
     *  containing new love numbers at current degree.
     */
    void resetLoveNumbersOfDegree( const std::vector< std::complex< double > > loveNumbers,
                                   const int degree )
    {
        if( loveNumbers_.count( degree ) == 0 )
        {
            std::string errorMessage = "Error, tried to set love numbers at degree " +
                    std::to_string( degree ) +
                    " in BasicSolidBodyTideGravityFieldVariations: not available";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            if( loveNumbers.size( ) <= static_cast< unsigned int >( degree + 1 ) )
            {
                loveNumbers_[ degree ] = loveNumbers;
            }
            else
            {
                std::string errorMessage = "Error, tried to set love numbers at degree " +
                        std::to_string( degree ) + " in BasicSolidBodyTideGravityFieldVariations with" +
                        std::to_string( loveNumbers.size( ) ) + " orders";
                throw std::runtime_error( errorMessage );
            }
        }
    }

protected:


    //! Calculates basic solid body gravity field corrections due to single body.
    /*!
     *  Calculates basic solid body gravity field corrections for all degrees and orders set.
     *  The arguments are modified as they are passed by reference, through which the corrections
     *  are returned.
     *  Class variables denoting properties of currently considered body must have been set before
     *  this function is called.
     *  \param cTermCorrections Corrections to cosine terms
     *  (passed by reference; correction added to input value).
     *  \param sTermCorrections Corrections to sine terms.
     *
     *  (passed by reference; correction added to input value).
     */
    virtual void addBasicSolidBodyTideCorrections(
            Eigen::MatrixXd& cTermCorrections, Eigen::MatrixXd& sTermCorrections );


    // LOVE NUMBER TODO: FIX DOCS
    //! List of love numbers for each degree and order
    /*!
     *  List of love numbers for each degree and order. First vector level denotes degree
     *  (index 0 = degree 2), second vector level must be of size (loveNumbers.size( ) + 2, i.e.
     *  maximum degree == maximum order.
     */
    std::map< int, std::vector< std::complex< double > > > loveNumbers_;

};

int getModeCoupledMaximumResponseDegree(
    const std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > >& loveNumbers );

int getModeCoupledMaximumResponseOrder(
    const std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > >& loveNumbers );


class ModeCoupledSolidBodyTideGravityFieldVariations: public SolidBodyTideGravityFieldVariations
{
public:

    ModeCoupledSolidBodyTideGravityFieldVariations(
        const std::function< Eigen::Vector6d( const double ) > deformedBodyStateFunction,
        const std::function< Eigen::Quaterniond( const double ) > deformedBodyOrientationFunction,
        const std::vector< std::function< Eigen::Vector6d( const double ) > > deformingBodyStateFunctions,
        const double deformedBodyReferenceRadius,
        const std::function< double( ) > deformedBodyMass,
        const std::vector< std::function< double( ) > > deformingBodyMasses,
        const std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > loveNumbers,
        const std::vector< std::string > deformingBodies ):SolidBodyTideGravityFieldVariations(
        deformedBodyStateFunction, deformedBodyOrientationFunction, deformingBodyStateFunctions, deformedBodyReferenceRadius,
        deformedBodyMass, deformingBodyMasses, deformingBodies,
        getModeCoupledMaximumResponseDegree( loveNumbers ),
        getModeCoupledMaximumResponseOrder( loveNumbers ) ), loveNumbers_( loveNumbers )
    {
        // Set basic deformation functon as function to be evaluated when requesting variations.
        correctionFunctions.push_back(
            std::bind(
                &ModeCoupledSolidBodyTideGravityFieldVariations::addBasicSolidBodyTideCorrections,
                this, std::placeholders::_1, std::placeholders::_2 ) );
    }

    //! Destructor
    /*!
     *  Destructor
     */
    virtual ~ModeCoupledSolidBodyTideGravityFieldVariations( ){ }

    std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > getLoveNumbers( )
    {
        return loveNumbers_;
    }

    void resetLoveNumber(
        const std::pair< int, int > forcingIndices, const std::pair< int, int > responseIndices, const double loveNumber )
    {
        if( loveNumbers_.count( forcingIndices ) == 0 )
        {
            throw std::runtime_error( "Error when resetting mode-coupled Love number, no number at forcing D/O " +
                                      std::to_string( forcingIndices.first ) + "/" +
                                      std::to_string( forcingIndices.second ) + " found " );
        }
        else
        {
            if( loveNumbers_.at( forcingIndices ).count( responseIndices ) == 0  )
            {
                throw std::runtime_error( "Error when resetting mode-coupled Love number, no number at forcing D/O " +
                                          std::to_string( forcingIndices.first ) + "/" + std::to_string( forcingIndices.second ) +
                                          " and response D/O " +  std::to_string( responseIndices.first ) + "/" + std::to_string( responseIndices.second ) +" found ");
            }
        }
        loveNumbers_[ forcingIndices ][ responseIndices ] = loveNumber;
    }

protected:

    virtual void addBasicSolidBodyTideCorrections(
        Eigen::MatrixXd& cTermCorrections, Eigen::MatrixXd& sTermCorrections );

    std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > loveNumbers_;

};
} // namespace gravitation

} // namespace tudat
#endif // TUDAT_BASICSOLIDBODYTIDEGRAVITYFIELDVARIATIONS_H
