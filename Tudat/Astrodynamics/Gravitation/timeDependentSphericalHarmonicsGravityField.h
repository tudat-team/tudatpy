/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TIMEDEPENDENTSPHERICALHARMONICSGRAVITYFIELD_H
#define TUDAT_TIMEDEPENDENTSPHERICALHARMONICSGRAVITYFIELD_H

#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <vector>

#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"


namespace tudat
{

namespace gravitation
{

//! Class for time dependent spherical harmonic gravity field.
/*!
 *  Class for time dependent spherical harmonic gravity field, i.e where the sine and cosine
 *  coefficients are functions of time. This class combines nominal values with time-dependent
 *  variations, which are calculated by object of GravityFieldVariation derived classes (one object
 *  per variation). All spherical harmonic coefficients used in this class, as well as the
 *  variations, are implicitly assumed to be geodesy-normalized.
 */
class TimeDependentSphericalHarmonicsGravityField: public SphericalHarmonicsGravityField
{
public:

    //! Semi-dummy constructor, used for setting up gravity field and variations.
    /*!
     *  Semi-dummy constructor, used for setting up gravity field and variations.  This constructor
     *  is neede, as some gravity field variations need to be linked to properties of the nominal
     *  gravity field (for instance gravitational parameter), but the complete object of this type
     *  cannot be created until all variations are created, causing a circular dependency. The
     *  object is fully created when subsequently calling the setFieldVariationSettings function and
     *  setting the field variation objects.
     *  \param gravitationalParameter Gravitational parameter of massive body.
     *  \param referenceRadius Reference radius of spherical harmonic field expansion.
     *  \param nominalCosineCoefficients Nominal (i.e. with zero variation) cosine spherical
     *  harmonic coefficients.
     *  \param nominalSineCoefficients Nominal (i.e. with zero variation) sine spherical harmonic
     *  coefficients.
     *  \param fixedReferenceFrame Identifier for body-fixed reference frame to which the field is
     *  fixed (optional).
     */
    TimeDependentSphericalHarmonicsGravityField(
            const double gravitationalParameter, const double referenceRadius,
            const Eigen::MatrixXd& nominalCosineCoefficients,
            const Eigen::MatrixXd& nominalSineCoefficients,
            const std::string& fixedReferenceFrame = "" ):
        SphericalHarmonicsGravityField(
            gravitationalParameter, referenceRadius, nominalCosineCoefficients,
            nominalSineCoefficients, fixedReferenceFrame ),
        nominalSineCoefficients_( nominalSineCoefficients ),
        nominalCosineCoefficients_( nominalCosineCoefficients )
    { }

    //! Full class constructor.
    /*!
     *  Full class constructor.
     *  \param gravitationalParameter Gravitational parameter of massive body.
     *  \param referenceRadius Reference radius of spherical harmonic field expansion.
     *  \param nominalCosineCoefficients Nominal (i.e. with zero variation) cosine spherical
     *  harmonic coefficients.
     *  \param nominalSineCoefficients Nominal (i.e. with zero variation) sine spherical harmonic
     *  coefficients.
     *  \param gravityFieldVariationUpdateSettings Object containing all gravity field variations
     *  and related settings.
     *  \param fixedReferenceFrame Identifier for body-fixed reference frame to which the field is
     *  fixed (optional).
     */
    TimeDependentSphericalHarmonicsGravityField(
            const double gravitationalParameter, const double referenceRadius,
            const Eigen::MatrixXd& nominalCosineCoefficients,
            const Eigen::MatrixXd& nominalSineCoefficients,
            const boost::shared_ptr< GravityFieldVariationsSet > gravityFieldVariationUpdateSettings,
            const std::string& fixedReferenceFrame = "" ):
        SphericalHarmonicsGravityField(
            gravitationalParameter, referenceRadius,
            nominalCosineCoefficients, nominalSineCoefficients, fixedReferenceFrame ),
        nominalSineCoefficients_( nominalSineCoefficients ),
        nominalCosineCoefficients_( nominalCosineCoefficients ),
        gravityFieldVariationsSet_( gravityFieldVariationUpdateSettings )
    {
        updateCorrectionFunctions( );
    }

    //! Destructor
    /*!
     *  Destructor
     */
    ~TimeDependentSphericalHarmonicsGravityField( ){ }

    //! Update gravity field to current time.
    /*!
     *  Update gravity field coefficient corrections to current time. All correction functions are
     *  called and subsequently added to the nominal value.
     *  \param time Current time.
     */
    void update( const double time );

    //! Update correction functions.
    /*!
     *  Update correction functions, for instance to account for changed changed environmental
     *  parameters.
     */
    void updateCorrectionFunctions( )
    {
        // Check if field variation set exists.
        if( gravityFieldVariationsSet_ == NULL )
        {
            throw std::runtime_error( "Warning, gravity field coefficient update functions are NULL when requesting update" );
        }
        else
        {
            // Reset correction functions.
            correctionFunctions_ = gravityFieldVariationsSet_->getVariationFunctions( );
        }

    }

    //! Function to (re)set the gravity field variations
    /*!
     *  Function to (re)set the gravity field variations object. An option is provided for
     *  determining whether or not the variations should be immediately recalculated.
     *  \param gravityFieldVariationUpdateSettings Object storing all variation models and
     *  associated settings.
     *  \param updateCorrections Flag to determine whether the gravity field variation functions
     *  should be immediately updated with new settings.
     */
    void setFieldVariationSettings(
           const boost::shared_ptr< GravityFieldVariationsSet > gravityFieldVariationUpdateSettings,
           const bool updateCorrections = 1 );

    //! Function to clear all gravity field variations
    /*!
     *  Function to clear all gravity field variations, the gravityFieldVariationsSet_ is set to
     *  NULL, and the correctionFunctions_ list is cleared.
     */
    void clearVariations( );

    //! Get nominal (i.e. with zero variations) cosine coefficients.
    /*!
     *  Function to get nominal (i.e. with zero variations) cosine coefficients.
     *  \return Nominal cosine coefficients.
     */
    Eigen::MatrixXd getNominalCosineCoefficients( )
    {
        return nominalCosineCoefficients_;
    }

    //! Set nominal (i.e. with zero variations) cosine coefficients.
    /*!
     *  Function to set nominal (i.e. with zero variations) cosine coefficients.
     *  \param nominalCosineCoefficients New nominal cosine coefficients.
     */
    void setNominalCosineCoefficients( Eigen::MatrixXd nominalCosineCoefficients )
    {
        nominalCosineCoefficients_ = nominalCosineCoefficients;
    }

    //! Set nominal (i.e. with zero variations) cosine coefficient of given degree and order.
    /*!
     *  Set nominal (i.e. with zero variations) cosine coefficients of given degree and order.
     *  \param degree Spherical harmonic degree.
     *  \param order Spherical harmonic order.
     *  \param coefficient New cosine coefficient for given degree and order.
     */
    void setNominalCosineCoefficient( const int degree, const int order, const double coefficient )
    {
        if( degree <= nominalCosineCoefficients_.rows( ) &&
                order <= nominalCosineCoefficients_.cols( ) )
        {
            nominalCosineCoefficients_( degree, order ) = coefficient;
        }
        else
        {
            throw std::runtime_error( "Error when resetting nominal cosine coefficient" );
        }
    }

    //! Get nominal (i.e. with zero variations) sine coefficients.
    /*!
     *  Function to get nominal (i.e. with zero variations) sine coefficients.
     *  \return Nominal sine coefficients.
     */
    Eigen::MatrixXd getNominalSineCoefficients( )
    {
        return nominalSineCoefficients_;
    }

    //! Set nominal (i.e. with zero variations) sine coefficients.
    /*!
     *  Function to set nominal (i.e. with zero variations) sine coefficients.
     *  \param nominalSineCoefficients New nominal sine coefficients.
     */
    void setNominalSineCoefficients( const Eigen::MatrixXd& nominalSineCoefficients )
    {
        nominalSineCoefficients_ = nominalSineCoefficients;
    }

    //! Set nominal (i.e. with zero variations) sine coefficient of given degree and order.
    /*!
     *  Set nominal (i.e. with zero variations) sine coefficients of given degree and order.
     *  \param degree Spherical harmonic degree.
     *  \param order Spherical harmonic order.
     *  \param coefficient New sine coefficient for given degree and order.
     */
    void setNominalSineCoefficient( const int degree, const int order, const double coefficient )
    {
        if( degree <= nominalSineCoefficients_.rows( ) &&
                order <= nominalSineCoefficients_.cols( ) )
        {
            nominalSineCoefficients_( degree, order ) = coefficient;
        }
        else
        {
            throw std::runtime_error( "Error when resetting nominal sine coefficient" );
        }
    }

    //! Function to get object containing all gravity field variations and related settings
    /*!
     *  Function to get object containing all gravity field variations and related settings
     */
    boost::shared_ptr< GravityFieldVariationsSet > getGravityFieldVariationsSet( )
    {
        return gravityFieldVariationsSet_;
    }

private:

    //! Nominal (i.e. with zero variations) cosine coefficients.
    /*!
     *  Nominal (i.e. with zero variations) cosine coefficients. When calling the update function,
     *  all corrections are calculated and the sum of these corrections and this nominal value is
     *  set as cosineCoefficients_ base class member.
     */
    Eigen::MatrixXd nominalSineCoefficients_;

    //! Nominal (i.e. with zero variations) sine coefficients.
    /*!
     *  Nominal (i.e. with zero variations) sine coefficients. When calling the update function,
     *  all corrections are calculated and the sum of these corrections and this nominal value is
     *  set as sineCoefficients_ base class member.
     */
    Eigen::MatrixXd nominalCosineCoefficients_;

    //! List of update functions which are called when calculating current gravity field variations.
    /*!
     *  List of update functions which are called when calculating current gravity field variations.
     *  Functions are either linked to PairInterpolationInterface, which contains an interpolator
     *  created from a GravityFieldVariations object, or the addSphericalHarmonicsCorrections of the
     *  GravityFieldVariations object directly.
     */
    std::vector< boost::function< void( const double, Eigen::MatrixXd&, Eigen::MatrixXd& ) > >
        correctionFunctions_;

    //! Object containing all GravityFieldVariations objects and update settings.
    /*!
     *  Object containing all GravityFieldVariations objects and update settings
     *  (i.e. time settings for interpolator)
     */
    boost::shared_ptr< GravityFieldVariationsSet > gravityFieldVariationsSet_;

};

} // namespace gravitation

} // namespace tudat

#endif // TUDAT_TIMEDEPENDENTSPHERICALHARMONICSGRAVITYFIELD_H
