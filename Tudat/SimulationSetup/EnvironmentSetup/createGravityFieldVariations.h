/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <boost/assign/list_of.hpp>


#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/Mathematics/Interpolators/createInterpolator.h"

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
            const boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            boost::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ) ):
        interpolatorSettings_( interpolatorSettings ), initialTime_( initialTime ),
        finalTime_( finalTime ), timeStep_( timeStep ){ }

    boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings_;

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
     * is NULL, in which no interpolation is used, and the model is evaluated during propagation.
     */
    GravityFieldVariationSettings( const gravitation::BodyDeformationTypes bodyDeformationType,
                                   const boost::shared_ptr< ModelInterpolationSettings > interpolatorSettings = NULL ):
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
     * variations (if not NULL)
     */
    boost::shared_ptr< ModelInterpolationSettings > getInterpolatorSettings( ){ return interpolatorSettings_; }

protected:

    //! Type of gravity field variation to be used.
    gravitation::BodyDeformationTypes bodyDeformationType_;

    //! Settings that are to be used to create an interpolator for the gravity field variations
    /*!
     * Settings that are to be used to create an interpolator for the gravity field variations
     * immediately upon creation. If NULL, no interpolation is used, and the model is evaluated
     * during propagation.
     */
    boost::shared_ptr< ModelInterpolationSettings > interpolatorSettings_;

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
     * \param bodyReferenceRadius Reference (typically equatorial) radius of body being deformed
     * \param interpolatorSettings Settings that are to be used to create an interpolator for the
     * gravity field variations immediately upon creation (to be used during propagation). Default
     * is NULL, in which no interpolation is used, and the model is evaluated during propagation.
     */
    BasicSolidBodyGravityFieldVariationSettings(
            const std::vector< std::string > deformingBodies,
            const std::vector< std::vector< std::complex< double > > > loveNumbers,
            const double bodyReferenceRadius,
            const boost::shared_ptr< ModelInterpolationSettings > interpolatorSettings = NULL ):
        GravityFieldVariationSettings( gravitation::basic_solid_body, interpolatorSettings ),
        deformingBodies_( deformingBodies ), loveNumbers_( loveNumbers ),
                bodyReferenceRadius_( bodyReferenceRadius ){ }

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
    std::vector< std::vector< std::complex< double > > > getLoveNumbers( ){ return loveNumbers_; }

    //! Function to retrieve reference (typically equatorial) radius of body being deformed
    /*!
     * \brief Function to retrieve reference (typically equatorial) radius of body being deformed
     * \return Reference (typically equatorial) radius of body being deformed
     */
    double getBodyReferenceRadius( ){ return bodyReferenceRadius_; }

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
    std::vector< std::vector< std::complex< double > > > loveNumbers_;

    //! Reference (typically equatorial) radius of body being deformed
    double bodyReferenceRadius_;

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
            const boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings ):
        GravityFieldVariationSettings(
            gravitation::tabulated_variation, boost::make_shared< ModelInterpolationSettings >(
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

//! Function to create a set of gravity field variations, stored in the associated interface class
/*!
 * Function to create a set of gravity field variations, stored in the associated interface class of
 * type GravityFieldVariationsSet
 * \param body Body for which gravity field variations are createad
 * \param bodyMap List of body objects in simulations.
 * \param gravityFieldVariationSettings List of settings for gravity field variations
 * \return Interface class containing list of GravityFieldVariations.
 */
boost::shared_ptr< gravitation::GravityFieldVariationsSet > createGravityFieldModelVariationsSet(
        const std::string& body,
        const NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< GravityFieldVariationSettings > >&
            gravityFieldVariationSettings );

//! Function to create a single gravity field variation object.
/*!
 * Function to create a single gravity field variation object.
 * \param gravityFieldVariationSettings Settings for the gravity field variation object that is to
 * be created. Depending on the settings and type of the gravity field variations, spherical
 * harmonic corrections are calculated a priori and handled by an interpolator during propagation,
 * or they are directly calculated from teh current state during numerical propagation.
 * \param body Body for which gravity field variations are createad
 * \param bodyMap List of body objects in simulations.
 * \return Single gravity field variation object.
 */
boost::shared_ptr< gravitation::GravityFieldVariations > createGravityFieldVariationsModel(
        const boost::shared_ptr< GravityFieldVariationSettings > gravityFieldVariationSettings,
        const std::string body,
        const NamedBodyMap& bodyMap );

} // namespace simulation_setup
} // namespace tudat
#endif // TUDAT_CREATEGRAVITYFIELDVARIATIONS_H
