/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GRAVITYFIELDVARIATIONS_H
#define TUDAT_GRAVITYFIELDVARIATIONS_H

#include <functional>
#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/interpolators/createInterpolator.h"

namespace tudat
{

namespace gravitation
{

enum BodyDeformationTypes
{
    basic_solid_body,
    mode_coupled_solid_body,
    tabulated_variation,
    periodic_variation,
    polynomial_variation,
    iers_2010,
    pole_tide,
    ocean_tide
};


//! Interface class between GravityFieldVariations objects that are interpolated and
//! TimeDependentSphericalHarmonicsGravityField.
/*!
 *  Interface class between GravityFieldVariations objects that are interpolated and
 *  TimeDependentSphericalHarmonicsGravityField, the getCosineSinePair function mimics the
 *  addSphericalHarmonicsCorrections function of GravityFieldVariations.
 *  All correction coefficients are calculated as rectnagular blocks in both the cosine and sine
 *  matrices (of equal size)
 */
class PairInterpolationInterface
{
public:

    //! Class constructor
    /*!
     *  Constructor, receives an interpolator (typically created from a GravityFieldVariations
     *  object by the createInterpolatedSphericalHarmonicCorrectionFunctions function), which is
     *  used to approximate the spherical harmonic coefficient corrections at any given time
     * (inside the interpolation window). All correction coefficients are calculated as rectangular
     *  blocks in both the cosine and sine matrices (of equal size)
     *  \param cosineSineInterpolator Interpolator object for approximating coefficient corrections.
     *  \param startDegree Degree where the rectangular correction block starts.
     *  \param startOrder Order where the rectangular correction block starts.
     *  \param numberOfDegrees Size of the rectangular correction block in the degree direction.
     *  \param numberOfOrders Size of the rectangular correction block in the order direction.
     */
    PairInterpolationInterface(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator<
            double, Eigen::MatrixXd > > cosineSineInterpolator,
            const int startDegree, const int startOrder,
            const int numberOfDegrees, const int numberOfOrders ):
        cosineSineInterpolator_( cosineSineInterpolator ),
        startDegree_( startDegree ), startOrder_( startOrder ),
        numberOfDegrees_( numberOfDegrees ), numberOfOrders_( numberOfOrders )
    { }

    //! Function to add sine and cosine corrections at given time to coefficient matrices.
    /*!
     *  Function to add sine and cosine corrections at given time to coefficient matrices.
     *  The current sine and cosine matrices are passed by reference, the corrections are
     *  calculated internally and added to them.
     *  \param time Time at which corrections are to be evaluated.
     *  \param sineCoefficients Current spherical harmonic sine coefficients, calculated
     *  corrections are added and returned by reference
     *  \param cosineCoefficients Current spherical harmonic cosine coefficients, calculated
     *  corrections are added and returned by reference
     */
    void getCosineSinePair( const double time,
                            Eigen::MatrixXd& sineCoefficients,
                            Eigen::MatrixXd& cosineCoefficients );

private:

    //! Interpolator object for approximating coefficient corrections.
    /*!
     *  Interpolator object for approximating coefficient corrections, interpolates the cosine and
     *  sine corrections, concatenated next to each other (i.e. as [C S] row 'vector' of
     *  cosine and sine corrections C and S, respectively).
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
    cosineSineInterpolator_;

    //! Degree where the rectangular correction block starts.
    /*!
     *  Degree where the rectangular correction block starts.
     */
    int startDegree_;

    //! Order where the rectangular correction block starts.
    /*!
     *  Order where the rectangular correction block starts.
     */
    int startOrder_;

    //! Size of the rectangular correction block in the degree direction.
    /*!
     *  Size of the rectangular correction block in the degree direction.
     */
    int numberOfDegrees_;

    //! Size of the rectangular correction block in the order direction.
    /*!
     *  Size of the rectangular correction block in the order direction.
     */
    int numberOfOrders_;
};

//! Virual base class for spherical harmonic gravity field variations
//! (i.e. time-dependencies of sine and cosine coefficients)
/*!
 *  Virual base class for spherical harmonic gravity field variations
 *  (i.e. time-dependencies of sine and cosine coefficients). The interface of derived classes
 *  with TimeDependentSphericalHarmonicsGravityField, which combines all corrections and adds them
 *  to nominal values, can be performed either directly, or through an interpolator
 *  (see PairInterpolationInterface and createInterpolatedSphericalHarmonicCorrectionFunctions)
 *  to save on computation time when evaluating slowly changing functions at very short intervals.
 *  All correction coefficients are calculated as rectangular blocks in both the cosine and sine
 *  matrices (of equal size).
 */
class GravityFieldVariations
{
public:

    //! Base class constructor
    /*!
     *  Base class constructor, input defines the size and position of correction blocks in sine and
     *  cosine matrices
     *  \param minimumDegree Degree where the rectangular correction blocks start.
     *  \param minimumOrder Order where the rectangular correction blocks start.
     *  \param maximumDegree Degree where the rectangular correction blocks end.
     *  \param maximumOrder Order where the rectangular correction blocks end.
     */
    GravityFieldVariations( const int minimumDegree, const int minimumOrder,
                            const int maximumDegree, const int maximumOrder ):
        minimumDegree_( minimumDegree ), minimumOrder_( minimumOrder ),

                maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder )
    {
        numberOfDegrees_ = maximumDegree_ - minimumDegree_ + 1;
        numberOfOrders_ = maximumOrder_ - minimumOrder_ + 1;
        lastCosineCorrection_.setZero( maximumDegree_ + 1, maximumOrder_ + 1 );
        lastSineCorrection_.setZero( maximumDegree_ + 1, maximumOrder_ + 1 );
    }

    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    virtual ~GravityFieldVariations( ){ }

    //! Pure virtual function for calculating corrections.
    /*!
     *  Pure virtual function for calculating corrections at given time.
     *  \param time Time at which variations are to be calculated.
     *  \return Pair of  matrices containing variations in (cosine, sine) coefficients at block
     *  positions in total matrices defined by minimumDegree_, minimumOrder_, numberOfDegrees_,
     *  numberOfOrders_;
     */
    virtual std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections(
            const double time ) = 0;

    //! Function to add sine and cosine corrections at given time to coefficient matrices.
    /*!
     *  Function to add sine and cosine corrections at given time to coefficient matrices.
     *  The current sine and cosine matrices are passed by reference, the corrections are calculated
     *  internally and added to them.
     *  \param time Time at which corrections are to be evaluated.
     *  \param sineCoefficients Current spherical harmonic sine coefficients, calculated
     *  corrections are added and returned by reference
     *  \param cosineCoefficients Current spherical harmonic cosine coefficients, calculated
     *  corrections are added and returned by reference
     */
    void addSphericalHarmonicsCorrections(
            const double time,
            Eigen::MatrixXd& sineCoefficients,
            Eigen::MatrixXd& cosineCoefficients );

    //! Function to return the maximum degree of the corrections.
    /*!
     *  Function to return the maximum degree of the corrections.
     */
    int getMaximumDegree( )
    {
        return maximumDegree_;
    }

    //! Function to return the maximum order of the corrections.
    /*!
     *  Function to return the maximum order of the corrections.
     */
    int getMaximumOrder( )
    {
        return maximumOrder_;
    }

    //! Function to return the minimum degree of the corrections.
    /*!
     *  Function to return the minimum degree of the corrections.
     */
    int getMinimumDegree( )
    {
        return minimumDegree_;
    }

    //! Function to return the minimum order of the corrections.
    /*!
     *  Function to return the minimum order of the corrections.
     */
    int getMinimumOrder( )
    {
        return minimumOrder_;
    }

    //! Function to return the number of degrees of the corrections.
    /*!
     *  Function to return the number of degrees of the corrections.
     */
    int getNumberOfDegrees( )
    {
        return numberOfDegrees_;
    }

    //! Function to return the number of orders of the corrections.
    /*!
     *  Function to return the number of orders of the corrections.
     */
    int getNumberOfOrders( )
    {
        return numberOfOrders_;
    }

    //! Function to retrieve correction to cosine coefficients, as computed by last call to addSphericalHarmonicsCorrections
    /*!
     *  Function to retrieve correction to cosine coefficients, as computed by last call to addSphericalHarmonicsCorrections
     * \return Correction to cosine coefficients, as computed by last call to addSphericalHarmonicsCorrections
     */
    Eigen::MatrixXd getLastCosineCorrection( )
    {
        return lastCosineCorrection_;
    }

    //! Function to retrieve correction to sine coefficients, as computed by last call to addSphericalHarmonicsCorrections
    /*!
     *  Function to retrieve correction tosine coefficients, as computed by last call to addSphericalHarmonicsCorrections
     * \return Correction to sine coefficients, as computed by last call to addSphericalHarmonicsCorrections
     */
    Eigen::MatrixXd getLastSineCorrection( )
    {
        return lastSineCorrection_;
    }

protected:

    //! Minimum degree of variations
    /*!
     *  Minimum degree of variations
     */
    int minimumDegree_;

    //! Minimum order of variations
    /*!
     *  Minimum order of variations
     */
    int minimumOrder_;

    //! Maximum degree of variations
    /*!
     *  Maximum degree of variations
     */
    int maximumDegree_;

    //! Maximum order of variations
    /*!
     *  Maximum order of variations
     */
    int maximumOrder_;

    //! Number of degrees of variations
    /*!
     *  Number of degrees of variations
     */
    int numberOfDegrees_;

    //! Number of orders of variations
    /*!
     *  Number of orders of variations
     */
    int numberOfOrders_;    

    //! Latest correction to cosine coefficients, as computed by last call to addSphericalHarmonicsCorrections
    Eigen::MatrixXd lastCosineCorrection_;

    //! Latest correction to sine coefficients, as computed by last call to addSphericalHarmonicsCorrections
    Eigen::MatrixXd lastSineCorrection_;
};

//! Function to create a function linearly interpolating the sine and cosine correction coefficients
//! produced by an object of GravityFieldVariations type.
/*!
 *  Function to create a function linearly interpolating the sine and cosine correction coefficients
 *  produced by an object of GravityFieldVariations type.
 *  The function creates a function pointer to the getCosineSinePair function of
 *  PairInterpolationInterface. The object of type PairInterpolationInterface is created by
 *  generating an interpolator for sine/cosine coefficients from the variationObject object and
 *  the initial/final time and time step that are passed.
 *  \param variationObject Object generating cosine and sine coefficient corrections.
 *  \param initialTime Start time of interpolator.
 *  \param finalTime End time of interpolator.
 *  \param timeStep Time step between subsequent evaluations of coefficient corrections
 *  \return Function pointer to function mimicing the addSphericalHarmonicsCorrections
 *  function of GravityFieldVariations.
 */
std::function< void( const double, Eigen::MatrixXd&, Eigen::MatrixXd& ) >
createInterpolatedSphericalHarmonicCorrectionFunctions(
        std::shared_ptr< GravityFieldVariations > variationObject,
        const double initialTime,
        const double finalTime,
        const double timeStep,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        std::make_shared< interpolators::InterpolatorSettings >(
            interpolators::linear_interpolator, interpolators::huntingAlgorithm ) );

//! Container class containing all gravity field variations for a single Body
//! (and TimeDependentSphericalHarmonicsGravityField).
/*!
 *  Container class containing all gravity field variations for a single Body
 *  (and TimeDependentSphericalHarmonicsGravityField). Also contains information on whether an
 *  interpolator is used for calculation of corrections by
 *  TimeDependentSphericalHarmonicsGravityField and, if so, the associated interpolator settings.
 *  Multiple corrections of a single type may be contained in this class, in which case each of them
 *  must be supplied with a unique identifier (string).
 */
class GravityFieldVariationsSet
{
public:

    //! Class constructor.
    /*!
     *  Class contructor, requires set of correction objects (and associated properties).
     *  \param variationObjects List of GravityFieldVariations objects denoting the complete set of
     *  variations to take into account.
     *  \param variationType List of type identifiers of variationObjects (prevents use of dynamic
     *  casts), must be of same size as variationObjects.
     *  \param variationIdentifier Name of variation object for each entry of variationObjects, must
     *  be of same size as variationObjects.
     *  Entries only required to be non-empty if multiple variations objects of same type are
     *  included in variationObjects list.
     *  \param createInterpolator List of booleans denoting whether to interpolate a given
     *  entry of variationObjects or not, must be of same size as variationObjects.
     *  \param initialTimes Initial times for interpolation, must contain an entry for each
     *  variation (map key denotes index of variationObjects) for which createInterpolator is true.
     *  \param finalTimes Final times for interpolation, must contain an entry for each variation
     *  (map key denotes index of variationObjects) for which createInterpolator is true.
     *  \param timeSteps Time steps for interpolation, must contain an entry for each variation
     *  (map key denotes index of variationObjects) for which createInterpolator is true.
     */
    GravityFieldVariationsSet(
            const std::vector< std::shared_ptr< GravityFieldVariations > > variationObjects,
            const std::vector< BodyDeformationTypes > variationType,
            const std::vector< std::string > variationIdentifier,
            const std::map< int, std::shared_ptr< interpolators::InterpolatorSettings > >
            createInterpolator =
            std::map< int, std::shared_ptr< interpolators::InterpolatorSettings > >( ),
            const std::map< int, double > initialTimes = std::map< int, double >( ),
            const std::map< int, double > finalTimes = std::map< int, double >( ),
            const std::map< int, double > timeSteps = std::map< int, double >( ) );

    //! Function to retrieve a variation object of given type (and name if necessary).
    /*!
     *  Function to retrieve a variation object of given type (and name if necessary).
     *  Name must be provided only if if multiple variation objects of same type are included in
     *  variationObjects_ list.
     *  \param deformationType Type of gravity field variation.
     *  \param identifier Name of gravity field variation (only required if if multiple variation
     *  objects of same type are included in variationObjects_ list (ignored otherwise).
     *  \return Pair containing boolean (true if requested variation found, false otherwise) and
     *  pointer to variation object (only if requested variation found).
     */
    std::pair< bool, std::shared_ptr< gravitation::GravityFieldVariations > >
     getGravityFieldVariation(
            const BodyDeformationTypes deformationType,
            const std::string identifier = "" );

    //! Function to retrieve list of variation functions.
    /*!
     *  Function to retrieve list of variation functions, entries are either created using function
     *  pointer binding to PairInterpolationInterface (if given variation is to be interpolated)
     *  or to GravityFieldVariations directly (if no interpolation requested).
     *  \return List of gravity field coefficient variation functions, matching the interface of
     *  GravityFieldVariations::addSphericalHarmonicsCorrections
     */
    std::vector< std::function< void( const double, Eigen::MatrixXd&, Eigen::MatrixXd& ) > >
    getVariationFunctions( );

    //! Function to retrieve the complete set of variations to take nto account.
    /*!
     * Function to retrieve the complete set of variations to take nto account.
     * \return Complete set of variations to take nto account.
     */
    std::vector< std::shared_ptr< GravityFieldVariations > > getVariationObjects( )
    {
        return variationObjects_;
    }

    //! Function to retrieve the tidal gravity field variation with the specified bodies causing deformation
    /*!
     * Function to retrieve the tidal gravity field variation with the specified bodies causing deformation. If the
     * deformingBodies list is empty, and only one tidal gravity field variation exists, this object is returned. Function
     * throws an exception in no object correspond to input is found
     * \param deformingBodies List of objects that cause tidal gravity field variation
     * \return The tidal gravity field variation with the specified bodies causing deformation
     */
    std::shared_ptr< GravityFieldVariations > getDirectTidalGravityFieldVariation(
            const std::vector< std::string >& deformingBodies,
            const BodyDeformationTypes tideType = basic_solid_body );

    //! Function to retrieve the tidal gravity field variations
    /*!
     * Function to retrieve the tidal gravity field variations
     * \return List of tidal gravity field variations objects
     */
    std::vector< std::shared_ptr< GravityFieldVariations > > getDirectTidalGravityFieldVariations( );

private:

    //! List of GravityFieldVariations objects denoting the complete set of variations to take nto account.
    /*!
     *  List of GravityFieldVariations objects denoting the complete set of variations to take into account.
     */
    std::vector< std::shared_ptr< GravityFieldVariations > > variationObjects_;

    //! List of type identifiers of variationObjects.
    /*!
     *  List of type identifiers of variationObjects (prevents use of dynamic casts), must be of
     *  same size as variationObjects.
     */
    std::vector< BodyDeformationTypes > variationType_;

    //! List of model identifiers per deformation type.
    std::map< BodyDeformationTypes, std::vector< std::string > > identifierPerType_;

    //! Name of variation object for each entry of variationObjects
    /*!
     *  Name of variation object for each entry of variationObjects, must be of same size as
     *  variationObjects. Used for discriminating between different variation objects of same type.
     *  Entries only required to be non-empty if multiple variations objects of same type are
     *  included in variationObjects list.
     */
    std::vector< std::string > variationIdentifier_;

    //! List of booleans denoting whether to interpolate a given entry of variationObjects or not
    /*!
     *  List of booleans denoting whether to interpolate a given entry of variationObjects or not,
     *  must be of same size as variationObjects.
     */
    std::map< int, std::shared_ptr< interpolators::InterpolatorSettings > > createInterpolator_;

    //! Initial times for interpolation,
    /*!
     *  Initial times for interpolation, must contain an entry for each variation
     *  (map key denotes index of variationObjects) for which createInterpolator is true.
     */
    std::map< int, double > initialTimes_;

    //! Final times for interpolation,
    /*!
     *  Final times for interpolation, must contain an entry for each variation
     *  (map key denotes index of variationObjects) for which createInterpolator is true.
     */
    std::map< int, double > finalTimes_;

    //! Time steps for interpolation,
    /*!
     *  Time steps for interpolation, must contain an entry for each variation
     *  (map key denotes index of variationObjects) for which createInterpolator is true.
     */
    std::map< int, double > timeSteps_;
};

} // namespace gravitation

} // namespace tudat

#endif // TUDAT_GRAVITYFIELDVARIATIONS_H
