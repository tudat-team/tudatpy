/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TIDALLOVENUMBERPARTIALINTERFACE_H
#define TUDAT_TIDALLOVENUMBERPARTIALINTERFACE_H

#include <boost/math/special_functions/factorials.hpp>

#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"

#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicPartialFunctions.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/tidalLoveNumber.h"

namespace tudat
{

namespace orbit_determination
{

void getMaximumUsedDegreeAndOrder(
        const int maximumDegree, const int maximumOrder, const int evaluationDegree,
        int& maximumUsedDegree, int& maximumUsedOrder );

//! (Base) class for calculating the partial of a spherical harmonic acceleration w.r.t. a tidal love number.
/*!
 *  (Base) class for calculating the partial of a spherical harmonic acceleration w.r.t. a tidal love number. These calculations
 *  are implemented separately from the SphericalHarmonicsGravityPartial class as its functionality is quite specific and only
 *  used by the SphericalHarmonicsGravityPartial for certain cases. Also, this architecture allows different implementations for
 *  different tidal models. This base class assumes no direct tidal lag in the calculations, i.e. state and rotation functions are
 *  all evaluated at same time and lag is implemented by the possibility of complex love numbers (for non-zonal terms).
 */
class TidalLoveNumberPartialInterface
{
public:

    //! Constructor from objects.
    /*!
     *  Constructor from objects.
     *  \param gravityFieldVariations Gravity field variation object which calculates spherical harmonic coefficient variations
     *  due to considered love number(s).
     *  \param deformedBodyPositionFunction Function to retrieve current state of body being tidally deformed
     *  \param deformingBodyStateFunctions Functions for retrieving states at current time of bodies causing deformation.
     *  \param rotationToDeformedBodyFrameFrameFunction Function returning the rotation from inertial to body-fixed frame of
     *  deformed body.
     *  \param deformedBody Name of body being tidally deformed.
     */
    TidalLoveNumberPartialInterface(
            const std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariations,
            const std::function< Eigen::Vector3d( ) > deformedBodyPositionFunction,
            const std::vector< std::function< Eigen::Vector3d( ) > > deformingBodyStateFunctions,
            const std::function< Eigen::Quaterniond( ) > rotationToDeformedBodyFrameFrameFunction,
            const std::string& deformedBody ):
        deformedBodyPositionFunction_( deformedBodyPositionFunction ),
        deformingBodyStateFunctions_( deformingBodyStateFunctions ),
        rotationToDeformedBodyFrameFrameFunction_( rotationToDeformedBodyFrameFrameFunction ),
        deformedBody_( deformedBody )
    {
        // Get members from input objects.
        deformingBodyGravitationalParameters_ = gravityFieldVariations->getDeformingBodyMasses( );
        positionsOfDeformingBodies_.resize( deformingBodyStateFunctions_.size( ) );
        deformedBodyReferenceRadius_ = gravityFieldVariations->getDeformedBodyReferenceRadius( );
        deformedBodyGravitationalParameterFunction_ = gravityFieldVariations->getDeformedBodyMassFunction( );
        deformingBodies_  = gravityFieldVariations->getDeformingBodies( );


        realLoveNumberScaler_ =
                std::make_pair( ( Eigen::Vector2d( ) << 1.0, 0.0 ).finished( ), ( Eigen::Vector2d( ) << 0.0, 1.0 ).finished( ) );
        complexLoveNumberScaler_ =
                std::make_pair( ( Eigen::Matrix2d( ) << 1.0, 0.0, 0.0, -1.0  ).finished( ),
                                                   ( Eigen::Matrix2d( ) << 0.0, 1.0, 1.0, 0.0 ).finished( ) );
        for( unsigned int i = 0; i < deformingBodyStateFunctions_.size( ); i++ )
        {
            allDeformingBodyIndices_.push_back( i );
        }

    }

    //! Destructor
    virtual ~TidalLoveNumberPartialInterface( ){ }

    //! Function to obtain the indices of given list of body names in deformingBodies_ member vector
    /*!
     * Function to obtain the indices of given list of body names in deformingBodies_ member vector. Function throws an
     * exception if the input vector contains names of bodies not contained in the deformingBodies_ vector
     * \param selectedBodyNames List of bodies for which the indices in the deformingBodies_ vector are to be obtained
     * \return The indices of given list of body names in deformingBodies_ member vector
     */
    std::vector< int > getSelectedDeformingBodyIds( const std::vector< std::string >& selectedBodyNames );

    //! Function to calculate the partial of spherical harmonic acceleration w.r.t. complex tidal love numbers.
    /*!
     *  Function to calculate the partial of spherical harmonic acceleration w.r.t. complex tidal love numbers at a single degree.
     *  The orders of the love numbers w.r.t. which the partials to be taken are required as input.
     *  \param degree Degree of love numbers w.r.t. which partials are to be taken
     *  \param orders Spherical harmonic orders in current degree at which the partials are to be taken.
     *  \param deformingBodyIndices Indices of deforming bodies (in deformingBodies_ member vector) due to which the deformation
     *  is to be taken into account
     *  \param maximumDegree Maximum degree of the acceleration for which the partial derivatives are to be calculated.
     *  \param maximumOrder Maximum order of the acceleration for which the partial derivatives are to be calculated.
     *  \return Partial derivatives of spherical harmonic coefficients. Stl vector entries denote values at requested orders,
     *  Eigen vector rows denote C and S coefficient, respectively; columns denote real and complex tidal Love number components.
     */
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >
    calculateSphericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumbers(
            const int degree,
            const std::vector< int >& orders,
            const std::vector< int >& deformingBodyIndices,
            const int maximumDegree,
            const int maximumOrder );


    //! Function to calculate the partial of spherical harmonic acceleration w.r.t. complex tidal love numbers.
    /*!
     *  Function to calculate the partial of spherical harmonic acceleration w.r.t. complex tidal love numbers at all orders of
     *  a single degree.
     *  \param degree Degree of love numbers w.r.t. which partials are to be taken
     *  \param deformingBodyIndices Indices of deforming bodies (in deformingBodies_ member vector) due to which the deformation
     *  is to be taken into account
     *  \param maximumDegree Maximum degree of the acceleration for which the partial derivatives are to be calculated.
     *  \param maximumOrder Maximum order of the acceleration for which the partial derivatives are to be calculated.
     *  \return Partial derivatives of spherical harmonic coefficients. Stl vector entries denote values at all orders at degree,
     *  Eigen vector rows denote C and S coefficient, respectively; columns denote real and complex tidal Love number components.
     */
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >
    calculateSphericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumber(
            const int degree,
            const std::vector< int >& deformingBodyIndices,
            const int maximumDegree,
            const int maximumOrder )
    {
        return calculateSphericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumbers(
                        degree, estimatable_parameters::fullDegreeOrders[ degree - 2 ], deformingBodyIndices,
                maximumDegree, maximumOrder );
    }

    //! Function to calculate the partial of spherical harmonic acceleration w.r.t. real tidal love numbers.
    /*!
     *  Function to calculate the partial of spherical harmonic acceleration w.r.t. real tidal love numbers at a single degree.
     *  The orders of the love numbers w.r.t. which the partials to be taken are required as input.
     *  \param degree Degree of love numbers w.r.t. which partials are to be taken
     *  \param orders Spherical harmonic orders in current degree at which the partials are to be taken.
     *  \param deformingBodyIndices Indices of deforming bodies (in deformingBodies_ member vector) due to which the deformation
     *  is to be taken into account
     *  \param maximumDegree Maximum degree of the acceleration for which the partial derivatives are to be calculated.
     *  \param maximumOrder Maximum order of the acceleration for which the partial derivatives are to be calculated.
     *  \return Partial derivatives of spherical harmonic coefficients. Stl vector entries denote values at requested orders,
     *  Eigen vector rows denote C and S coefficient, respectively; Eigen vector has one column, denoting real part of tidal
     *  Love number.
     */
    virtual std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >
    calculateSphericalHarmonicCoefficientsPartialWrtRealTidalLoveNumbers(
            const int degree,
            const std::vector< int >& orders,
            const std::vector< int >& deformingBodyIndices,
            const int maximumDegree,
            const int maximumOrder );

    //! Function to calculate the partial of spherical harmonic acceleration w.r.t. real tidal love numbers.
    /*!
     *  Function to calculate the partial of spherical harmonic acceleration w.r.t. real tidal love numbers at all orders of
     *  a single degree.
     *  \param degree Degree of love numbers w.r.t. which partials are to be taken
     *  \param deformingBodyIndices Indices of deforming bodies (in deformingBodies_ member vector) due to which the deformation
     *  is to be taken into account
     *  \param maximumDegree Maximum degree of the acceleration for which the partial derivatives are to be calculated.
     *  \param maximumOrder Maximum order of the acceleration for which the partial derivatives are to be calculated.
     *  \return Partial derivatives of spherical harmonic coefficients. Stl vector entries denote values at all orders at degree,
     *  Eigen vector rows denote C and S coefficient, respectively; Eigen vector has one column, denoting real part of tidal
     *  Love number.
     */
    virtual std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >
    calculateSphericalHarmonicCoefficientsPartialWrtRealTidalLoveNumber(
            const int degree,
            const std::vector< int >& deformingBodyIndices,
            const int maximumDegree,
            const int maximumOrder )
    {
        return calculateSphericalHarmonicCoefficientsPartialWrtRealTidalLoveNumbers(
                    degree, estimatable_parameters::fullDegreeOrders[ degree - 2 ], deformingBodyIndices,
                maximumDegree, maximumOrder );
    }

    //! Update tidal Love number interface to current time
    /*!
     *  Update tidal Love number interface to current time. NOTE: This function currently has very limited functionality,
     *  resulting in reduced computational efficiency of the class.
     *  \param currentTime Current time to which interface is to be set.
     */
    void update( const double currentTime )
    {
        currentTime_ = currentTime;
        rotationToTidallyDeformedBody_ = rotationToDeformedBodyFrameFrameFunction_( );
    }

    //! Reset current time of tidal Love number interface.
    /*!
     *  Reset current time of tidal Love number interface, typically used to set NaN time, signalling that an update of the object
     *  is required
     *  \param currentTime Current time to which interface is to be set.
     */
    void resetCurrentTime( )
    {
        currentDoubleParameterPartials_.clear( );
        currentVectorParameterPartials_.clear( );
        currentTime_ = TUDAT_NAN;
    }

    //! Function to update the values of the partial derivatives to current state and time.
    /*!
     *  Function to update the values of the partial derivatives to current state and time. The list of functions computing the
     *  partial derivatives is determined by the calls to the setParameterPartialFunction functions of this class.
     */
    void updateParameterPartials( )
    {
        // Update double parameters
        for( parameterDoublePartialFunctionIterator_ = parameterDoublePartialFunctions_.begin( );
             parameterDoublePartialFunctionIterator_ != parameterDoublePartialFunctions_.end( );
             parameterDoublePartialFunctionIterator_++ )
        {
            if( currentDoubleParameterPartials_.count( parameterDoublePartialFunctionIterator_->first ) == 0 )
            {
                currentDoubleParameterPartials_[ parameterDoublePartialFunctionIterator_->first ] =
                        parameterDoublePartialFunctionIterator_->second( );
            }
        }

        // Update vector parameters
        for( parameterVectorPartialFunctionIterator_ = parameterVectorPartialFunctions_.begin( );
             parameterVectorPartialFunctionIterator_ != parameterVectorPartialFunctions_.end( );
             parameterVectorPartialFunctionIterator_++ )
        {
            if( currentVectorParameterPartials_.count( parameterVectorPartialFunctionIterator_->first ) == 0 )
            {
                currentVectorParameterPartials_[ parameterVectorPartialFunctionIterator_->first ] =
                        parameterVectorPartialFunctionIterator_->second( );
            }
        }
    }

    //! Function to set a dependency of this partial object w.r.t. a given double parameter.
    /*!
     *  Function to set a dependency of this partial object w.r.t. a given double parameter. If a dependency exists, the given
     *  partial is recomputed on every call of updateParameterPartials.
     *  \param parameter Partial w.r.t. which dependency is to be checked and set.
     *  \param maximumDegree Maximum degree of the acceleration for which the partial derivatives are to be calculated.
     *  \param maximumOrder Maximum order of the acceleration for which the partial derivatives are to be calculated.
     *  \return Pair, with first: Size (number of columns) of parameter partial. Zero if no dependency, 1 otherwise. Pair second:
     *  maximum degree and order of gravity field that are affected by the partials (NaN if no dependency).
     */
    virtual std::pair< int, std::pair< int, int > > setParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
            const int maximumDegree,
            const int maximumOrder );

    //! Function to set a dependency of this partial object w.r.t. a given vector parameter.
    /*!
     *  Function to set a dependency of this partial object w.r.t. a given vector parameter. If a dependency exists, the given
     *  partial is recomputed on every call of updateParameterPartials.
     *  \param parameter Partial w.r.t. which dependency is to be checked and set.
     *  \param maximumDegree Maximum degree of the acceleration for which the partial derivatives are to be calculated.
     *  \param maximumOrder Maximum order of the acceleration for which the partial derivatives are to be calculated.
     *  \return Pair, with first: Size (number of columns) of parameter partial. Zero if no dependency. Pair second:
     *  maximum degree and order of gravity field that are affected by the partials (NaN if no dependency).
     */
    virtual std::pair< int, std::pair< int, int > > setParameterPartialFunction(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
            const int maximumDegree,
            const int maximumOrder );

    //! Function to retrieve a partial w.r.t. a double parameter
    /*!
     * Function to retrieve a partial w.r.t. a double parameter. An error is thrown if there is no dependency w.r.t.
     * the requested parameter. A warning is printed if the dependency exists, but has not yet been computed for the
     * current time step.
     * \param parameter Partial w.r.t. which a parameter is to be computed
     * \param maximumDegreeAndOrder Maximum degree and order that are to be used.
     * \return Partial of state derivative w.r.t. given parameter
     */
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >& getCurrentParameterPartial(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
            const std::pair< int, int > maximumDegreeAndOrder )
    {
        if( currentDoubleParameterPartials_.count( std::make_pair( parameter, maximumDegreeAndOrder ) ) == 0 )
        {
            if( parameterDoublePartialFunctions_.count( std::make_pair( parameter, maximumDegreeAndOrder ) ) == 0 )
            {
                std::string errorMessage =
                        "Parameter of type " + std::to_string( parameter->getParameterName( ).first ) + ", " +
                           parameter->getParameterName( ).second.first + ", " +
                           parameter->getParameterName( ).second.second + " not found in list of existing partials";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                std::cerr << "Warning, double partial should already be calculatedz in Love number interface" << std::endl;
                currentDoubleParameterPartials_[ std::make_pair( parameter, maximumDegreeAndOrder ) ] =
                        parameterDoublePartialFunctions_.at( std::make_pair( parameter, maximumDegreeAndOrder ) )( );
            }
        }
        return currentDoubleParameterPartials_.at( std::make_pair( parameter, maximumDegreeAndOrder ) );
    }

    //! Function to retrieve a partial w.r.t. a double parameter
    /*!
     * Function to retrieve a partial w.r.t. a double parameter. An error is thrown if there is no dependency w.r.t.
     * the requested parameter. A warning is printed if the dependency exists, but has not yet been computed for the
     * current time step.
     * \param parameter Partial w.r.t. which a parameter is to be computed
     * \param maximumDegreeAndOrder Maximum degree and order that are to be used.
     * \return Partial of state derivative w.r.t. given parameter
     */
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >& getCurrentDoubleParameterPartial(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
            const std::pair< int, int > maximumDegreeAndOrder )
    {
        return getCurrentParameterPartial( parameter, maximumDegreeAndOrder );
    }


    //! Function to retrieve a partial w.r.t. a vector parameter
    /*!
     * Function to retrieve a partial w.r.t. a vector parameter. An error is thrown if there is no dependency w.r.t.
     * the requested parameter. A warning is printed if the dependency exists, but has not yet been computed for the
     * current time step.
     * \param parameter Partial w.r.t. which a parameter is to be computed
     * \param maximumDegreeAndOrder Maximum degree and order that are to be used.
     * \return Partial of state derivative w.r.t. given parameter
     */
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >& getCurrentParameterPartial(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
            const std::pair< int, int > maximumDegreeAndOrder  )
    {

        if( currentVectorParameterPartials_.count( std::make_pair( parameter, maximumDegreeAndOrder ) ) == 0 )
        {
            if( parameterVectorPartialFunctions_.count( std::make_pair( parameter, maximumDegreeAndOrder ) ) == 0 )
            {
                std::string errorMessage =
                        "Parameter of type " + std::to_string( parameter->getParameterName( ).first ) + ", " +
                           parameter->getParameterName( ).second.first + ", " +
                           parameter->getParameterName( ).second.second + " not found in list of existing partials";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                std::cerr << "Warning, vector partial should already be calculated in Love number interface" << std::endl;
                currentVectorParameterPartials_[ std::make_pair( parameter, maximumDegreeAndOrder ) ] =
                        parameterVectorPartialFunctions_.at( std::make_pair( parameter, maximumDegreeAndOrder ) )( );
            }
        }

        return currentVectorParameterPartials_.at( std::make_pair( parameter, maximumDegreeAndOrder ) );
    }

    //! Function to retrieve a partial w.r.t. a vector parameter
    /*!
     * Function to retrieve a partial w.r.t. a vector parameter. An error is thrown if there is no dependency w.r.t.
     * the requested parameter. A warning is printed if the dependency exists, but has not yet been computed for the
     * current time step.
     * \param parameter Partial w.r.t. which a parameter is to be computed
     * \param maximumDegreeAndOrder Maximum degree and order that are to be used.
     * \return Partial of state derivative w.r.t. given parameter
     */
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >& getCurrentVectorParameterPartial(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
            const std::pair< int, int > maximumDegreeAndOrder )
    {
        return getCurrentParameterPartial( parameter, maximumDegreeAndOrder );
    }



protected:

    //! Function to compute Love number partials from pre-computed partials and provided scaling values
    /*!
     * Function to compute Love number partials from pre-computed partials and provided scaling values. This function is
     * used to compute, for instance, the partials of complex Love numbers from the partials of the real components, which
     * are obtained from a simple sclaing
     * \param coefficientPartialsPerOrder Partial derivatives of C and S coefficients (row 0 and 1) that are to be scaled by
     * coefficientPartialScalers.
     * \param coefficientPartialScalers Pair of values by which to scale C and S coefficients (first and second of pair,
     * respectively).
     * \return Love number partials obtained from pre-computed partials and provided scaling values
     */
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > calculateSphericalHarmonicCoefficientPartialMatrix(
            const std::vector< Eigen::Vector2d >& coefficientPartialsPerOrder,
            const std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, Eigen::Matrix< double, 2, Eigen::Dynamic > >&
            coefficientPartialScalers );

    //! Function to calculate the partial of spherical harmonic coefficients w.r.t. real part of tidal love numbers.
    /*!
     *  Function to calculate the partial of spherical harmonic coefficients w.r.t. real part of tidal love numbers.
     *  Partials at a single degree are calculated,  with the desired orders required as input.
     *  \param degree Spherical harmonic degree at which partials are to be taken.
     *  \param orders Spherical harmonic orders at given degree at which partials are to be taken (if vector is empty,
     *  all orders at given degree are used,  in ascending order).
     *  \param deformingBodyIndices Indices of deforming bodies (in deformingBodies_ member vector) due to which the deformation
     *  is to be taken into account
     *  \param maximumDegree Maximum degree of the acceleration for which the partial derivatives are to be calculated.
     *  \param maximumOrder Maximum order of the acceleration for which the partial derivatives are to be calculated.
     *  \return Vector of partials of Spherical harmonic coefficients ([C_{nm};S_{nm}]). w.r.t. real part of tidal love number
     *  Entries of return vector correspond to orders given by corresponding entry in orders vector.
     */
    std::vector< Eigen::Vector2d > calculateCoefficientPartialWrtRealTidalLoveNumber(
            const int degree,
            const std::vector< int >& orders,
            const std::vector< int >& deformingBodyIndices,
            const int maximumDegree,
            const int maximumOrder );

    //! Function to pre-calculate all states of bodies involved.
    /*!
     *  Function to pre-calculate all states of bodies involved.
     *  Implemented here assuming no time lag, derived classes could provide modified implementations.
     *  \param deformingBodiesToUpdate Indices of deforming bodies (in deformingBodies_ member vector) due to which the
     *  deformationis to be taken into accountused.
     */
    virtual void updateCurrentTidalBodyStates( const std::vector< int >& deformingBodiesToUpdate );

    //! Function to set current state of single body and derived quantities
    /*!
     *  Function to set current state of single body and derived quantities
     *  Sets the variables used in calculation loop from pre-calculations of updateCurrentTidalBodyStates function.
     *  \param body Index of body for which variables are to be computed.
     *  \param order Degree of field at which updated is to be performed (not used; may be used in derived class function)
     *  \param degree Degree of field at which updated is to be performed (not used; may be used in derived class function)
     */
    virtual void setCurrentTidalBodyStates( const int degree, const int order, const int body );



    //! Reference radius used by tidal deformation model.
    double deformedBodyReferenceRadius_;



    //! Function to retrieve current state of body being tidally deformed
    /*!
     *  Function to retrieve current state of  body being tidally deformed. Typically linked to spherical harmonic acceleration
     *  model. Function is also used for tidal spherical harmonic variation calculations in case of no tidal time delay.
     */
    std::function< Eigen::Vector3d( ) > deformedBodyPositionFunction_;

    //! Function to retrieve current state of body being accelerated.
    std::function< Eigen::Vector3d( ) > acceleratingBodyPositionFunction_;



    //! Function to retrieve current gravitational parameter of body being deformed.
    std::function< double( ) > deformedBodyGravitationalParameterFunction_;


    //! Current gravitational parameter of body being deformed.
    double deformedBodyGravitationalParameter_;

    //! Vector of function to retrieve current gravitational parameters of bodies causing deformation.
    std::vector< std::function< double( ) > > deformingBodyGravitationalParameters_;

    //! Functions for retrieving states at current time of bodies causing deformation.
    std::vector< std::function< Eigen::Vector3d( ) > > deformingBodyStateFunctions_;

    //! Function returning the rotation from inertial to body-fixed frame of deformed body.
    std::function< Eigen::Quaterniond( ) > rotationToDeformedBodyFrameFrameFunction_;

    //! Current time to which object has been updated.
    double currentTime_;



    //! Current relative position in inertial frame of body causing deformation, as set by setCurrentTidalBodyStates function
    Eigen::Vector3d relativeDeformingBodyPosition_;

    //! Current relative position of body causing deformation in body-fixed frame of deformed body, as set by
    //! setCurrentTidalBodyStates function
    Eigen::Vector3d relativeDeformingBodySphericalPosition_;

    //! Current mass ratio of body causing, w.r.t. body undergoing, deformation, as set by setCurrentTidalBodyStates function
    double massRatio_;

    //! Current radius ratio of body causing, w.r.t. body undergoing, deformation, as set by setCurrentTidalBodyStates function
    double radiusRatio_;

    //! Current sine of latitude of body causing deformation, in frame fixed to body being deformed, as set by
    //! setCurrentTidalBodyStates function
    double sineOfLatitude_;

    //! Current sine of longitude of body causing deformation, in frame fixed to body being deformed, times i (sqrt(-1)) as set by
    //! setCurrentTidalBodyStates function
    std::complex< double > iLongitude_;


    //! Values used to scale real Love number partials to themselves (trivial, but required for interface)
    std::pair< Eigen::Matrix< double, 2, 1 >, Eigen::Matrix< double, 2, 1 > > realLoveNumberScaler_;

    //! Values used to scale real Love number partials to complex Love number partials
    std::pair< Eigen::Matrix< double, 2, 2 >, Eigen::Matrix< double, 2, 2 > > complexLoveNumberScaler_;

    std::vector< int > allDeformingBodyIndices_;

    //! Name of body being deformed
    std::string deformedBody_;

    //! List of bodies causing deformation
    std::vector< std::string > deformingBodies_;

    //! Current rotation from inertial frame to body-fixed frame of deformed body
    Eigen::Quaterniond rotationToTidallyDeformedBody_;

    //! Current position of deformed body
    Eigen::Vector3d positionOfTidallyDeformedBody_;

    //! List of current positions of bodies causing deformation
    std::vector< Eigen::Vector3d > positionsOfDeformingBodies_;


    //! List of functions to compute values of partials w.r.t. double parameter partials
    std::map< std::pair< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >,
    std::pair< int, int > >, std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > > currentDoubleParameterPartials_;

    //! List of functions to compute values of partials w.r.t. double parameter partials
    std::map< std::pair< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >,
    std::pair< int, int > >, std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >
    parameterDoublePartialFunctions_;

    //! Iterator over list of functions to compute values of partials w.r.t. double parameter partials
    std::map< std::pair< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >,
    std::pair< int, int > >, std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >::iterator
    parameterDoublePartialFunctionIterator_;


    //! List of current values of partials w.r.t. double parameter values (emptied at beginning of every time step).
    std::map< std::pair< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
    std::pair< int, int > >, std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > > currentVectorParameterPartials_;

    //! List of functions to compute values of partials w.r.t. vector parameter partials
    std::map< std::pair< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
    std::pair< int, int > >, std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >
    parameterVectorPartialFunctions_;

    //! Iterator over list of functions to compute values of partials w.r.t. vector parameter partials
    std::map< std::pair< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
    std::pair< int, int > >, std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >::iterator
    parameterVectorPartialFunctionIterator_;
};

}

}


#endif // TUDAT_TIDALLOVENUMBERPARTIALINTERFACE_H
