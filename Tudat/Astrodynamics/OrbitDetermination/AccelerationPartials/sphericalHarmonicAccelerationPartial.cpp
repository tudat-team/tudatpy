#include <omp.h>

#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicPartialFunctions.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/centralGravityAccelerationPartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

//! Contructor.
SphericalHarmonicsGravityPartial::SphericalHarmonicsGravityPartial(
        const std::string& acceleratedBody,
        const std::string& acceleratingBody,
        const boost::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > accelerationModel,
        const std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
        boost::shared_ptr< observation_partials::RotationMatrixPartial > >& rotationMatrixPartials ):
    AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::spherical_harmonic_gravity ),
    gravitationalParameterFunction_( accelerationModel->getGravitationalParameterFunction( ) ),
    bodyReferenceRadius_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getReferenceRadius, accelerationModel ) ),
    cosineCoefficients_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getCosineCoefficients, accelerationModel ) ),
    sineCoefficients_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getSineCoefficients, accelerationModel ) ),
    sphericalHarmonicCache_( accelerationModel->getSphericalHarmonicsCache( ) ),
    positionFunctionOfAcceleratedBody_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::
                                                     getCurrentPositionOfBodySubjectToAcceleration, accelerationModel ) ),
    positionFunctionOfAcceleratingBody_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::
                                                      getCurrentPositionOfBodyExertingAcceleration, accelerationModel ) ),
    toBodyFixedFrameRotation_( boost::bind( &reference_frames::QuaternionRotationWrapper::getRotationMatrix,
                                            boost::make_shared< reference_frames::QuaternionRotationWrapper >(
                                                accelerationModel->getToLocalFrameTransformation( ) ) ) ),
    accelerationFunction_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getAcceleration, accelerationModel ) ),
    updateFunction_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::updateMembers, accelerationModel, _1 ) ),
    rotationMatrixPartials_( rotationMatrixPartials ),
    accelerationUsesMutualAttraction_( accelerationModel->getIsMutualAttractionUsed( ) )
{

    // Update number of degrees and orders in legendre cache for calculation of position partials

    maximumDegree_ = accelerationModel->getCosineCoefficientFunction( )( ).rows( ) - 1;
    maximumOrder_ = accelerationModel->getCosineCoefficientFunction( )( ).cols( ) - 1;

    if( sphericalHarmonicCache_->getMaximumDegree( ) < maximumDegree_ || sphericalHarmonicCache_->getMaximumOrder( ) < maximumOrder_ + 2 )
    {
        sphericalHarmonicCache_->resetMaximumDegreeAndOrder( maximumDegree_, maximumOrder_ + 2 );
    }
}

//! Function to create a function returning a partial w.r.t. a double parameter.
std::pair< boost::function< Eigen::MatrixXd( ) >, int > SphericalHarmonicsGravityPartial::
getParameterPartialFunction( boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    // Declare return variables, default number of rows = 0 (i.e. no dependency)
    boost::function< Eigen::MatrixXd( ) > partialFunction;
    int numberOfRows = 0;


    // Check properties of body exerting acceleration.
    if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
    {
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunctionPair = getGravitationalParameterPartialFunction( parameter->getParameterName( ) );
        partialFunction = partialFunctionPair.first;
        numberOfRows = partialFunctionPair.second;
    }
    else if( parameter->getParameterName( ).second.first == acceleratingBody_ )
    {
        // Check if partial is a rotational property of body exerting acceleration.
        if( estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
        {
            // Check if required rotation matrix partial exists.
            if( rotationMatrixPartials_.count( std::make_pair( parameter->getParameterName( ).first,
                                                               parameter->getSecondaryIdentifier( ) ) ) != 0 )
            {
                // Get partial function.
                partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtRotationModelParameter,
                                               this, parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) );
                numberOfRows = 1;
            }
            else
            {
                std::cerr<<"Warning, not taking partial of sh acceleration wrt "<<parameter->getParameterName( ).first<<" of "<<
                           parameter->getParameterName( ).second.first<<" "<<std::endl;
            }
        }
        else if( estimatable_parameters::isParameterTidalProperty( parameter->getParameterName( ).first ) )
        {
            std::pair< int, std::pair< int, int > > currentTidalPartialOutput;

            boost::shared_ptr< estimatable_parameters::TidalLoveNumber< double > > tidalLoveNumber =
                    boost::dynamic_pointer_cast< estimatable_parameters::TidalLoveNumber< double >  >( parameter );
            if( tidalLoveNumber == NULL )
            {
                std::cerr<<"Error when getting tidal Love number vector parameter, object is NULL"<<std::endl;
            }

            int degree = tidalLoveNumber->getDegree( );
            std::vector< int > orders = tidalLoveNumber->getOrders( );
            int sumOrders = tidalLoveNumber->getSumOrders( );

            for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
            {
                currentTidalPartialOutput = tidalLoveNumberPartialInterfaces_.at( i )->setParameterPartialFunction(
                            parameter, maximumDegree_, maximumOrder_ );
                if( numberOfRows != 0 && currentTidalPartialOutput.first > 0 )
                {
                    std::cerr<<"Error A when getting tidal parameter partial, inconsistent output"<<
                               numberOfRows<<" "<<currentTidalPartialOutput.first<<std::endl;
                }
                else
                {
                    if( currentTidalPartialOutput.first > 0 )
                    {
                        boost::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunction =
                                boost::bind( &TidalLoveNumberPartialInterface::getCurrentDoubleParameterPartial, tidalLoveNumberPartialInterfaces_.at( i ),
                                             parameter, currentTidalPartialOutput.second );
                        partialFunction = boost::bind(
                                    &SphericalHarmonicsGravityPartial::wrtTidalModelParameter, this, coefficientPartialFunction, degree, orders, sumOrders,
                                    parameter->getParameterSize( ) );
                        numberOfRows = currentTidalPartialOutput.first;
                    }
                }
            }

        }
    }

    // Return partial function and partial size.
    return std::make_pair( partialFunction, numberOfRows );
}

//! Function to create a function returning a partial w.r.t. a vector parameter.
std::pair< boost::function< Eigen::MatrixXd( ) >, int > SphericalHarmonicsGravityPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace tudat::estimatable_parameters;

    // Declare return variables, default number of rows = 0 (i.e. no dependency)
    boost::function< Eigen::MatrixXd( ) > partialFunction;
    int numberOfRows = 0;

    // Check properties of body exerting acceleration.
    if( parameter->getParameterName( ).second.first == acceleratingBody_ )
    {
        // Check if partial is a rotational property of body exerting acceleration.
        if( estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
        {
            // Check if required rotation matrix partial exists.
            if( rotationMatrixPartials_.count( std::make_pair( parameter->getParameterName( ).first,
                                                               parameter->getSecondaryIdentifier( ) ) ) != 0 )
            {
                // Get partial function.
                partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtRotationModelParameter,
                                               this, parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) );
                numberOfRows = parameter->getParameterSize( );
            }
            else
            {
                std::cerr<<"Warning, not taking partial of sh acceleration wrt "<<parameter->getParameterName( ).first<<" of "<<
                           parameter->getParameterName( ).second.first<<std::endl;
            }
        }
        else if( estimatable_parameters::isParameterTidalProperty( parameter->getParameterName( ).first ) )
        {
            std::pair< int, std::pair< int, int > > currentTidalPartialOutput;

            boost::shared_ptr< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd > > tidalLoveNumber =
                    boost::dynamic_pointer_cast< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd >  >( parameter );
            if( tidalLoveNumber == NULL )
            {
                std::cerr<<"Error when getting tidal Love number vector parameter, object is NULL"<<std::endl;
            }

            int degree = tidalLoveNumber->getDegree( );
            std::vector< int > orders = tidalLoveNumber->getOrders( );
            int sumOrders = tidalLoveNumber->getSumOrders( );

            for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
            {
                currentTidalPartialOutput = tidalLoveNumberPartialInterfaces_.at( i )->setParameterPartialFunction( parameter, maximumDegree_, maximumOrder_ );
                if( numberOfRows != 0 && currentTidalPartialOutput.first > 0 )
                {
                    std::cerr<<"Error B when getting tidal parameter partial, inconsistent output "<<
                               numberOfRows<<" "<<currentTidalPartialOutput.first<<std::endl;
                }
                else
                {
                    if( currentTidalPartialOutput.first > 0 )
                    {
                        boost::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunction =
                                boost::bind( &TidalLoveNumberPartialInterface::getCurrentVectorParameterPartial, tidalLoveNumberPartialInterfaces_.at( i ), parameter, currentTidalPartialOutput.second );
                        partialFunction = boost::bind(
                                    &SphericalHarmonicsGravityPartial::wrtTidalModelParameter, this, coefficientPartialFunction, degree, orders, sumOrders,
                                    parameter->getParameterSize( ) );
                        numberOfRows = currentTidalPartialOutput.first;

                    }
                }
            }

        }
        // Check non-rotational parameters.
        else
        {
            switch( parameter->getParameterName( ).first )
            {
            case spherical_harmonics_cosine_coefficient_block:
            {
                // Cast parameter object to required type.
                boost::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                        boost::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );

                partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtCosineCoefficientBlock, this,
                                               coefficientsParameter->getBlockIndices( ), coefficientsParameter->getParameterSize( ) );
                numberOfRows = coefficientsParameter->getParameterSize( );

                break;
            }
            case spherical_harmonics_sine_coefficient_block:
            {
                // Cast parameter object to required type.

                boost::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                        boost::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );

                partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtSineCoefficientBlock, this,
                                               coefficientsParameter->getBlockIndices( ), coefficientsParameter->getParameterSize( ) );
                numberOfRows = coefficientsParameter->getParameterSize( );

                break;
            }
            case tabulated_gravity_field_coefficient_variations:
            {
                // Cast parameter object to required type.
                boost::shared_ptr< TabulatedGravityFieldCoefficientVariationsParameter > tabulatedCoefficientParameter =
                        boost::dynamic_pointer_cast< TabulatedGravityFieldCoefficientVariationsParameter >( parameter );

                partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtTabulatedSphericalHarmonicCoefficients,
                                               this, tabulatedCoefficientParameter );
                numberOfRows = tabulatedCoefficientParameter->getParameterSize( );

                break;
            }
            case periodic_gravity_field_variations:
            {
                boost::shared_ptr< PeriodicGravityFieldVariationsParameter > periodicVariationsParameter =
                        boost::dynamic_pointer_cast< PeriodicGravityFieldVariationsParameter >( parameter );

                partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtPeriodicCoefficientVariations,
                                               this, periodicVariationsParameter );
                numberOfRows = periodicVariationsParameter->getParameterSize( );
                break;
            }
            default:
                break;
            }
        }
    }

    // Return partial function and partial size.
    return std::make_pair( partialFunction, numberOfRows );
}

std::pair< boost::function< Eigen::MatrixXd( ) >, int >
SphericalHarmonicsGravityPartial::getGravitationalParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterId )
{
    boost::function< Eigen::MatrixXd( ) > partialFunction;
    int numberOfColumns = 0;

    if( parameterId.first ==  estimatable_parameters::gravitational_parameter )
    {
        // Check for dependency
        if( parameterId.second.first == acceleratingBody_ )
        {
            partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtGravitationalParameterOfCentralBody,
                                           this );
            numberOfColumns = 1;

        }

        if( parameterId.second.first == acceleratedBody_ )
        {
            if( accelerationUsesMutualAttraction_ )
            {
                partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtGravitationalParameterOfCentralBody,
                                               this );
                numberOfColumns = 1;
            }
        }
    }
    return std::make_pair( partialFunction, numberOfColumns );
}

//! Function for updating the partial object to current state and time.
void SphericalHarmonicsGravityPartial::update( const double currentTime )
{
    using namespace tudat::coordinate_conversions;


    if( !( currentTime_ == currentTime ) )
    {
        // Update acceleration model
        updateFunction_( currentTime );

        // Calculate Cartesian position in frame fixed to body exerting acceleration
        Eigen::Matrix3d currentRotationToBodyFixedFrame_ = toBodyFixedFrameRotation_( );
        Eigen::Vector3d bodyFixedPosition =
                currentRotationToBodyFixedFrame_ * ( positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( ) );

        // Calculate spherical position in frame fixed to body exerting acceleration
        bodyFixedSphericalPosition_ = convertCartesianToSpherical( bodyFixedPosition );
        bodyFixedSphericalPosition_( 1 ) = mathematical_constants::PI / 2.0 - bodyFixedSphericalPosition_( 1 );

        // Get spherical harmonic coefficients
        currentCosineCoefficients_ = cosineCoefficients_( );
        currentSineCoefficients_ = sineCoefficients_( );

        // Update trogonometric functions of multiples of longitude.
        sphericalHarmonicCache_->update(
                    bodyFixedSphericalPosition_( 0 ), std::sin( bodyFixedSphericalPosition_( 1 ) ),
                    bodyFixedSphericalPosition_( 2 ), bodyReferenceRadius_( ) );


        // Calculate partial of acceleration wrt position of body undergoing acceleration.
        currentBodyFixedPartialWrtPosition_ = computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
                    bodyFixedPosition, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
                    currentCosineCoefficients_, currentSineCoefficients_, sphericalHarmonicCache_ );

        currentPartialWrtVelocity_.setZero( );
        currentPartialWrtPosition_ =
                currentRotationToBodyFixedFrame_.inverse( ) * currentBodyFixedPartialWrtPosition_ * currentRotationToBodyFixedFrame_;

        currentTime_ = currentTime;
    }
}

//! Function to calculate the partial of the acceleration wrt a set of cosine coefficients.
Eigen::Matrix< double, 3, Eigen::Dynamic > SphericalHarmonicsGravityPartial::wrtCosineCoefficientBlock(
        const std::map< int, std::pair< int, int > >& blockIndices, const int numberOfParameters )
{

    // Initialize partial matrix to zero values.
    Eigen::Matrix< double, 3, Eigen::Dynamic > partialMatrix =
            Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, numberOfParameters );

    // Pre-calculate multiplicative term found in all partial terms.
    double commonMultiplier = gravitationalParameterFunction_( ) / ( bodyReferenceRadius_( ) * bodyReferenceRadius_( ) );

    // Get rotation to global frame.
    Eigen::Matrix3d toGlobalFrame = toBodyFixedFrameRotation_( ).inverse( );


    //std::cout<<bodyFixedSphericalPosition_.transpose( )<<std::endl;

    // Iterate over all coefficient degree and orders for which a partial is to be calculated.
    int currentIndex = 0;
    for( std::map< int, std::pair< int, int > >::const_iterator blockIterator = blockIndices.begin( );
         blockIterator != blockIndices.end( ); blockIterator++ )
    {
        // Iterate over all required orders in current degree.
        for( int order = blockIterator->second.first; order < blockIterator->second.second + blockIterator->second.first; order++ )
        {
            if( blockIterator->first <= maximumDegree_ && order <= maximumOrder_ )
            {
                // Calculate and set partial of current degree and order.
                partialMatrix.block( 0, currentIndex, 3, 1 ) = toGlobalFrame * calculateSphericalHarmonicGravityWrtCCoefficient(
                            commonMultiplier, blockIterator->first, order, bodyFixedSphericalPosition_, bodyReferenceRadius_( ),
                            legendreCache_ );
            }
            currentIndex++;
        }
    }

    return partialMatrix;
}

//! Function to calculate the partial of the acceleration wrt a set of sine coefficients.
Eigen::Matrix< double, 3, Eigen::Dynamic > SphericalHarmonicsGravityPartial::wrtSineCoefficientBlock(
        const std::map< int, std::pair< int, int > >& blockIndices, const int numberOfParameters )
{
    // Initialize partial matrix to zero values.
    Eigen::Matrix< double, 3, Eigen::Dynamic > partialMatrix =
            Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, numberOfParameters );

    // Pre-calculate multiplicative term found in all partial terms.
    double commonMultiplier = gravitationalParameterFunction_( ) / ( bodyReferenceRadius_( ) * bodyReferenceRadius_( ) );

    // Get rotation to global frame.
    Eigen::Matrix3d toGlobalFrame = toBodyFixedFrameRotation_( ).inverse( );

    // Iterate over all coefficient degree and orders for which a partial is to be calculated.
    int currentIndex = 0;
    for( std::map< int, std::pair< int, int > >::const_iterator blockIterator = blockIndices.begin( );
         blockIterator != blockIndices.end( ); blockIterator++ )
    {
        // Iterate over all required orders in current degree.
        for( int order = blockIterator->second.first; order < blockIterator->second.second +
             blockIterator->second.first; order++ )
        {
            if( blockIterator->first <= maximumDegree_ && order <= maximumOrder_ )
            {
                // Calculate and set partial of current degree and order.
                partialMatrix.block( 0, currentIndex, 3, 1 ) = toGlobalFrame * calculateSphericalHarmonicGravityWrtSCoefficient(
                            commonMultiplier, blockIterator->first, order, bodyFixedSphericalPosition_, bodyReferenceRadius_( ),
                            legendreCache_ );
            }

            currentIndex++;
        }
    }

    return partialMatrix;
}

//! Function to calculate an acceleration partial wrt a rotational parameter.
Eigen::Matrix< double, 3, Eigen::Dynamic > SphericalHarmonicsGravityPartial::wrtRotationModelParameter(
        const estimatable_parameters::EstimatebleParametersEnum parameterType,
        const std::string& secondaryIdentifier )
{
    // Calculate distance vector between bodies.
    Eigen::Vector3d distanceVector = positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( );

    // Get rotation matrix partial(s) wrt requested parameter
    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
            rotationMatrixPartials_.at( std::make_pair( parameterType, secondaryIdentifier ) )->
            calculatePartialOfRotationMatrixWrParameter( currentTime_ );

    // Initialze return partial matrix to zero.
    Eigen::Matrix< double, 3, Eigen::Dynamic > accelerationPartial = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero(
                3, rotationMatrixPartials.size( ) );

    // Iterate for each single parameter entry partial.
    for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
    {
        // Calculate acceleration partial for current parameter entry.
        accelerationPartial.block( 0, i, 3, 1 ) = rotationMatrixPartials[ i ] * toBodyFixedFrameRotation_( ) * accelerationFunction_( ) +
                toBodyFixedFrameRotation_( ).inverse( ) * currentBodyFixedPartialWrtPosition_* rotationMatrixPartials[ i ].transpose( ) * distanceVector;
    }

    return accelerationPartial;
}

}

}

}

