#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/tidalLoveNumberPartialInterface.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/tidalLoveNumber.h"

namespace tudat
{

namespace orbit_determination
{

void getMaximumUsedDegreeAndOrder(
        const int maximumDegree, const int maximumOrder, const int evaluationDegree,
        int& maximumUsedDegree, int& maximumUsedOrder )
{
    if( maximumDegree >= evaluationDegree )
    {
        maximumUsedDegree = evaluationDegree;
    }
    else
    {
        maximumUsedDegree = 0;
    }

    if( maximumOrder >= evaluationDegree )
    {
        maximumUsedOrder = evaluationDegree;
    }
    else
    {
        maximumUsedOrder = maximumOrder;
    }
}

std::vector< int > TidalLoveNumberPartialInterface::getSelectedDeformingBodyIds( const std::vector< std::string >& selectedBodyNames )
{
    std::vector< int > selectedBodyIds;
    for( unsigned int i = 0; i < selectedBodyNames.size( ); i++ )
    {
        std::vector< std::string >::iterator findIterator = std::find( deformingBodies_.begin( ), deformingBodies_.end( ), selectedBodyNames.at( i ) );
        if( findIterator != deformingBodies_.end( ) )
        {
            selectedBodyIds.push_back( std::distance( deformingBodies_.begin( ), findIterator ) );
        }
    }
    return selectedBodyIds;
}

std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >  TidalLoveNumberPartialInterface::calculateShericalHarmonicCoefficientPartialMatrix(
        const std::vector< Eigen::Vector2d >& coefficientPartialsPerOrder,
        const std::vector< int >& orders,
        const std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, Eigen::Matrix< double, 2, Eigen::Dynamic > >& coefficientPartialScalers )
{
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > coefficientsPartial;

    if( coefficientPartialScalers.first.rows( ) !=  coefficientPartialScalers.second.rows( ) )
    {
        std::cerr<<"Error when making sh acceleration partial from coefficient partial, scale matrices not of same size."<<std::endl;
    }


    for( unsigned int m = 0; m < coefficientPartialsPerOrder.size( ); m++ )
    {
        coefficientsPartial.push_back(
                    ( coefficientPartialScalers.first * coefficientPartialsPerOrder[ m ].x( ) +
                      coefficientPartialScalers.second * coefficientPartialsPerOrder[ m ].y( ) ) );
    }

    return coefficientsPartial;
}

//! Function to calculate the partial of sh acceleration wrt complex tidal love numbers.
std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > TidalLoveNumberPartialInterface::calculateShericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumbers(
        const int degree,
        const std::vector< int >& orders,
        const std::vector< int >& deformingBodyIndices,
        const int maximumDegree,
        const int maximumOrder )
{
    //std::cout<<"Calculating partial w.r.t. complex Love number"<<std::endl;

    deformedBodyGravitationalParameter_ = deformedBodyGravitationalParameterFunction_( );

    // Set and calculate states and rotation needed for calculation of partial.
    updateCurrentTidalBodyStates( degree, orders, deformingBodyIndices );

    // Calculate partial of the coefficients wrt love number.
    std::vector< Eigen::Vector2d > realCoefficientPartials =
            calculateCoefficientPartialWrtRealTidalLoveNumber( degree, orders, deformingBodyIndices, maximumDegree, maximumOrder );

    return calculateShericalHarmonicCoefficientPartialMatrix(
                realCoefficientPartials, orders, complexLoveNumberScaler_ );
}

//! Function to calculate the partial of sh acceleration wrt real tidal love numbers.
std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > TidalLoveNumberPartialInterface::calculateShericalHarmonicCoefficientsPartialWrtRealTidalLoveNumbers(
        const int degree,
        const std::vector< int >& orders,
        const std::vector< int >& deformingBodyIndices,
        const int maximumDegree,
        const int maximumOrder )
{
    deformedBodyGravitationalParameter_ = deformedBodyGravitationalParameterFunction_( );

    // Set and calculate states and rotation needed for calculation of partial.
    updateCurrentTidalBodyStates( degree, orders, deformingBodyIndices );

    // Calculate partial of the coefficients wrt love number.
    std::vector< Eigen::Vector2d > realCoefficientPartials = calculateCoefficientPartialWrtRealTidalLoveNumber(
                degree, orders, deformingBodyIndices, maximumDegree, maximumOrder );

    return calculateShericalHarmonicCoefficientPartialMatrix(
                realCoefficientPartials, orders, realLoveNumberScaler_ );
}

std::pair< int, std::pair< int, int > > TidalLoveNumberPartialInterface::setParameterPartialFunction(
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const int maximumDegree,
        const int maximumOrder )
{
    return std::make_pair( 0, std::make_pair( 0, 0 ) );
}

std::pair< int, std::pair< int, int > > TidalLoveNumberPartialInterface::setParameterPartialFunction(
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const int maximumDegree,
        const int maximumOrder )
{
    using namespace estimatable_parameters;

    int maximumUsedDegree = TUDAT_NAN;
    int maximumUsedOrder = TUDAT_NAN;

    int numberOfRows = 0;

    if( parameter->getParameterName( ).second.first == deformedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {

        case full_degree_tidal_love_number:
        {
            // Cast parameter object to required type.
            boost::shared_ptr< FullDegreeTidalLoveNumber > coefficientsParameter =
                    boost::dynamic_pointer_cast< FullDegreeTidalLoveNumber >( parameter );

            getMaximumUsedDegreeAndOrder( maximumDegree, maximumOrder, coefficientsParameter->getDegree( ), maximumUsedDegree,  maximumUsedOrder );

            if( maximumUsedDegree > 0 )
            {

                std::vector< int > selectedDeformingBodies = getSelectedDeformingBodyIds(
                            coefficientsParameter->getDeformingBodies( ) );

                if( selectedDeformingBodies.size( ) == coefficientsParameter->getDeformingBodies( ).size( ) &&
                        coefficientsParameter->getDeformingBodies( ).size( ) == deformingBodies_.size( ) )
                {
                    if( parameterVectorPartialFunctions_.count( std::make_pair( parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ) == 0 )
                    {
                        if( coefficientsParameter->useComplexComponents( ) )
                        {
                            // Calculate partials for complex love number
                            parameterVectorPartialFunctions_[ std::make_pair( parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ] =
                                    boost::bind( &TidalLoveNumberPartialInterface::calculateShericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumber, this,
                                                 coefficientsParameter->getDegree( ), selectedDeformingBodies, maximumUsedDegree, maximumUsedOrder );
                        }
                        else
                        {
                            // Calculate partial for real love number
                            parameterVectorPartialFunctions_[ std::make_pair( parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ] =
                                    boost::bind( &TidalLoveNumberPartialInterface::calculateShericalHarmonicCoefficientsPartialWrtRealTidalLoveNumber, this,
                                                 coefficientsParameter->getDegree( ), selectedDeformingBodies, maximumUsedDegree, maximumUsedOrder );
                        }
                    }

                    numberOfRows = coefficientsParameter->getParameterSize( );
                }
            }

            break;
        }
        case single_degree_variable_tidal_love_number:
        {
            // Cast parameter object to required type.
            boost::shared_ptr< SingleDegreeVariableTidalLoveNumber > coefficientsParameter =
                    boost::dynamic_pointer_cast< SingleDegreeVariableTidalLoveNumber >( parameter );

            getMaximumUsedDegreeAndOrder( maximumDegree, maximumOrder, coefficientsParameter->getDegree( ), maximumUsedDegree,  maximumUsedOrder );

            if( maximumUsedDegree > 0 )
            {
                std::vector< int > selectedDeformingBodies = getSelectedDeformingBodyIds(
                            coefficientsParameter->getDeformingBodies( ) );
                if( selectedDeformingBodies.size( ) == coefficientsParameter->getDeformingBodies( ).size( ) &&
                        coefficientsParameter->getDeformingBodies( ).size( ) == deformingBodies_.size( ) )
                {
                    if( parameterVectorPartialFunctions_.count( std::make_pair( parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ) == 0 )
                    {
                        if( coefficientsParameter->useComplexComponents( ) )
                        {
                            // Calculate partials for complex love number
                            parameterVectorPartialFunctions_[ std::make_pair( parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ] =
                                    boost::bind( &TidalLoveNumberPartialInterface::calculateShericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumbers, this,
                                                 coefficientsParameter->getDegree( ), coefficientsParameter->getOrders( ), selectedDeformingBodies,
                                                 maximumUsedDegree, maximumUsedOrder );
                        }
                        else
                        {
                            // Calculate partial for real love number
                            parameterVectorPartialFunctions_[ std::make_pair( parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ] =
                                    boost::bind( &TidalLoveNumberPartialInterface::calculateShericalHarmonicCoefficientsPartialWrtRealTidalLoveNumbers, this,
                                                 coefficientsParameter->getDegree( ), coefficientsParameter->getOrders( ), selectedDeformingBodies,
                                                 maximumUsedDegree, maximumUsedOrder );
                        }
                    }

                    numberOfRows = coefficientsParameter->getParameterSize( );
                }
            }

            break;
        }
        default:
            break;
        }
    }
    return std::make_pair( numberOfRows, std::make_pair( maximumUsedDegree, maximumUsedOrder ) );

}


void TidalLoveNumberPartialInterface::updateCurrentTidalBodyStates( const int degree,
                                                                    const std::vector< int >& orders,
                                                                    const std::vector< int >& deformingBodiesToUpdate )
{
    positionOfTidallyDeformedBody_ = deformedBodyPositionFunction_( );

    for( unsigned int i = 0; i < deformingBodiesToUpdate.size( ); i++ )
    {
        positionsOfDeformingBodies_[ deformingBodiesToUpdate[ i ] ] = deformingBodyStateFunctions_[ deformingBodiesToUpdate[ i ] ]( );
    }
}

void TidalLoveNumberPartialInterface::setCurrentTidalBodyStates( const int degree, const int order, const int body )
{
    {
        // Calculate properties of currently considered body
        massRatio_ = deformingBodyGravitationalParameters_[ body ]( ) / deformedBodyGravitationalParameter_;
        relativeDeformingBodyPosition_ = rotationToTidallyDeformedBody_ * (
                    positionsOfDeformingBodies_[ body ] - positionOfTidallyDeformedBody_ );
        relativeDeformingBodySphericalPosition_  =coordinate_conversions::
                convertCartesianToSpherical( relativeDeformingBodyPosition_ );
        radiusRatio_ = deformedBodyReferenceRadius_ / relativeDeformingBodySphericalPosition_.x( );
        iLongitude_ = std::complex< double >( 0.0, relativeDeformingBodySphericalPosition_.z( ) );
        sineOfLatitude_ = std::sin( mathematical_constants::PI / 2.0 - relativeDeformingBodySphericalPosition_.y( ) );
    }
}

std::vector< Eigen::Vector2d > TidalLoveNumberPartialInterface::calculateCoefficientPartialWrtRealTidalLoveNumber(
        const int degree,
        const std::vector< int >& orders,
        const std::vector< int >& deformingBodyIndices,
        const int maximumDegree,
        const int maximumOrder )
{
    std::vector< Eigen::Vector2d > realCoefficientPartials;
    realCoefficientPartials.resize( orders.size( ) );

    for( unsigned int i = 0; i < orders.size( ); i ++ )
    {
        realCoefficientPartials[ i ] = Eigen::Vector2d::Zero( );
    }

    if( degree <= maximumDegree )
    {
        std::complex< double > unitLoveNumberCoefficientVariations;

        // Iterate over all bodies causing deformation and calculate and add associated corrections
        for( unsigned int i = 0; i < deformingBodyIndices.size( ); i++ )
        {
            for( unsigned int m = 0; m < orders.size( ); m++ )
            {
                if( orders.at( m ) <= maximumOrder )
                {
                    setCurrentTidalBodyStates( degree, orders.at( m ), deformingBodyIndices.at( i ) );
                    // Calculate and add coefficients.
                    unitLoveNumberCoefficientVariations = gravitation::calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                                std::complex< double >( 1.0, 0.0 ), massRatio_,
                                basic_mathematics::raiseToIntegerPower( radiusRatio_, degree + 1 ) , sineOfLatitude_,
                                iLongitude_, degree, orders.at( m ) );
                    realCoefficientPartials[ m ].x( ) += unitLoveNumberCoefficientVariations.real( );
                    realCoefficientPartials[ m ].y( ) -= unitLoveNumberCoefficientVariations.imag( );
                }
            }
        }
    }
    return realCoefficientPartials;
}

}

}
