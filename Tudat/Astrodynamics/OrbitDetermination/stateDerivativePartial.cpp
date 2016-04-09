#include "Tudat/Astrodynamics/OrbitDetermination/stateDerivativePartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{


void evaluateNegativeParameterPartialFunction( const boost::function< void( Eigen::MatrixXd& ) > parameterPartialFunction,
                                               Eigen::MatrixXd& partial )
{
    parameterPartialFunction( partial );
    partial *= -1.0;
}

void evaluateSubtractedParameterPartialFunction( const boost::function< void( Eigen::MatrixXd& ) > firstParameterPartialFunction,
                                                 const boost::function< void( Eigen::MatrixXd& ) > parameterPartialFunctionToSubtract,
                                                 Eigen::MatrixXd& partial )
{
    firstParameterPartialFunction( partial );

    Eigen::MatrixXd subtractedPartial = Eigen::MatrixXd::Zero( partial.rows( ), partial.cols( ) );
    parameterPartialFunctionToSubtract( subtractedPartial );

    partial -= subtractedPartial;
}

void evaluateAddedParameterPartialFunction( const boost::function< void( Eigen::MatrixXd& ) > firstParameterPartialFunction,
                                            const boost::function< void( Eigen::MatrixXd& ) > parameterPartialFunctionToAdd,
                                            Eigen::MatrixXd& partial )
{
    firstParameterPartialFunction( partial );

    Eigen::MatrixXd addedPartial = Eigen::MatrixXd::Zero( partial.rows( ), partial.cols( ) );
    parameterPartialFunctionToAdd( addedPartial );

    partial += addedPartial;
}

std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > createMergedParameterPartialFunction(
        const std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >& partialFunctionOfAccelerationToAdd,
        const std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >& partialFunctionOfAccelerationToSubtract )
{
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >  parameterPartialFunction;

    if( ( partialFunctionOfAccelerationToAdd.second == 0 ) &&
            ( partialFunctionOfAccelerationToSubtract.second == 0 ) )
    {
        parameterPartialFunction = std::make_pair( boost::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }
    else if( partialFunctionOfAccelerationToSubtract.second == 0 )
    {
        parameterPartialFunction = partialFunctionOfAccelerationToAdd;
    }
    else if( partialFunctionOfAccelerationToAdd.second == 0 )
    {
        parameterPartialFunction = std::make_pair( boost::bind(
                                                       &evaluateNegativeParameterPartialFunction,
                                                       partialFunctionOfAccelerationToSubtract.first, _1 ),
                                                   partialFunctionOfAccelerationToSubtract.second );
    }
    else if( partialFunctionOfAccelerationToSubtract.second !=
             partialFunctionOfAccelerationToAdd.second )
    {
        throw std::runtime_error( "Error when making merged parameter partial function, separate functions are incompatible" );
    }
    else
    {
        parameterPartialFunction = std::make_pair( boost::bind( &evaluateSubtractedParameterPartialFunction,
                                                               partialFunctionOfAccelerationToAdd.first,
                                                               partialFunctionOfAccelerationToSubtract.first, _1 ),
                                                  partialFunctionOfAccelerationToSubtract.second );
    }
    return parameterPartialFunction;
}

boost::function< void( Eigen::MatrixXd& ) > getCombinedCurrentDoubleParameterFunction(
        const boost::shared_ptr< StateDerivativePartial > firstPartial,
        const boost::shared_ptr< StateDerivativePartial > secondPartial,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterObject,
        const int firstPartialSize, const int secondPartialSize,
        const bool subtractPartials )
{
    boost::function< void( Eigen::MatrixXd& ) > partialFunction;

    boost::function< void( Eigen::MatrixXd& ) > firstPartialFunction =
            boost::bind( &StateDerivativePartial::getCurrentDoubleParameterPartial, firstPartial, parameterObject, _1 );
    boost::function< void( Eigen::MatrixXd& ) > secondPartialFunction =
            boost::bind( &StateDerivativePartial::getCurrentDoubleParameterPartial, secondPartial, parameterObject, _1 );

    if( firstPartialSize == 0 && secondPartialSize == 0 )
    {
        throw std::runtime_error( "Error when getting combined current partial size, both partials have size zero " );
    }
    else if( secondPartialSize == 0 )
    {
        partialFunction = firstPartialFunction;
    }
    else if( firstPartialSize == 0 )
    {
        if( subtractPartials )
        {
            partialFunction = boost::bind( &evaluateNegativeParameterPartialFunction,
                                           secondPartialFunction, _1 );
        }
        else
        {
            partialFunction = secondPartialFunction;

        }
    }
    else if( firstPartialSize == secondPartialSize )
    {
        if( subtractPartials )
        {
            partialFunction = boost::bind( &evaluateSubtractedParameterPartialFunction,
                                           firstPartialFunction,
                                           secondPartialFunction, _1 );
        }
        else
        {
            partialFunction = boost::bind( &evaluateAddedParameterPartialFunction,
                                           firstPartialFunction,
                                           secondPartialFunction, _1 );
        }

    }
    else
    {
        throw std::runtime_error( "Error when getting combined current partial size, both partials have different non-zero size." );
    }
    return partialFunction;
}

boost::function< void( Eigen::MatrixXd& ) > getCombinedCurrentVectorParameterFunction(
        const boost::shared_ptr< StateDerivativePartial > firstPartial,
        const boost::shared_ptr< StateDerivativePartial > secondPartial,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterObject,
        const int firstPartialSize, const int secondPartialSize,
        const bool subtractPartials )
{
    boost::function< void( Eigen::MatrixXd& ) > partialFunction;

    boost::function< void( Eigen::MatrixXd& ) > firstPartialFunction =
            boost::bind( &StateDerivativePartial::getCurrentVectorParameterPartial, firstPartial, parameterObject, _1 );
    boost::function< void( Eigen::MatrixXd& ) > secondPartialFunction =
            boost::bind( &StateDerivativePartial::getCurrentVectorParameterPartial, secondPartial, parameterObject, _1 );


    if( firstPartialSize == 0 && secondPartialSize == 0 )
    {
        throw std::runtime_error( "Error when getting combined current partial size, both partials have size zero " );
    }
    else if( secondPartialSize == 0 )
    {
        partialFunction = firstPartialFunction;
    }
    else if( firstPartialSize == 0 )
    {
        if( subtractPartials )
        {
            partialFunction = boost::bind( &evaluateNegativeParameterPartialFunction,
                                           secondPartialFunction, _1 );
        }
        else
        {
            partialFunction = secondPartialFunction;

        }
    }
    else if( firstPartialSize == secondPartialSize )
    {
        if( subtractPartials )
        {
            partialFunction = boost::bind( &evaluateSubtractedParameterPartialFunction,
                                           firstPartialFunction,
                                           secondPartialFunction, _1 );
        }
        else
        {
            partialFunction = boost::bind( &evaluateAddedParameterPartialFunction,
                                           firstPartialFunction,
                                           secondPartialFunction, _1 );
        }

    }
    else
    {
        throw std::runtime_error( "Error when getting combined current partial size, both partials have different non-zero size." );
    }
    return partialFunction;
}

}

}

}

