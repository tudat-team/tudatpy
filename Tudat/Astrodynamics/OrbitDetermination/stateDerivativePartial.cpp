#include "Tudat/Astrodynamics/OrbitDetermination/stateDerivativePartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{


Eigen::MatrixXd evaluateNegativeParameterPartialFunction( const boost::function< Eigen::MatrixXd( ) > parameterPartialFunction )
{
    return -parameterPartialFunction( );
}

Eigen::MatrixXd evaluateSubtractedParameterPartialFunction( const boost::function< Eigen::MatrixXd( ) > firstParameterPartialFunction,
                                                            const boost::function< Eigen::MatrixXd( ) > parameterPartialFunctionToSubtract )
{
    return firstParameterPartialFunction( ) - parameterPartialFunctionToSubtract( );
}

Eigen::MatrixXd evaluateAddedParameterPartialFunction( const boost::function< Eigen::MatrixXd( ) > firstParameterPartialFunction,
                                                       const boost::function< Eigen::MatrixXd( ) > parameterPartialFunctionToAdd )
{
    return firstParameterPartialFunction( ) +  parameterPartialFunctionToAdd( );
}

std::pair< boost::function< Eigen::MatrixXd( ) >, int > createMergedParameterPartialFunction(
        const std::pair< boost::function< Eigen::MatrixXd( ) >, int >& partialFunctionOfAccelerationToAdd,
        const std::pair< boost::function< Eigen::MatrixXd( ) >, int >& partialFunctionOfAccelerationToSubtract )
{
    std::pair< boost::function< Eigen::MatrixXd( ) >, int >  parameterPartialFunction;

    if( ( partialFunctionOfAccelerationToAdd.second == 0 ) &&
            ( partialFunctionOfAccelerationToSubtract.second == 0 ) )
    {
        parameterPartialFunction = std::make_pair( boost::function< Eigen::MatrixXd( ) >( ), 0 );
    }
    else if( partialFunctionOfAccelerationToSubtract.second == 0 )
    {
        parameterPartialFunction = partialFunctionOfAccelerationToAdd;
    }
    else if( partialFunctionOfAccelerationToAdd.second == 0 )
    {
        parameterPartialFunction = std::make_pair( boost::bind( &evaluateNegativeParameterPartialFunction,
                                                                partialFunctionOfAccelerationToSubtract.first ),
                                                   partialFunctionOfAccelerationToSubtract.second );
    }
    else if( partialFunctionOfAccelerationToSubtract.second !=
             partialFunctionOfAccelerationToAdd.second )
    {
        std::cerr<<"Error when making merged parameter partial function, ";
        std::cout<<"separate functions are incompatible"<<std::endl;
    }
    else
    {
        parameterPartialFunction = std::make_pair( boost::bind( &evaluateSubtractedParameterPartialFunction,
                                                               partialFunctionOfAccelerationToAdd.first,
                                                               partialFunctionOfAccelerationToSubtract.first ),
                                                  partialFunctionOfAccelerationToSubtract.second );
    }
    return parameterPartialFunction;
}

boost::function< Eigen::MatrixXd( ) > getCombinedCurrentDoubleParameterFunction(
        const boost::shared_ptr< StateDerivativePartial > firstPartial,
        const boost::shared_ptr< StateDerivativePartial > secondPartial,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterObject,
        const int firstPartialSize, const int secondPartialSize,
        const bool subtractPartials )
{
    boost::function< Eigen::MatrixXd( ) > partialFunction;

    boost::function< Eigen::MatrixXd( ) > firstPartialFunction =
            boost::bind( &StateDerivativePartial::getCurrentDoubleParameterPartial, firstPartial, parameterObject );
    boost::function< Eigen::MatrixXd( ) > secondPartialFunction =
            boost::bind( &StateDerivativePartial::getCurrentDoubleParameterPartial, secondPartial, parameterObject );

    if( firstPartialSize == 0 && secondPartialSize == 0 )
    {
        std::cerr<<"Error when getting combined current partial size, both partials have size zero "<<std::endl;
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
                                           secondPartialFunction );
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
                                           secondPartialFunction );
        }
        else
        {
            partialFunction = boost::bind( &evaluateAddedParameterPartialFunction,
                                           firstPartialFunction,
                                           secondPartialFunction );
        }

    }
    else
    {
        std::cerr<<"Error when getting combined current partial size, both partials have different non-zero size."<<std::endl;
    }
    return partialFunction;
}

boost::function< Eigen::MatrixXd( ) > getCombinedCurrentVectorParameterFunction(
        const boost::shared_ptr< StateDerivativePartial > firstPartial,
        const boost::shared_ptr< StateDerivativePartial > secondPartial,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterObject,
        const int firstPartialSize, const int secondPartialSize,
        const bool subtractPartials )
{
    boost::function< Eigen::MatrixXd( ) > partialFunction;

    boost::function< Eigen::MatrixXd( ) > firstPartialFunction =
            boost::bind( &StateDerivativePartial::getCurrentVectorParameterPartial, firstPartial, parameterObject );
    boost::function< Eigen::MatrixXd( ) > secondPartialFunction =
            boost::bind( &StateDerivativePartial::getCurrentVectorParameterPartial, secondPartial, parameterObject );


    if( firstPartialSize == 0 && secondPartialSize == 0 )
    {
        std::cerr<<"Error when getting combined current partial size, both partials have size zero "<<std::endl;
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
                                           secondPartialFunction );
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
                                           secondPartialFunction );
        }
        else
        {
            partialFunction = boost::bind( &evaluateAddedParameterPartialFunction,
                                           firstPartialFunction,
                                           secondPartialFunction );
        }

    }
    else
    {
        std::cerr<<"Error when getting combined current partial size, both partials have different non-zero size."<<std::endl;
    }
    return partialFunction;
}

}

}

}

