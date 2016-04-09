#ifndef TUDAT_STATEDERIVATIVEPARTIAL_H
#define TUDAT_STATEDERIVATIVEPARTIAL_H

#include <string>
#include <map>
#include <Eigen/Core>

#include <boost/bind.hpp>
#include <boost/assign/list_of.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"


namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{


class StateDerivativePartial
{

public:

    StateDerivativePartial( const propagators::IntegratedStateType integratedStateType,
                            const std::pair< std::string, std::string >& integrationReferencePoint ):
        integratedStateType_( integratedStateType ), integrationReferencePoint_( integrationReferencePoint )
    {
        currentTime_ = TUDAT_NAN;
    }

    virtual ~StateDerivativePartial( ) { }

//    Eigen::MatrixXd wrtStateOfIntegratedBody(
//            const std::pair< std::string, std::string >& stateReferencePoint,
//            const propagators::IntegratedStateType integratedStateType )
//    {
//        return getDerivativeFunctionWrtStateOfIntegratedBody( stateReferencePoint, integratedStateType ).first( );
//    }

    virtual std::pair< boost::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >  getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) = 0;

    virtual bool isStateDerivativeDependentOnIntegratedState(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) = 0;

    Eigen::MatrixXd wrtParameter(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        Eigen::MatrixXd partial;

        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunction = getParameterPartialFunction( parameter );
        if( partialFunction.second > 0 )
        {
            partial = partialFunction.first( );
        }
        else
        {
            partial = Eigen::MatrixXd::Zero( propagators::getSingleIntegrationSize( integratedStateType_ ), 1 );
        }
        return partial;
    }


    virtual std::pair< boost::function< Eigen::MatrixXd( ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        boost::function< Eigen::MatrixXd( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    Eigen::MatrixXd wrtParameter( boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        Eigen::MatrixXd partial;
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > partialFunction = getParameterPartialFunction( parameter );
        if( partialFunction.second > 0 )
        {
            partial = partialFunction.first( );
        }
        else
        {
            partial = Eigen::MatrixXd::Zero( propagators::getSingleIntegrationSize( integratedStateType_ ), parameter->getParameterSize( ) );
        }
        return partial;
    }


    virtual std::pair< boost::function< Eigen::MatrixXd( ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        boost::function< Eigen::MatrixXd( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Pure virtual for updating the partial object to current state.
    /*!
     *  Pure virtual for updating the partial object to current state.
     */
    virtual void update( const double currentTime ) = 0;

    propagators::IntegratedStateType getIntegratedStateType( )
    {
        return integratedStateType_;
    }

    std::pair< std::string, std::string > getIntegrationReferencePoint( )
    {
        return integrationReferencePoint_;
    }

    void resetTime( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime  ) )
        {
            resetCurrentParameterValues( );
            currentTime_ = currentTime;
        }        

        resetTimeOfMemberObjects( );

    }

    Eigen::MatrixXd& getCurrentParameterPartial(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        if( currentDoubleParameterPartials_.count( parameter ) == 0 )
        {
            if( parameterDoublePartialFunctions_.count( parameter ) == 0 )
            {
                throw std::runtime_error(
                            "Parameter of type " +
                            boost::lexical_cast< std::string >( parameter->getParameterName( ).first ) + ", " +
                            parameter->getParameterName( ).second.first + ", " +
                            parameter->getParameterName( ).second.second + ", " + " not found in list of existing partials" );
            }
            else
            {
                std::cerr<<"Warning, double partial should already be calculated"<<std::endl;
                currentDoubleParameterPartials_[ parameter ] = parameterDoublePartialFunctions_.at( parameter )( );
            }
        }
        return currentDoubleParameterPartials_[ parameter ];
    }

    Eigen::MatrixXd& getCurrentDoubleParameterPartial(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        return getCurrentParameterPartial( parameter );
    }

    Eigen::MatrixXd& getCurrentParameterPartial(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        if( currentVectorParameterPartials_.count( parameter ) == 0 )
        {
            if( parameterVectorPartialFunctions_.count( parameter ) == 0 )
            {
                std::cerr<<"Parameter of type "<<parameter->getParameterName( ).first<<", "<<
                           parameter->getParameterName( ).second.first<<", "<<
                           parameter->getParameterName( ).second.second<<" not found in list of existing partials"<<std::endl;
            }
            else
            {
                std::cerr<<"Warning, vector partial should already be calculated"<<std::endl;
                currentVectorParameterPartials_[ parameter ] = parameterVectorPartialFunctions_.at( parameter )( );
            }
        }
        return currentVectorParameterPartials_[ parameter ];
    }

    Eigen::MatrixXd& getCurrentVectorParameterPartial(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        return getCurrentParameterPartial( parameter );
    }

    virtual int setParameterPartialUpdateFunction(
                boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > parameterPartialFunction = getParameterPartialFunction( parameter );
        if( parameterPartialFunction.second > 0 && parameterDoublePartialFunctions_.count( parameter ) == 0 )
        {
            if( parameterDoublePartialFunctions_.count( parameter ) == 0 )
            {
                parameterDoublePartialFunctions_[ parameter ] = parameterPartialFunction.first;
                doesCurrentDoubleParameterPartialExist_[ parameter ] = 0;
            }
        }

        return parameterPartialFunction.second;
    }

    virtual int setParameterPartialUpdateFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::pair< boost::function< Eigen::MatrixXd( ) >, int > parameterPartialFunction = getParameterPartialFunction( parameter );
        if( parameterPartialFunction.second > 0 && currentVectorParameterPartials_.count( parameter ) == 0 )
        {
            if( parameterVectorPartialFunctions_.count( parameter ) == 0 )
            {
                parameterVectorPartialFunctions_[ parameter ] = parameterPartialFunction.first;
                doesCurrentVectorParameterPartialExist_[ parameter ] = 0;
            }
        }

        return parameterPartialFunction.second;
    }

    void updateParameterPartials( )
    {
        updateParameterPartialsOfMemberObjects( );

        for( parameterDoublePartialFunctionIterator_ = parameterDoublePartialFunctions_.begin( );
             parameterDoublePartialFunctionIterator_ != parameterDoublePartialFunctions_.end( );
             parameterDoublePartialFunctionIterator_++ )
        {
            if( doesCurrentDoubleParameterPartialExist_.at( parameterDoublePartialFunctionIterator_->first ) == 0 )
            {
                currentDoubleParameterPartials_[ parameterDoublePartialFunctionIterator_->first ] = parameterDoublePartialFunctionIterator_->second( );
                doesCurrentDoubleParameterPartialExist_[ parameterDoublePartialFunctionIterator_->first ] = 1;
            }
        }

        for( parameterVectorPartialFunctionIterator_ = parameterVectorPartialFunctions_.begin( );
             parameterVectorPartialFunctionIterator_ != parameterVectorPartialFunctions_.end( );
             parameterVectorPartialFunctionIterator_++ )
        {
            //std::cout<<"New vec. partial: "<<parameterVectorPartialFunctionIterator_->first->getParameterName( ).first<<std::endl;
            if( doesCurrentVectorParameterPartialExist_.at( parameterVectorPartialFunctionIterator_->first ) == 0 )
            {
                currentVectorParameterPartials_[ parameterVectorPartialFunctionIterator_->first ] = parameterVectorPartialFunctionIterator_->second( );
                doesCurrentVectorParameterPartialExist_[  parameterVectorPartialFunctionIterator_->first ] = 1;
            }
        }
    }

protected:

    virtual void updateParameterPartialsOfMemberObjects( )
    {

    }

    virtual void resetTimeOfMemberObjects( )
    {

    }

    void resetCurrentParameterValues( )
    {
        for( doesDoubleParameterExistIterator_ = doesCurrentDoubleParameterPartialExist_.begin( );
             doesDoubleParameterExistIterator_ !=  doesCurrentDoubleParameterPartialExist_.end( );
             doesDoubleParameterExistIterator_++ )
        {
            doesCurrentDoubleParameterPartialExist_[ doesDoubleParameterExistIterator_->first ] = 0;
        }

        for( doesVectorParameterExistIterator_ = doesCurrentVectorParameterPartialExist_.begin( );
             doesVectorParameterExistIterator_ !=  doesCurrentVectorParameterPartialExist_.end( );
             doesVectorParameterExistIterator_++ )
        {
            doesCurrentVectorParameterPartialExist_[ doesVectorParameterExistIterator_->first ] = 0;
        }
    }


    propagators::IntegratedStateType integratedStateType_;

    std::pair< std::string, std::string > integrationReferencePoint_;


    std::map< boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > >, bool > doesCurrentDoubleParameterPartialExist_;

    std::map< boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > >, bool >::iterator doesDoubleParameterExistIterator_;


    std::map< boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > >, Eigen::MatrixXd > currentDoubleParameterPartials_;

    std::map< boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > >, boost::function< Eigen::MatrixXd( ) > > parameterDoublePartialFunctions_;

    std::map< boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > >, boost::function< Eigen::MatrixXd( ) > >::iterator parameterDoublePartialFunctionIterator_;


    std::map< boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >, bool > doesCurrentVectorParameterPartialExist_;

    std::map< boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >, bool >::iterator doesVectorParameterExistIterator_;


    std::map< boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >, Eigen::MatrixXd > currentVectorParameterPartials_;

    std::map< boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >, boost::function< Eigen::MatrixXd( ) > > parameterVectorPartialFunctions_;

    std::map< boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >, boost::function< Eigen::MatrixXd( ) > >::iterator parameterVectorPartialFunctionIterator_;


    double currentTime_;

};

typedef std::vector< std::vector< boost::shared_ptr< partial_derivatives::StateDerivativePartial > > > StateDerivativePartialsMap;

Eigen::MatrixXd evaluateNegativeParameterPartialFunction( const boost::function< Eigen::MatrixXd( ) > parameterPartialFunction );

Eigen::MatrixXd evaluateSubtractedParameterPartialFunction( const boost::function< Eigen::MatrixXd( ) > firstParameterPartialFunction,
                                                            const boost::function< Eigen::MatrixXd( ) > parameterPartialFunctionToSubtract );

Eigen::MatrixXd evaluateAddedParameterPartialFunction( const boost::function< Eigen::MatrixXd( ) > firstParameterPartialFunction,
                                                       const boost::function< Eigen::MatrixXd( ) > parameterPartialFunctionToAdd );

std::pair< boost::function< Eigen::MatrixXd( ) >, int > createMergedParameterPartialFunction(
        const std::pair< boost::function< Eigen::MatrixXd( ) >, int >& partialFunctionOfAccelerationToAdd,
        const std::pair< boost::function< Eigen::MatrixXd( ) >, int >& partialFunctionOfAccelerationToSubtract );

boost::function< Eigen::MatrixXd( ) > getCombinedCurrentDoubleParameterFunction(
        const boost::shared_ptr< StateDerivativePartial > firstPartial,
        const boost::shared_ptr< StateDerivativePartial > secondPartial,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterObject,
        const int firstPartialSize, const int secondPartialSize,
        const bool subtractPartials = 0 );

boost::function< Eigen::MatrixXd( ) > getCombinedCurrentVectorParameterFunction(
        const boost::shared_ptr< StateDerivativePartial > firstPartial,
        const boost::shared_ptr< StateDerivativePartial > secondPartial,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterObject,
        const int firstPartialSize, const int secondPartialSize,
        const bool subtractPartials = 0 );

}

}

}

#endif // TUDAT_STATEDERIVATIVEPARTIAL_H
