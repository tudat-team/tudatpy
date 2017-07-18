#ifndef TIDALLOVENUMBERPARTIALINTERFACE_H
#define TIDALLOVENUMBERPARTIALINTERFACE_H

#include <boost/math/special_functions/factorials.hpp>

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"

#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicPartialFunctions.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/tidalLoveNumber.h"

namespace tudat
{

namespace orbit_determination
{

void getMaximumUsedDegreeAndOrder(
        const int maximumDegree, const int maximumOrder, const int evaluationDegree,
        int& maximumUsedDegree, int& maximumUsedOrder );

//! (Base) class for calculating the partial of a spherical harmonic acceleration w.r.t. a tidal love number.
/*!
 *  (Base) class for calculating the partial of a spherical harmonic acceleration w.r.t. a tidal love number. These calculations are
 *  implemented separately from the SphericalHarmonicsGravityPartial class as its functionality is quite specific and only used by
 *  the SphericalHarmonicsGravityPartial for certain cases. Also, this architecture allows different implementations for
 *  different tidal models. This base class assumes no direct tidal lag in teh calculations, i.e. state and rotation functions are all
 *  evaluated at same time and lag is implemented by the possibility of complex love numbers (for non-zonal terms).
 */
class TidalLoveNumberPartialInterface
{
public:

    //! Constructor from objects.
    /*!
     *  Constructor from objects.
     *  \param gravityFieldVariations Gravity field variation object which calculates spherical harmonic coefficient variations
     *  due to considered love number(s).
     *  \param accelerationModel Spherical harmonic acceleration for which the partial is to be calculated.
     *  \param deformingBodyStateFunctions Position functions of bodies causing deformation. These should be taken directly from gravityFieldVariations.
     */
    TidalLoveNumberPartialInterface(
            const boost::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariations,
            const boost::function< Eigen::Vector3d( ) > deformedBodyPositionFunction,
            const std::vector< boost::function< Eigen::Vector3d( ) > > deformingBodyStateFunctions,
            const boost::function< Eigen::Quaterniond( ) > rotationToDeformedBodyFrameFrameFunction,
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

        realLoveNumberScaler_ = std::make_pair( ( Eigen::Vector2d( )<< 1.0, 0.0 ).finished( ), ( Eigen::Vector2d( )<< 0.0, 1.0 ).finished( ) );
        complexLoveNumberScaler_ = std::make_pair( ( Eigen::Matrix2d( )<< 1.0, 0.0, 0.0, -1.0  ).finished( ),
                                                   ( Eigen::Matrix2d( )<< 0.0, 1.0, 1.0, 0.0 ).finished( ) );

        for( unsigned int i = 0; i < deformingBodyStateFunctions_.size( ); i++ )
        {
            allDeformingBodyIndices_.push_back( i );
        }

    }

    //! Destructor
    /*!
     *  Destructor
     */
    virtual ~TidalLoveNumberPartialInterface( ){ }

    std::vector< int > getSelectedDeformingBodyIds( const std::vector< std::string >& selectedBodyNames );

    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > calculateShericalHarmonicCoefficientPartialMatrix(
            const std::vector< Eigen::Vector2d >& coefficientPartialsPerOrder,
            const std::vector< int >& orders,
            const std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, Eigen::Matrix< double, 2, Eigen::Dynamic > >& coefficientPartialScalers );

    //! Function to calculate the partial of sh acceleration wrt complex tidal love numbers.
    /*!
     *  Function to calculate the partial of sh acceleration wrt complex tidal love numbers at a single degree. The orders of
     *  the love numbers wrt which the partials to be taken are required as input. The output can be either summed over all orders,
     *  or per order.
     *  \param degree Degree of love numbers wrt which partials are to be taken
     *  \param orders SH orders in current degree at which the partials are to be taken.
     */
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > calculateShericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumbers(
            const int degree,
            const std::vector< int >& orders,
            const std::vector< int >& deformingBodyIndices,
            const int maximumDegree,
            const int maximumOrder );

    //! Function to calculate the partial of sh acceleration wrt complex tidal love numbers.
    /*!
     *  Function to calculate the partial of sh acceleration wrt complex tidal love numbers at a single degree, summed over all orders in degree.
     *  \param degree Degree of love numbers wrt which partials are to be taken
     *  \return Partial of sh acceleration wrt complex love number, two columns (first real, then complex part).
     */
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > calculateShericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumber(
            const int degree,
            const std::vector< int >& deformingBodyIndices,
            const int maximumDegree,
            const int maximumOrder )
    {
        return calculateShericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumbers(
                        degree, fullDegreeOrders[ degree - 2 ], deformingBodyIndices, maximumDegree, maximumOrder );
    }

    //! Function to calculate the partial of sh acceleration wrt real tidal love numbers.
    /*!
     *  Function to calculate the partial of sh acceleration wrt real tidal love numbers at a single degree. The orders of
     *  the love numbers wrt which the partials to be taken are required as input. The output can be either summed over all orders,
     *  or per order.
     *  \param degree Degree of love numbers wrt which partials are to be taken
     *  \param orders SH orders in current degree at which the partials are to be taken.
     */
    virtual std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > calculateShericalHarmonicCoefficientsPartialWrtRealTidalLoveNumbers(
            const int degree,
            const std::vector< int >& orders,
            const std::vector< int >& deformingBodyIndices,
            const int maximumDegree,
            const int maximumOrder );

    //! Function to calculate the partial of sh acceleration wrt real tidal love numbers.
    /*!
     *  Function to calculate the partial of sh acceleration wrt real tidal love numbers at a single degree, summed over all orders in degree.
     *  \param degree Degree of love numbers wrt which partials are to be taken
     *  \return Partial of sh acceleration wrt real love number, two columns..
     */
    virtual std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > calculateShericalHarmonicCoefficientsPartialWrtRealTidalLoveNumber(
            const int degree,
            const std::vector< int >& deformingBodyIndices,
            const int maximumDegree,
            const int maximumOrder )
    {
        return calculateShericalHarmonicCoefficientsPartialWrtRealTidalLoveNumbers(
                    degree, fullDegreeOrders[ degree - 2 ], deformingBodyIndices, maximumDegree, maximumOrder );
    }

    void update( const double currentTime )
    {
        currentTime_ = currentTime;
        rotationToTidallyDeformedBody_ = rotationToDeformedBodyFrameFrameFunction_( );

    }

    void resetTime( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime  ) )
        {
            currentDoubleParameterPartials_.clear( );
            currentVectorParameterPartials_.clear( );
        }
        currentTime_ = currentTime;
    }

    void updateParameterPartials( )
    {
        for( parameterDoublePartialFunctionIterator_ = parameterDoublePartialFunctions_.begin( );
             parameterDoublePartialFunctionIterator_ != parameterDoublePartialFunctions_.end( );
             parameterDoublePartialFunctionIterator_++ )
        {
            if( currentDoubleParameterPartials_.count( parameterDoublePartialFunctionIterator_->first ) == 0 )
            {
                currentDoubleParameterPartials_[ parameterDoublePartialFunctionIterator_->first ] = parameterDoublePartialFunctionIterator_->second( );
            }
        }

        for( parameterVectorPartialFunctionIterator_ = parameterVectorPartialFunctions_.begin( );
             parameterVectorPartialFunctionIterator_ != parameterVectorPartialFunctions_.end( );
             parameterVectorPartialFunctionIterator_++ )
        {
            //std::cout<<"New vec. partial: "<<parameterVectorPartialFunctionIterator_->first->getParameterName( ).first<<std::endl;
            if( currentVectorParameterPartials_.count( parameterVectorPartialFunctionIterator_->first ) == 0 )
            {
                currentVectorParameterPartials_[ parameterVectorPartialFunctionIterator_->first ] = parameterVectorPartialFunctionIterator_->second( );
            }
        }
    }

    virtual std::pair< int, std::pair< int, int > > setParameterPartialFunction(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
            const int maximumDegree,
            const int maximumOrder );

    virtual std::pair< int, std::pair< int, int > > setParameterPartialFunction(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
            const int maximumDegree,
            const int maximumOrder );

    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >& getCurrentParameterPartial(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
            const std::pair< int, int > maximumDegreeAndOrder )
    {
        if( currentDoubleParameterPartials_.count( std::make_pair( parameter, maximumDegreeAndOrder ) ) == 0 )
        {
            if( parameterDoublePartialFunctions_.count( std::make_pair( parameter, maximumDegreeAndOrder ) ) == 0 )
            {
                std::cerr<<"Parameter of type "<<parameter->getParameterName( ).first<<", "<<
                           parameter->getParameterName( ).second.first<<", "<<
                           parameter->getParameterName( ).second.second<<" not found in list of existing partials"<<std::endl;
            }
            else
            {
                std::cerr<<"Warning, double partial should already be calculated"<<std::endl;
                currentDoubleParameterPartials_[ std::make_pair( parameter, maximumDegreeAndOrder ) ] =
                        parameterDoublePartialFunctions_.at( std::make_pair( parameter, maximumDegreeAndOrder ) )( );
            }
        }
        return currentDoubleParameterPartials_.at( std::make_pair( parameter, maximumDegreeAndOrder ) );
    }

    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >& getCurrentDoubleParameterPartial(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
            const std::pair< int, int > maximumDegreeAndOrder )
    {
        return getCurrentParameterPartial( parameter, maximumDegreeAndOrder );
    }

    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >& getCurrentParameterPartial(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
            const std::pair< int, int > maximumDegreeAndOrder  )
    {

        if( currentVectorParameterPartials_.count( std::make_pair( parameter, maximumDegreeAndOrder ) ) == 0 )
        {
            if( parameterVectorPartialFunctions_.count( std::make_pair( parameter, maximumDegreeAndOrder ) ) == 0 )
            {
                std::cerr<<"Parameter of type "<<parameter->getParameterName( ).first<<", "<<
                           parameter->getParameterName( ).second.first<<", "<<
                           parameter->getParameterName( ).second.second<<" not found in list of existing partials"<<std::endl;
            }
            else
            {
                std::cerr<<"Warning, vector partial should already be calculated in Love number interface"<<std::endl;
                currentVectorParameterPartials_[ std::make_pair( parameter, maximumDegreeAndOrder ) ] =
                        parameterVectorPartialFunctions_.at( std::make_pair( parameter, maximumDegreeAndOrder ) )( );
            }
        }

        return currentVectorParameterPartials_.at( std::make_pair( parameter, maximumDegreeAndOrder ) );
    }

    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >& getCurrentVectorParameterPartial(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
            const std::pair< int, int > maximumDegreeAndOrder )
    {
        return getCurrentParameterPartial( parameter, maximumDegreeAndOrder );
    }



protected:

    //! Function to calculate the partial of sh coefficients wrt real part of tidal love numbers.
    /*!
     *  Function to calculate the partial of sh coefficients wrt real part of tidal love numbers. Partials at a single degree are calculated,
     *  with the desired orders required as input.
     *  \param degree SH degree at which partials are to be taken.
     *  \param orders SH orders at given degree at which partials are to be taken (if vector is empty, all orders at given degree are used,
     *  in ascending order).
     *  \return Vector of partials of SH coefficients wrt real part of tidal love number ([C_{nm};S_{nm}]). Entries of return vector correspond
     *  to orders given by corresponding entry in orders vector.
     */
    std::vector< Eigen::Vector2d > calculateCoefficientPartialWrtRealTidalLoveNumber(
            const int degree,
            const std::vector< int >& orders,
            const std::vector< int >& deformingBodyIndices,
            const int maximumDegree,
            const int maximumOrder );

    //! Function to pre-calculate all states of bodies (+ rotation to local frame) at given degree and set of orders.
    /*!
     *  Function to pre-calculate all states of bodies (+ rotation to local frame) at given degree and set of orders.
     *  Implemented here assuming no time lag, derived classes can provide modified implementations.
     *  \param degree Degree of deformation for calculations of which calculated states will be used
     *  \param orders SH orders in current degree of deformation for calculations of which calculated states will be used.
     */
    virtual void updateCurrentTidalBodyStates( const int degree,
                                               const std::vector< int >& orders,
                                               const std::vector< int >& deformingBodiesToUpdate );

    //! Function to set current states of bodies (+ rotation of local frame) and derived quantities
    /*!
     *  Function to set current states of bodies (+ rotation of local frame) and derived quantities for a given body and order.
     *  Sets the 'Variables used in calculation loop' from pre-calculations of updateCurrentTidalBodyStates function.
     *  \param order SH order of calculations to use, corresponding to entry of orders varaiable passed to updateCurrentTidalBodyStates.
     *  \param body Index of body causing deformation for which to set 'Variables used in calculation loop'.
     */
    virtual void setCurrentTidalBodyStates( const int degree, const int order, const int body );

    //! Reference radius used by tidal deformation model.
    /*!
     *  Reference radius used by tidal deformation model.
     */
    double deformedBodyReferenceRadius_;



    // Functions for retrieving states at current time, for acceleration derivative.
    //! Function to retrieve current state of body causing acceleration (= body being tidally deformed).
    /*!
     *  Function to retrieve current state of body causing acceleration (= body being tidally deformed). Typically linked to
     *  SH acceleration model. Function is also used for tidal sh variation calculations in case of no tidal time delay.
     */
    boost::function< Eigen::Vector3d( ) > deformedBodyPositionFunction_;

    //! Function to retrieve current state of body being accelerated.
    /*!
     *  Function to retrieve current state of body being accelerated.
     */
    boost::function< Eigen::Vector3d( ) > acceleratingBodyPositionFunction_;



    //! Function to retrieve current gravitational parameter of body being deformed.
    /*!
     *  Function to retrieve current gravitational parameter of body being deformed.
     */
    boost::function< double( ) > deformedBodyGravitationalParameterFunction_;


    //! Current gravitational parameter of body being deformed.
    /*!
     *  Current gravitational parameter of body being deformed.
     */
    double deformedBodyGravitationalParameter_;

    //! Vector of function to retrieve current gravitational parameters of bodies causing deformation.
    /*!
     *  Vector of function to retrieve current gravitational parameters of bodies causing deformation.
     */
    std::vector< boost::function< double( ) > > deformingBodyGravitationalParameters_;

    // Functions for retrieving states at current time of bodies causing deformation.
    std::vector< boost::function< Eigen::Vector3d( ) > > deformingBodyStateFunctions_;

    boost::function< Eigen::Quaterniond( ) > rotationToDeformedBodyFrameFrameFunction_;

    double currentTime_;



    // Variables used in calculation loop
    Eigen::Vector3d relativeDeformingBodyPosition_;

    Eigen::Vector3d relativeDeformingBodySphericalPosition_;

    double massRatio_, radiusRatio_, sineOfLatitude_;

    std::complex< double > iLongitude_;



    std::pair< Eigen::Matrix< double, 2, 1 >, Eigen::Matrix< double, 2, 1 > > realLoveNumberScaler_;

    std::pair< Eigen::Matrix< double, 2, 2 >, Eigen::Matrix< double, 2, 2 > > complexLoveNumberScaler_;

    std::vector< int > allDeformingBodyIndices_;

    std::string deformedBody_;

    std::vector< std::string > deformingBodies_;


    // Positions + rotation of deforming/deformed bodies.
    Eigen::Quaterniond rotationToTidallyDeformedBody_;

    Eigen::Vector3d positionOfTidallyDeformedBody_;

    std::vector< Eigen::Vector3d > positionsOfDeformingBodies_;


    std::map< std::pair< boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > >,
    std::pair< int, int > >, std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > > currentDoubleParameterPartials_;

    std::map< std::pair< boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > >,
    std::pair< int, int > >, boost::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > > parameterDoublePartialFunctions_;

    std::map< std::pair< boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > >,
    std::pair< int, int > >, boost::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >::iterator parameterDoublePartialFunctionIterator_;


    std::map< std::pair< boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
    std::pair< int, int > >, std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > > currentVectorParameterPartials_;

    std::map< std::pair< boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
    std::pair< int, int > >, boost::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > > parameterVectorPartialFunctions_;

    std::map< std::pair< boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
    std::pair< int, int > >, boost::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >::iterator parameterVectorPartialFunctionIterator_;
};

}

}


#endif // TIDALLOVENUMBERPARTIALINTERFACE_H
