#ifndef FIRSTORDERRELATIVISTICLIGHTTIMECORRECTION_H
#define FIRSTORDERRELATIVISTICLIGHTTIMECORRECTION_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include <cmath>
#include <vector>
#include <iostream>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "Astrodynamics/ObservationModels/ObservableCorrections/lightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

//! Class to calculate first order relativistic light time correction due to a set of gravitating point masses.
/*!
 *  Class to calculate first order relativistic light time correction due to a gravitating point mass. Has the properties (mass and position)
 *  of the gravitating bodies as member variables; receives the states and times of receiver and transmitter as input to
 *  calculateCumulativeFirstOrderCorrections function each time a light time is to be calculated.
 */
class FirstOrderLightTimeCorrectionCalculator: public LightTimeCorrection
{
public:

    //! Constructor, takes and sets gravitating body properties.
    /*!
     *  Constructor, takes and sets gravitating body properties.
     *  \param perturbingBodyStateFunctions Set of function returning the state of the gravitating bodies as a function of time.
     *  \param perturbingBodyGravitationalParameterFunctions Set of functions returning the gravitational parameters of the gravitating bodies.
     *  \param perturbingBodyStateFunctions Function returning the parametric post-Newtonian parameter gamma, a measure for the
     *  space-time curvature due to a unit rest mass (1.0 in GR)
     */
    FirstOrderLightTimeCorrectionCalculator(
            const std::vector< boost::function< basic_mathematics::Vector6d( const double ) > >& perturbingBodyStateFunctions,
            const std::vector< boost::function< double( ) > >& perturbingBodyGravitationalParameterFunctions,
            const std::vector< std::string > perturbingBodyNames,
            const boost::function< double( ) >& ppnParameterGammaFunction = boost::lambda::constant( 1.0 ) ):
        LightTimeCorrection( first_order_relativistic ),
        perturbingBodyStateFunctions_( perturbingBodyStateFunctions ),
        perturbingBodyGravitationalParameterFunctions_( perturbingBodyGravitationalParameterFunctions ),
        perturbingBodyNames_( perturbingBodyNames ),
        ppnParameterGammaFunction_( ppnParameterGammaFunction )
    {
        currentTotalLightTimeCorrection_ = 0.0;
        currentLighTimeCorrectionComponents_.resize( perturbingBodyNames_.size( ) );
    }

    ~FirstOrderLightTimeCorrectionCalculator( ){ }

    //! Function to calculate first order relativistic light time correction due to set of gravitating point masses.
    /*!
     *  Function to calculate first order relativistic light time correction due to set of gravitating point masses,
     *  according to Eq. (11.17) of 2010 IERS conventions. Calculation are performed by calling calculateFirstOrderLightTimeCorrectionFromCentralBody
     *  function for each gravitating body.
     *  \param transmitterState State of transmitter at transmission time.
     *  \param receiverState State of receiver at reception time
     *  \param transmissionTime Time of signal transmission
     *  \param receptionTime Time of signal reception
     *  \return Total light time correction due to gravitating masses defined by perturbingBodyStateFunctions_ and
     *  perturbingBodyGravitationalParameterFunctions_ member variables.
     */
    double calculateLightTimeCorrection( const basic_mathematics::Vector6d& transmitterState,
                                         const basic_mathematics::Vector6d& receiverState,
                                         const double transmissionTime,
                                         const double receptionTime );

    std::vector< std::string > getPerturbingBodyNames( )
    {
        return perturbingBodyNames_;
    }

    std::vector< boost::function< double( ) > > getPerturbingBodyGravitationalParameterFunctions( )
    {
        return perturbingBodyGravitationalParameterFunctions_;
    }

    double getCurrentTotalLightTimeCorrection( )
    {
        return currentTotalLightTimeCorrection_;
    }

    double getCurrentLightTimeCorrectionComponent( const int bodyIndex )
    {
        return currentLighTimeCorrectionComponents_.at( bodyIndex );
    }

    boost::function< double( ) > getPpnParameterGammaFunction_( )
    {
        return ppnParameterGammaFunction_;
    }

private:
    //! Set of function returning the state of the gravitating bodies as a function of time.
    /*!
     *  Set of function returning the state of the gravitating bodies as a function of time.
     */
    std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > perturbingBodyStateFunctions_;

    //! Set of functions returning the gravitational parameters of the gravitating bodies.
    /*!
     *  Set of functions returning the gravitational parameters of the gravitating bodies.
     */
    std::vector< boost::function< double( ) > > perturbingBodyGravitationalParameterFunctions_;

    //! Function returning the parametric post-Newtonian parameter gamma
    /*!
     *  Function returning the parametric post-Newtonian parameter gamma, a measure for the space-time curvature due to a unit rest mass (1.0 in GR)
     */
    boost::function< double( ) > ppnParameterGammaFunction_;

    std::vector< std::string > perturbingBodyNames_;

    std::vector< double > currentLighTimeCorrectionComponents_;

    double currentTotalLightTimeCorrection_;
};

class FirstOrderLightTimeCorrectionDerivativeCalculator: public LightTimeDerivativeCorrection
{
public:


    FirstOrderLightTimeCorrectionDerivativeCalculator(
            const std::vector< boost::function< basic_mathematics::Vector6d( const double ) > >& perturbingBodyStateFunctions,
            const std::vector< boost::function< double( ) > >& perturbingBodyGravitationalParameterFunctions,
            const std::vector< std::string > perturbingBodyNames,
            const boost::function< double( ) >& ppnParameterGammaFunction = boost::lambda::constant( 1.0 ) ):
        LightTimeDerivativeCorrection( first_order_relativistic ),
        perturbingBodyStateFunctions_( perturbingBodyStateFunctions ),
        perturbingBodyGravitationalParameterFunctions_( perturbingBodyGravitationalParameterFunctions ),
        perturbingBodyNames_( perturbingBodyNames ),
        ppnParameterGammaFunction_( ppnParameterGammaFunction )
    {
        currentTotalLightTimeDerivativeCorrection_ = 0.0;
        currentLighTimeCorrectionDerivativeComponents_.resize( perturbingBodyNames_.size( ) );
    }

    ~FirstOrderLightTimeCorrectionDerivativeCalculator( ){ }


    double calculateLightTimeDerivativeCorrection( const basic_mathematics::Vector6d& transmitterState,
                                                   const basic_mathematics::Vector6d& receiverState,
                                                   const double transmissionTime,
                                                   const double receptionTime,
                                                   const bool fixTransmissionTime );

    std::vector< std::string > getPerturbingBodyNames( )
    {
        return perturbingBodyNames_;
    }

    std::vector< boost::function< double( ) > > getPerturbingBodyGravitationalParameterFunctions( )
    {
        return perturbingBodyGravitationalParameterFunctions_;
    }

    double getCurrentTotalLightTimeCorrection( )
    {
        return currentTotalLightTimeDerivativeCorrection_;
    }

    double getCurrentLightTimeCorrectionComponent( const int bodyIndex )
    {
        return currentLighTimeCorrectionDerivativeComponents_.at( bodyIndex );
    }

    boost::function< double( ) > getPpnParameterGammaFunction_( )
    {
        return ppnParameterGammaFunction_;
    }

private:

    std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > perturbingBodyStateFunctions_;


    std::vector< boost::function< double( ) > > perturbingBodyGravitationalParameterFunctions_;


    boost::function< double( ) > ppnParameterGammaFunction_;

    std::vector< std::string > perturbingBodyNames_;

    std::vector< double > currentLighTimeCorrectionDerivativeComponents_;

    double currentTotalLightTimeDerivativeCorrection_;
};


}

}
#endif // FIRSTORDERRELATIVISTICLIGHTTIMECORRECTION_H
