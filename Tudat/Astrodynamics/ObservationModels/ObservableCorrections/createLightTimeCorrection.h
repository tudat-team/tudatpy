#ifndef LIGHTTIMECORRECTIONFUNCTION_H
#define LIGHTTIMECORRECTIONFUNCTION_H

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Astrodynamics/ObservationModels/ObservableCorrections/lightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

//! Typedef for function calculating light time correction in light time calculation loop.
typedef boost::function< double( const basic_mathematics::Vector6d, const basic_mathematics::Vector6d,
                                 const double, const double ) > LightTimeCorrectionFunction;

//! Base class for light-time correction settings.
/*!
 *  Base class for light-time correction settings. This class is not used for calculations of corrections,
 *  but is used for input purposes. The createLightTimeCorrections function produces the functions
 *  that calculate teh actual corrections.
 */
class LightTimeCorrectionSettings
{
public:

    //!  Constructor, takes light-time correction type.
    /*!
     *   Constructor, takes light-time correction type.
     */
    LightTimeCorrectionSettings( const LightTimeCorrectionType correctionType ):
        correctionType_( correctionType ){ }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    virtual ~LightTimeCorrectionSettings( ){ }

    //! Function returning the correction type.
    /*!
     *  Function returning the correction type.
     */
    LightTimeCorrectionType getCorrectionType( ){ return correctionType_; }

protected:

    //! Correction type.
    /*!
     *  Correction type.
     */
    LightTimeCorrectionType correctionType_;
};

typedef std::map< LinkEnds, std::vector< boost::shared_ptr< LightTimeCorrectionSettings > > > LightTimeCorrectionSettingsMap;

class TroposphericCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    TroposphericCorrectionSettings( const LightTimeCorrectionType correctionType,
                                    const std::string& bodyWithAtmosphere  ):
        LightTimeCorrectionSettings( correctionType ), bodyWithAtmosphere_( bodyWithAtmosphere ){ }

    std::string getBodyWithAtmosphere( ){ return bodyWithAtmosphere_; }

private:
    std::string bodyWithAtmosphere_;

};

class IonosphericCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    IonosphericCorrectionSettings( const LightTimeCorrectionType correctionType,
                                    const std::string& bodyWithAtmosphere ):
        LightTimeCorrectionSettings( correctionType ), bodyWithAtmosphere_( bodyWithAtmosphere ){ }

    std::string getBodyWithAtmosphere( ){ return bodyWithAtmosphere_; }

private:
    std::string bodyWithAtmosphere_;

};


class FirstOrderRelativisticLightTimeCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    FirstOrderRelativisticLightTimeCorrectionSettings( const std::vector< std::string >& perturbingBodies ):
        LightTimeCorrectionSettings( first_order_relativistic ), perturbingBodies_( perturbingBodies ){ }

    std::vector< std::string > getPerturbingBodies( ){ return perturbingBodies_; }

private:
    std::vector< std::string > perturbingBodies_;

};


class SecondOrderRelativisticLightTimeCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    SecondOrderRelativisticLightTimeCorrectionSettings( const std::vector< std::string >& perturbingBodies ):
        LightTimeCorrectionSettings( second_order_relativistic ), perturbingBodies_( perturbingBodies ){ }

    std::vector< std::string > getPerturbingBodies( ){ return perturbingBodies_; }

private:
    std::vector< std::string > perturbingBodies_;

};

class J2RelativisticLightTimeCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    J2RelativisticLightTimeCorrectionSettings( const std::vector< std::string >& perturbingBodies ):
        LightTimeCorrectionSettings( j2_relativistic ), perturbingBodies_( perturbingBodies ){ }

    std::vector< std::string > getPerturbingBodies( ){ return perturbingBodies_; }

private:
    std::vector< std::string > perturbingBodies_;

};

boost::shared_ptr< LightTimeCorrection > createLightTimeCorrections(
        const boost::shared_ptr< LightTimeCorrectionSettings >& correctionSettings,
        const NamedBodyMap& bodyMap,
        const std::pair< std::string, std::string >& transmitter,
        const std::pair< std::string, std::string >& receiver );

boost::shared_ptr< LightTimeDerivativeCorrection > createLightTimeDerivativeCorrections(
        const boost::shared_ptr< LightTimeCorrectionSettings > correctionSettings,
        const NamedBodyMap& bodyMap,
        const std::pair< std::string, std::string >& transmitter,
        const std::pair< std::string, std::string >& receiver );
}

}


#endif // LIGHTTIMECORRECTIONFUNCTION_H
