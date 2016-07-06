/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATELIGHTTIMECORRECTION_H
#define TUDAT_CREATELIGHTTIMECORRECTION_H

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/lightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

//! Typedef for function calculating light-time correction in light-time calculation loop.
typedef boost::function< double(
        const basic_mathematics::Vector6d&, const basic_mathematics::Vector6d&,
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


class FirstOrderRelativisticLightTimeCorrectionSettings: public LightTimeCorrectionSettings
{
public:
    FirstOrderRelativisticLightTimeCorrectionSettings( const std::vector< std::string >& perturbingBodies ):
        LightTimeCorrectionSettings( first_order_relativistic ), perturbingBodies_( perturbingBodies ){ }

    ~FirstOrderRelativisticLightTimeCorrectionSettings( ){ }

    std::vector< std::string > getPerturbingBodies( ){ return perturbingBodies_; }

private:
    std::vector< std::string > perturbingBodies_;

};


boost::shared_ptr< LightTimeCorrection > createLightTimeCorrections(
        const boost::shared_ptr< LightTimeCorrectionSettings > correctionSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::pair< std::string, std::string >& transmitter,
        const std::pair< std::string, std::string >& receiver );

}

}


#endif // TUDAT_CREATELIGHTTIMECORRECTION_H
