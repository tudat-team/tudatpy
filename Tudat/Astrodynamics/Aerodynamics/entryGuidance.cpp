#include "Tudat/Astrodynamics/Aerodynamics/entryGuidance.h"

namespace tudat
{

namespace aerodynamics
{

//! Function that must be called to link the EntryGuidance object to the simulation
void setGuidanceAnglesFunctions(
        const boost::shared_ptr< EntryGuidance > entryGuidance,
        const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator )
{
    angleCalculator->setOrientationAngleFunctions(
                boost::bind( &EntryGuidance::getCurrentAngleOfAttack, entryGuidance ),
                boost::bind( &EntryGuidance::getCurrentAngleOfSideslip, entryGuidance ),
                boost::bind( &EntryGuidance::getCurrentBankAngle, entryGuidance ),
                boost::bind( &EntryGuidance::updateGuidance, entryGuidance,_1 ) );
}

}

}
