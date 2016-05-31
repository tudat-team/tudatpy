#ifndef ENTRYGUIDANCE_H
#define ENTRYGUIDANCE_H

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"

namespace tudat
{

namespace aerodynamics
{

class EntryGuidance
{
public:
    EntryGuidance( ){ }

    virtual ~EntryGuidance( ){ }

    virtual void updateGuidance( const double currentTime ) = 0;

    double getCurrentAngleOfAttack( )
    {
        return currentAngleOfAttack_;
    }

    double getCurrentAngleOfSideslip( )
    {
        return currentAngleOfSideslip_;
    }

    double getCurrentBankAngle( )
    {
        return currentBankAngle_;
    }

protected:

    double currentAngleOfAttack_;

    double currentAngleOfSideslip_;

    double currentBankAngle_;

};

//class DerivedClassEntryGuidance: public EntryGuidance
//{
//public:
//    DerivedClassEntryGuidance( ){ }

//    ~DerivedClassEntryGuidance( ){ }

//    void updateGuidance( const double currentTime )
//    {
//        currentAngleOfAttack_ = ... performCalculationsHere...
//                currentAngleOfSideslip_ = ... performCalculationsHere...
//                currentBankAngle_ = ... performCalculationsHere...

//    }



//};


void setGuidanceAnglesFunctions(
        const boost::shared_ptr< EntryGuidance > entryGuidance,
        const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator );

}

}

#endif // ENTRYGUIDANCE_H
