#ifndef TORQUEMODEL_H
#define TORQUEMODEL_H

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace basic_astrodynamics
{

class TorqueModel
{
public:

    TorqueModel( ):
        currentTime_( TUDAT_NAN ){ }

virtual ~TorqueModel( ) { }

virtual Eigen::Vector3d getTorque( ) = 0;

virtual void updateMembers( const double currentTime ) = 0;

    //! Function to reset the current time
    /*!
     * Function to reset the current time of the acceleration model.
     * \param currentTime Current time (default NaN).
     */
    virtual void resetTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
    }

protected:
    double currentTime_;
private:

};


typedef std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::TorqueModel > > > SingleBodyTorqueModelMap;

typedef std::map< std::string, SingleBodyTorqueModelMap > TorqueModelMap;

}

}

#endif // TORQUEMODEL_H
