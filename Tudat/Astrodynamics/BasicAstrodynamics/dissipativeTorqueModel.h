#ifndef TUDAT_DISSIPATIVETORQUEMODEL_H
#define TUDAT_DISSIPATIVETORQUEMODEL_H

#include <iomanip>
#include <functional>

#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Class to compute a phenomenological dissipative torque model
/*!
 *  Class to compute a phenomenological dissipative torque model, with the torque proportional to the angular velocity deviation
 *  w.r.t. the mean rotation (about body-fixed z-axis).
 */
class DissipativeTorqueModel: public TorqueModel
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param bodyFixedRotationVectorFunction Function that returns body-fixed angular velocity vector of body
     * \param dampingMatrixFunction Function that returns the damping matrix that defines the dissipation strength
     * \param bodyMeanRotationRate Mean angular velocity vector about z-axis (w.r.t. which the deviation is computed)
     */
    DissipativeTorqueModel(
            const std::function< Eigen::Vector3d( ) >& bodyFixedRotationVectorFunction,
            const std::function< Eigen::Matrix3d( ) >& dampingMatrixFunction,
            const double bodyMeanRotationRate ):
        bodyFixedRotationVectorFunction_( bodyFixedRotationVectorFunction ),
        dampingMatrixFunction_( dampingMatrixFunction ),
    bodyMeanRotationRate_( bodyMeanRotationRate ){ }

    //! Destructor
    ~DissipativeTorqueModel( ) { }

    //! Get dissipative torque.
    /*!
     * Returns the dissipative torque.
     * \return Dissipative torque.
     */
    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }


    //! Update member variables used by the torque model.
    /*!
     * Updates member variables used by the torque model.
     * Function pointers to retrieve the current values of quantities from which the
     * torque is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current state.
     * \param currentTime Time at which torque model is to be updated.
     */
    void updateMembers( const double currentTime )
    {
        currentBodyRotationPerturbationVector_ = bodyFixedRotationVectorFunction_( );
        currentBodyRotationPerturbationVector_( 2 ) -= bodyMeanRotationRate_;
        currentTorque_ =  -dampingMatrixFunction_( ) * currentBodyRotationPerturbationVector_;
    }

    //! Function to modify the damping matrix
    /*!
     * Function to modify the damping matrix
     * \param dampingMatrix New damping matrix
     */
    void setDampingMatrixFunction( const Eigen::Matrix3d& dampingMatrix )
    {
        dampingMatrixFunction_ = [ = ]( ){ return dampingMatrix; };
    }

protected:

private:

    //! tion that returns body-fixed angular velocity vector of body
    std::function< Eigen::Vector3d( ) > bodyFixedRotationVectorFunction_;

    //! Function that returns the damping matrix that defines the dissipation strength
    std::function< Eigen::Matrix3d( ) > dampingMatrixFunction_;

    //! Mean angular velocity vector about z-axis (w.r.t. which the deviation is computed)
    double bodyMeanRotationRate_;

    //! Current deviation  of angular velocity vector from nominal rotation
    Eigen::Vector3d currentBodyRotationPerturbationVector_;

    //! Current torque, as computed by last call to updateMembers function
    Eigen::Vector3d currentTorque_;

};

}

}

#endif // TUDAT_DISSIPATIVETORQUEMODEL_H
