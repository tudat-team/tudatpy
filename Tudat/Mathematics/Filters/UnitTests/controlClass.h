/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_CONTROL_CLASS_H
#define TUDAT_CONTROL_CLASS_H

#include <Eigen/Core>
#include <boost/function.hpp>

namespace tudat
{

namespace filters
{

//! Class for control vector.
template< typename IndependentVariableType, typename DependentVariableType, int NumberOfElements >
class ControlWrapper
{
public:

    //! Typedef of the control vector.
    typedef Eigen::Matrix< DependentVariableType, NumberOfElements, 1 > DependentVector;

    //! Typedef of the function describing the system.
    typedef std::function< DependentVector( const IndependentVariableType,
                                            const DependentVector& ) > ControlFunction;

    //! Default constructor.
    /*!
     *  Default constructor.
     *  \param controlFunction Function to compute the control vector.
     */
    ControlWrapper( const ControlFunction& controlFunction ) :
        controlFunction_( controlFunction )
    {
        // Set control vector to zero
        currentControlVector_.setZero( );
    }

    //! Default destructor.
    ~ControlWrapper( ) { }

    //! Function to retireve the current control vector.
    /*!
     *  Function to retireve the current control vector. The function setCurrentControlVector needs to be called before the control
     *  vector is retrieved.
     *  \return Vector denoting the current control.
     */
    DependentVector getCurrentControlVector( )
    {
        return currentControlVector_;
    }

    //! Function to set the control vector.
    /*!
     *  Function to set the control vector, based on the current time and state. The control vector is then
     *  computed based on the input control function.
     *  \param currentTime Double denoting the current time.
     *  \param currentStateVector Vector denoting the current state.
     */
    void setCurrentControlVector( const IndependentVariableType currentTime,
                                  const DependentVector& currentStateVector )
    {
        currentControlVector_ = controlFunction_( currentTime, currentStateVector );
    }

private:

    //! Function to compute the control vector.
    ControlFunction controlFunction_;

    //! Vector denoting the current control.
    DependentVector currentControlVector_;

};

} // namespace filters

} // namespace tudat

#endif // TUDAT_CONTROL_CLASS_H
