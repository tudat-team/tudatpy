/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_BASESTATEINTERFACE_H
#define TUDAT_BASESTATEINTERFACE_H

#include <functional>
#include <string>

#include <Eigen/Core>

#include "tudat/basics/timeType.h"

namespace tudat
{

namespace simulation_setup
{

//! Base class used for the determination of the inertial state of a Body's ephemeris origin
/*!
 *  Base class used for the determination of the inertial state of a Body's ephemeris origin. This base class is used
 *  to provide an untemplated interface class through which to call the base frame state. The state may be defined
 *  in a templated manner in the derived class.
 */
class BaseStateInterface
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param baseFrameId Name of frame origin for which inertial state is computed by this class
     */
    BaseStateInterface( const std::string baseFrameId ): baseFrameId_( baseFrameId ) { }

    //! Destructor
    virtual ~BaseStateInterface( ) { }

    //! Function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Function through which the state of baseFrameId_ in the inertial frame can be determined
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    template< typename OutputTimeType, typename OutputStateScalarType >
    Eigen::Matrix< OutputStateScalarType, 6, 1 > getBaseFrameState( const OutputTimeType time );

    std::string getBaseFrameId( )
    {
        return baseFrameId_;
    }

protected:
    //! Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (double time and double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    virtual Eigen::Matrix< double, 6, 1 > getBaseFrameDoubleState( const double time ) = 0;

    //! Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (double time and double long state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    virtual Eigen::Matrix< long double, 6, 1 > getBaseFrameLongDoubleState( const double time ) = 0;

    //! Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (Time object time and double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    virtual Eigen::Matrix< double, 6, 1 > getBaseFrameDoubleState( const Time& time ) = 0;

    //! Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Pure virtual function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (Time object time and long double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    virtual Eigen::Matrix< long double, 6, 1 > getBaseFrameLongDoubleState( const Time& time ) = 0;

    //! Name of frame origin for which inertial state is computed by this class
    std::string baseFrameId_;
};

//! Class used for the determination of the inertial state of a Body's ephemeris origin
template< typename TimeType, typename StateScalarType >
class BaseStateInterfaceImplementation : public BaseStateInterface
{
public:
    //! Constructor using a default zero state function.
    /*!
     * Constructor using a default zero state function.
     * \param baseFrameId Name of frame origin for which inertial state is computed by this class
     */
    BaseStateInterfaceImplementation( const std::string baseFrameId ):
        BaseStateInterface( baseFrameId ), stateFunction_( &BaseStateInterfaceImplementation::getDefaultZeroState ),
        stateMultiplier_( 1 )
    { }

    //! Constructor
    /*!
     * Constructor
     * \param baseFrameId Name of frame origin for which inertial state is computed by this class
     * \param stateFunction Function returning frame's inertial state as a function of time.
     * \param subtractStateFunction Boolean denoting whether to subtract or add the state function (i.e. whether to multiply
     * result of stateFunction by -1).
     */
    BaseStateInterfaceImplementation( const std::string baseFrameId,
                                      const std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction,
                                      const bool subtractStateFunction = 0 ):
        BaseStateInterface( baseFrameId ), stateFunction_( stateFunction ), stateMultiplier_( ( subtractStateFunction == 0 ) ? 1.0 : -1.0 )
    { }

    //! Destructor
    ~BaseStateInterfaceImplementation( ) { }

protected:
    //! Function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (double time and double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    Eigen::Matrix< double, 6, 1 > getBaseFrameDoubleState( const double time )
    {
        return static_cast< double >( stateMultiplier_ ) * std::move( stateFunction_( time ) ).template cast< double >( );
    }

    //! Function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (double time and double long state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    Eigen::Matrix< long double, 6, 1 > getBaseFrameLongDoubleState( const double time )
    {
        return static_cast< long double >( stateMultiplier_ ) * stateFunction_( time ).template cast< long double >( );
    }

    //! Function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (Time object time and double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    Eigen::Matrix< double, 6, 1 > getBaseFrameDoubleState( const Time& time )
    {
        return static_cast< double >( stateMultiplier_ ) * stateFunction_( time ).template cast< double >( );
    }

    //! Function through which the state of baseFrameId_ in the inertial frame can be determined
    /*!
     *  Function through which the state of baseFrameId_ in the inertial frame can be determined
     *  (Time object time and long double state scalar).
     *  \param time Time at which state is to be computed
     *  \return Inertial state of frame origin at requested time
     */
    Eigen::Matrix< long double, 6, 1 > getBaseFrameLongDoubleState( const Time& time )
    {
        return static_cast< long double >( stateMultiplier_ ) * std::move( stateFunction_( time ) ).template cast< long double >( );
    }

private:
    //! Default state function returning a zero vector for all inputs.
    static Eigen::Matrix< StateScalarType, 6, 1 > getDefaultZeroState( const TimeType )
    {
        return Eigen::Matrix< StateScalarType, 6, 1 >::Zero( );
    }

    //! Function returning frame's inertial state as a function of time.
    std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction_;

    //! Value (1 or -1) by which to multiply the state returned by stateFunction_.
    int stateMultiplier_;
};

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_BASESTATEINTERFACE_H
