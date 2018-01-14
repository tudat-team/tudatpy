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

#ifndef TUDAT_CONVERGENCE_EXEPTION_H
#define TUDAT_CONVERGENCE_EXEPTION_H

#include <stdexcept>


namespace tudat
{
namespace basic_mathematics
{

//! Exception indicating that an algorithm could not complete convergence to a solution.
/*!
 * You can add extra information using;
 * void f()
 * {
 *   throw ConvergenceExeption() << my_info(42);
 * }
 *
 * \sa Boost.Exception.
 */
struct ConvergenceException : public virtual boost::exception, public virtual std::exception
{
public:

    //! Constructor that sets exception message to default string.
    ConvergenceException( )
        : message( "Failed to converge to a solution.")
    { }

    //! Constructor that sets exception message based on input character array.
    ConvergenceException( const char* errMessage )
        : message( errMessage )
    { }

    //! Constructor that sets exception message based on input string.
    ConvergenceException( const std::string& errMessage )
        : message( errMessage.c_str( ) )
    { }

    //! Return what the exception message stored is.
    const char* what( ) const throw( ) { return message; }

protected:

private:

    //! Stored exception message.
    const char* message;
};

} // namespace basic_mathematics
} // namespace tudat

#endif  // TUDAT_CONVERGENCE_EXEPTION_H
