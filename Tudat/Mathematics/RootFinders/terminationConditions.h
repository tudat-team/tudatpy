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

#ifndef TUDAT_TERMINATION_CONDITIONS_H
#define TUDAT_TERMINATION_CONDITIONS_H

#include <cmath>
#include <stdexcept>
#include <string>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

namespace tudat
{
namespace root_finders
{
namespace termination_conditions
{

//! Check if maximum number of iterations specified has been exceeded.
/*!
 * Checks if maximum number of iterations specified has been exceeded. If the
 * throwRunTimeException, a run-time exception will be triggered, otherwise, the function will
 * return a false value.
 * \param numberOfIterations Number of iterations that have been completed.
 * \param maximumNumberOfIterations Maximum number of iterations specified.
 * \param throwRunTimeException Flag that indicates if run-time error should be triggered if
 *          maximum number of iterations is exceeded (default=true).
 * \return Flag indicating if maximum number of iterations have been exceeded.
 */
inline bool checkMaximumIterationsExceeded( const unsigned int numberOfIterations,
                                            const unsigned int maximumNumberOfIterations,
                                            const bool throwRunTimeException = true )
{
    // Check if maximum number of iterations have been exceeded.
    const bool isMaximumIterationsExceeded = numberOfIterations > maximumNumberOfIterations;

    // If exceeded, and flag is set to throw run-time exception, proceed.
    if ( isMaximumIterationsExceeded && throwRunTimeException )
    {
        std::string errorMessage
                = "Root-finder did not converge within maximum number of iterations!";

        throw std::runtime_error( errorMessage );
    }

    // Else, simply return whether maximum number of iterations have been exceeded.
    return isMaximumIterationsExceeded;
}

//! Check if absolute tolerance for root value has been achieved.
/*!
 * Checks if absolute tolerance has been achieved for successive root values, compared against the
 * specified absolute tolerance.
 * \param currentRootGuess Current value of root.
 * \param previousRootGuess Previous value of root.
 * \param absoluteTolerance Absolute tolerance.
 * \return Flag indicating if absolute tolerance has been achieved.
 */
template< typename ScalarType >
inline bool checkRootAbsoluteTolerance( const ScalarType currentRootGuess,
                                        const ScalarType previousRootGuess,
                                        const ScalarType absoluteTolerance )
{
    return std::fabs( currentRootGuess - previousRootGuess ) < absoluteTolerance;
}

//! Check if relative tolerance for root value has been achieved.
/*!
 * Checks if relative tolerance has been achieved for successive root values, compared against the
 * specified relative tolerance.
 * \param currentRootGuess Current value of root.
 * \param previousRootGuess Previous value of root.
 * \param relativeTolerance Relative tolerance.
 * \return Flag indicating if relative tolerance has been achieved.
 */
template< typename ScalarType >
inline bool checkRootRelativeTolerance( const ScalarType currentRootGuess,
                                        const ScalarType previousRootGuess,
                                        const ScalarType relativeTolerance )
{
    return std::fabs( ( currentRootGuess - previousRootGuess )
                      / currentRootGuess ) < relativeTolerance;
}

//! Termination condition base class.
/*!
 * This base class should be used for all termination condition classes, as it specifies the
 * checkTerminationCondition()-function used by the root-finders in Tudat.
 */
template< typename ScalarType = double >
class TerminationConditionBase
{
public:

    //! Virtual destructor.
    virtual ~TerminationConditionBase( ) { }

    //! Check termination condition.
    /*!
     * Checks termination condition with standard input arguments passed on by root-finder that
     * this function is used in. This function returns a flag indicating whether the termination
     * condition has been reached (true if the termination condition is achieved, indicating that
     * the root-finder should stop its search).
     * \param currentRootGuess Current root value.
     * \param previousRootGuess Previous root value.
     * \param currentRootFunctionValue Current root function value.
     * \param previousRootFunctionValue Previous root function value.
     * \param numberOfIterations Number of iterations that have been completed.
     * \return Flag indicating if termination condition has been reached.
     */
    virtual bool checkTerminationCondition( const ScalarType currentRootGuess,
                                            const ScalarType previousRootGuess,
                                            const ScalarType currentRootFunctionValue,
                                            const ScalarType previousRootFunctionValue,
                                            const unsigned int numberOfIterations ) = 0;
};

//! Maximum iterations termination condition.
/*!
 * This class implements the TerminationCondition base class to check if the maximum number of
 * iterations has been reached by the root-finder that this class is used in conjunction with. This
 * class effectively serves as a wrapper for the checkMaximumIterationsExceeded()-function.
 * \sa checkMaximumIterationsExceeded().
 */
class MaximumIterationsTerminationCondition : public TerminationConditionBase< double >
{
public:

    //! Default constructor, taking maximum number of iterations and run-time exception flag.
    /*!
     * Default constructor, taking a specified maximum number of iterations (default=1000) and a
     * flag indicating if a run-time exception should be thrown if this number is exceeded
     * (default=true).
     * \param aMaximumNumberOfIterations Maximum number of iterations (default=1000).
     * \param aThrowRunTimeExceptionFlag Flag that indicates if run-time error should be triggered
     * if maximum number of iterations is exceeded (default=true).
     */
    MaximumIterationsTerminationCondition( const unsigned int aMaximumNumberOfIterations = 1000,
                                           const bool aThrowRunTimeExceptionFlag = true )
        : maximumNumberOfIterations( aMaximumNumberOfIterations ),
          throwRunTimeException( aThrowRunTimeExceptionFlag )
    { }

    //! Check termination condition (wrapper for checkMaximumIterationsExceeded()-function).
    bool checkTerminationCondition( const double currentRootGuess,
                                    const double previousRootGuess,
                                    const double currentRootFunctionValue,
                                    const double previousRootFunctionValue,
                                    const unsigned int numberOfIterations )
    {
        return checkMaximumIterationsExceeded( numberOfIterations, maximumNumberOfIterations,
                                               throwRunTimeException );
    }

private:

    //! Maximum number of iterations.
    const unsigned int maximumNumberOfIterations;

    //! Flag indicating if run-time exception should be thrown.
    const bool throwRunTimeException;
};

//! Absolute tolerance for root termination condition.
/*!
 * This class implements the TerminationCondition base class to check if the absolute tolerance for
 * the root value computed by the root-finder that this class is used in conjunction with has been
 * reached. Additionally, this termination condition also checks if the maximum number of
 * iterations has been exceeded. This class effectively serves as a wrapper for the
 * checkRootAbsoluteTolerance()- and checkMaximumIterationsExceeded()-functions.
 * \sa checkRootAbsoluteTolerance(), checkMaximumIterationsExceeded().
 */
template< typename ScalarType = double >
class RootAbsoluteToleranceTerminationCondition : public TerminationConditionBase< ScalarType >
{
public:

    //! Default constructor, taking root absolute tolerance, maximum number of iterations and
    //! run-time exception flag.
    /*!
     * Default constructor, taking a root absolute tolerance (default=1.0e-15), a specified maximum
     * number of iterations (default=1000) and a flag indicating if a run-time exception should be
     * thrown if this number is exceeded (default=true).
     * \param anAbsoluteTolerance An absolute tolerance (default=1.0e-15).
     * \param aMaximumNumberOfIterations Maximum number of iterations (default=1000).
     * \param aThrowRunTimeExceptionFlag Flag that indicates if run-time error should be triggered
     * if maximum number of iterations is exceeded (default=true).
     */
    RootAbsoluteToleranceTerminationCondition( const ScalarType anAbsoluteTolerance = 1.0e-15,
                                               const unsigned int aMaximumNumberOfIterations = 1000,
                                               const bool aThrowRunTimeExceptionFlag = true )
        : absoluteTolerance( anAbsoluteTolerance ),
          maximumNumberOfIterations( aMaximumNumberOfIterations ),
          throwRunTimeException( aThrowRunTimeExceptionFlag )
    { }

    //! Check termination condition (combined abs. tolerance and max. itetations)
    /*!
     * Check termination condition (wrapper for checkRootAbsoluteTolerance()- and
     * checkMaximumIterationsExceeded()-functions).
     * \param currentRootGuess Current root value.
     * \param previousRootGuess Previous root value.
     * \param currentRootFunctionValue Current root function value (not used, included for
     * base class compatibility).
     * \param previousRootFunctionValue Previous root function value (not used, included for
     * base class compatibility).
     * \param numberOfIterations Number of iterations that have been completed.
     * \return Flag indicating if termination condition has been reached.
     */
    bool checkTerminationCondition( const ScalarType currentRootGuess,
                                    const ScalarType previousRootGuess,
                                    const ScalarType currentRootFunctionValue,
                                    const ScalarType previousRootFunctionValue,
                                    const unsigned int numberOfIterations )
    {
        return checkRootAbsoluteTolerance< ScalarType >(
                    currentRootGuess, previousRootGuess, absoluteTolerance )
                || checkMaximumIterationsExceeded(
                    numberOfIterations, maximumNumberOfIterations, throwRunTimeException );
    }

private:

    //! Absolute tolerance.
    const ScalarType absoluteTolerance;

    //! Maximum number of iterations.
    const unsigned int maximumNumberOfIterations;

    //! Flag indicating if run-time exception should be thrown.
    const bool throwRunTimeException;
};

//! Relative tolerance for root termination condition.
/*!
 * This class implements the TerminationCondition base class to check if the relative tolerance for
 * the root value computed by the root-finder that this class is used in conjunction with has been
 * reached. Additionally, this termination condition also checks if the maximum number of
 * iterations has been exceeded. This class effectively serves as a wrapper for the
 * checkRootRelativeTolerance()- and checkMaximumIterationsExceeded()-functions.
 * \sa checkRootRelativeTolerance(), checkMaximumIterationsExceeded().
 */
template< typename ScalarType = double >
class RootRelativeToleranceTerminationCondition : public TerminationConditionBase< ScalarType >
{
public:

    //! Default constructor, taking root relative tolerance, maximum number of iterations and
    //! run-time exception flag.
    /*!
     * Default constructor, taking a root relative tolerance (default=1.0e-12), a specified maximum
     * number of iterations (default=1000) and a flag indicating if a run-time exception should be
     * thrown if this number is exceeded (default=true).
     * \param aRelativeTolerance A relative tolerance (default=1.0e-12).
     * \param aMaximumNumberOfIterations Maximum number of iterations (default=1000).
     * \param aThrowRunTimeExceptionFlag Flag that indicates if run-time error should be triggered
     *        if maximum number of iterations is exceeded (default=true).
     */
    RootRelativeToleranceTerminationCondition( const ScalarType aRelativeTolerance = 1.0e-12,
                                               const unsigned int aMaximumNumberOfIterations = 1000,
                                               const ScalarType aThrowRunTimeExceptionFlag = true )
        : relativeTolerance( aRelativeTolerance ),
          maximumNumberOfIterations( aMaximumNumberOfIterations ),
          throwRunTimeException( aThrowRunTimeExceptionFlag )
    { }

    //! Check termination condition (combined rel. tolerance and max. itetations)
    /*!
     * Check termination condition (wrapper for checkRootRelativeTolerance()- and
     * checkMaximumIterationsExceeded()-functions).
     * \param currentRootGuess Current root value.
     * \param previousRootGuess Previous root value.
     * \param currentRootFunctionValue Current root function value (not used, included for
     * base class compatibility).
     * \param previousRootFunctionValue Previous root function value (not used, included for
     * base class compatibility).
     * \param numberOfIterations Number of iterations that have been completed.
     * \return Flag indicating if termination condition has been reached.
     */
    bool checkTerminationCondition( const ScalarType currentRootGuess,
                                    const ScalarType previousRootGuess,
                                    const ScalarType currentRootFunctionValue,
                                    const ScalarType previousRootFunctionValue,
                                    const unsigned int numberOfIterations )
    {
        return checkRootRelativeTolerance( currentRootGuess, previousRootGuess, relativeTolerance )
                || checkMaximumIterationsExceeded(
                    numberOfIterations, maximumNumberOfIterations, throwRunTimeException );
    }

private:

    //! Relative tolerance.
    const ScalarType relativeTolerance;

    //! Maximum number of iterations.
    const unsigned int maximumNumberOfIterations;

    //! Flag indicating if run-time exception should be thrown.
    const bool throwRunTimeException;
};

//! Absolute or relative tolerance for root termination condition.
/*!
 * This class implements the TerminationCondition base class to check if the absolute or relative
 * tolerance for the root value computed by the root-finder that this class is used in conjunction
 * with has been reached. Additionally, this termination condition also checks if the maximum
 * number of iterations has been exceeded. This class effectively serves as a wrapper for the
 * checkRootAbsoluteTolerance()-, checkRootRelativeTolerance()- and
 * checkMaximumIterationsExceeded()-functions.
 * \sa checkRootAbsoluteTolerance(), checkRootRelativeTolerance(),
 *      checkMaximumIterationsExceeded().
 */
template< typename ScalarType = double >
class RootAbsoluteOrRelativeToleranceTerminationCondition :
        public TerminationConditionBase< ScalarType >
{
public:

    //! Default constructor, taking root absolute and relative tolerances, maximum number of
    //! iterations and run-time exception flag.
    /*!
     * Default constructor, taking a root absolute tolerance (default=1.0e-15), root relative
     * tolerance (default=1.0e-12), a specified maximum number of iterations (default=1000) and a
     * flag indicating if a run-time exception should be thrown if this number is exceeded
     * (default=true).
     * \param anAbsoluteTolerance A absolute tolerance (default=1.0e-12).
     * \param aRelativeTolerance A relative tolerance (default=1.0e-12).
     * \param aMaximumNumberOfIterations Maximum number of iterations (default=1000).
     * \param aThrowRunTimeExceptionFlag Flag that indicates if run-time error should be triggered
     *        if maximum number of iterations is exceeded (default=true).
     */
    RootAbsoluteOrRelativeToleranceTerminationCondition(
            const ScalarType anAbsoluteTolerance = 1.0e-15,
            const ScalarType aRelativeTolerance = 1.0e-12,
            const unsigned int aMaximumNumberOfIterations = 1000,
            const bool aThrowRunTimeExceptionFlag = true )
        : absoluteTolerance( anAbsoluteTolerance ),
          relativeTolerance( aRelativeTolerance ),
          maximumNumberOfIterations( aMaximumNumberOfIterations ),
          throwRunTimeException( aThrowRunTimeExceptionFlag )
    { }

    //! Check termination condition (combined ans. and rel. tolerance and max. itetations)
    /*!
     * Check termination condition (wrapper for checkRootRelativeTolerance()-
     * checkRootAbsoluteTolerance()- and checkMaximumIterationsExceeded()-functions).
     * \param currentRootGuess Current root value.
     * \param previousRootGuess Previous root value.
     * \param currentRootFunctionValue Current root function value (not used, included for
     * base class compatibility).
     * \param previousRootFunctionValue Previous root function value (not used, included for
     * base class compatibility).
     * \param numberOfIterations Number of iterations that have been completed.
     * \return Flag indicating if termination condition has been reached.
     */
    bool checkTerminationCondition( const ScalarType currentRootGuess,
                                    const ScalarType previousRootGuess,
                                    const ScalarType currentRootFunctionValue,
                                    const ScalarType previousRootFunctionValue,
                                    const unsigned int numberOfIterations )
    {
        return checkRootAbsoluteTolerance( currentRootGuess, previousRootGuess, absoluteTolerance )
                || checkRootRelativeTolerance( currentRootGuess, previousRootGuess,
                                               relativeTolerance )
                || checkMaximumIterationsExceeded(
                    numberOfIterations, maximumNumberOfIterations, throwRunTimeException );
    }

private:

    //! Absolute tolerance.
    const ScalarType absoluteTolerance;

    //! Relative tolerance.
    const ScalarType relativeTolerance;

    //! Maximum number of iterations.
    const unsigned int maximumNumberOfIterations;

    //! Flag indicating if run-time exception should be thrown.
    const bool throwRunTimeException;
};

//! Check termination condition (required maximum absolute value of root function)
/*!
 * Check termination condition (required maximum absolute value of root function)
 * \param currentRootGuess Current root value.
 * \param previousRootGuess Previous root value.
 * \param currentRootFunctionValue Current root function value (not used, included for
 * base class compatibility).
 * \param previousRootFunctionValue Previous root function value (not used, included for
 * base class compatibility).
 * \param numberOfIterations Number of iterations that have been completed.
 * \param rootToleranceValue Maximum allowed value fo root function
 * \return Flag indicating if termination condition has been reached.
 */
template< typename ScalarType = double >
bool checkRootFunctionValueCondition( const ScalarType currentRootGuess,
                                      const ScalarType previousRootGuess,
                                      const ScalarType currentRootFunctionValue,
                                      const ScalarType previousRootFunctionValue,
                                      const unsigned int numberOfIterations,
                                      const ScalarType rootToleranceValue )
{
    return (std::fabs( currentRootFunctionValue ) < rootToleranceValue );
}

//! Typedef for shared-pointer to MaximumIterationsTerminationCondition object.
typedef boost::shared_ptr< MaximumIterationsTerminationCondition >
MaximumIterationsTerminationConditionPointer;

} // namespace termination_conditions
} // namespace root_finders
} // namespace tudat

#endif // TUDAT_TERMINATION_CONDITIONS_H
