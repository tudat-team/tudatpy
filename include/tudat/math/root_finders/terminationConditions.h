/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <iostream>

#include <functional>
#include <memory>
#include <vector>

namespace tudat
{
namespace root_finders
{


enum MaximumIterationHandling
{
    accept_result,
    accept_result_with_warning,
    throw_exception
};


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
                                            const MaximumIterationHandling maximumIterationHandling )
{
    // Check if maximum number of iterations have been exceeded.
    const bool isMaximumIterationsExceeded = numberOfIterations > maximumNumberOfIterations;

    // If exceeded, and flag is set to throw run-time exception, proceed.
    if ( isMaximumIterationsExceeded )
    {
        static std::string errorMessage
                = "Root-finder did not converge within maximum number of iterations!";
        switch( maximumIterationHandling )
        {
        case accept_result:
            break;
        case accept_result_with_warning:
            std::cerr<<errorMessage<<std::endl;
            break;
        case throw_exception:
            throw std::runtime_error( errorMessage );
            break;
        }

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



//! Check termination condition (required maximum absolute value of root function)
/*!
 * Check termination condition (required maximum absolute value of root function)
 * \param currentRootFunctionValue Current root function value
 * \param rootToleranceValue Maximum allowed value fo root function
 * \return Flag indicating if termination condition has been reached.
 */
template< typename ScalarType = double >
bool checkRootFunctionValueCondition( const ScalarType currentRootFunctionValue,
                                      const ScalarType rootToleranceValue )
{
    return (std::fabs( currentRootFunctionValue ) < rootToleranceValue );
}

//! Termination condition base class.
/*!
 * This base class should be used for all termination condition classes, as it specifies the
 * checkTerminationCondition()-function used by the root-finders in Tudat.
 */
template< typename ScalarType = double >
class TerminationCondition
{
public:

    //! Virtual destructor.
    virtual ~TerminationCondition( ) { }

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
template< typename ScalarType = double >
class MaximumIterationsTerminationCondition : public TerminationCondition< ScalarType >
{
public:

    //! Default constructor, taking maximum number of iterations and run-time exception flag.
    /*!
     * Default constructor, taking a specified maximum number of iterations (default=1000) and a
     * flag indicating if a run-time exception should be thrown if this number is exceeded
     * (default=true).
     * \param maximumNumberOfIterations Maximum number of iterations (default=1000).
     */
    MaximumIterationsTerminationCondition( const unsigned int maximumNumberOfIterations = 1000,
                                           const MaximumIterationHandling maximumIterationHandling = throw_exception )
        : maximumNumberOfIterations_( maximumNumberOfIterations ),
          maximumIterationHandling_( maximumIterationHandling )
    { }

    //! Check termination condition (wrapper for checkMaximumIterationsExceeded()-function).
    bool checkTerminationCondition( const ScalarType currentRootGuess,
                                    const ScalarType previousRootGuess,
                                    const ScalarType currentRootFunctionValue,
                                    const ScalarType previousRootFunctionValue,
                                    const unsigned int numberOfIterations )
    {
        return checkMaximumIterationsExceeded( numberOfIterations, maximumNumberOfIterations_,
                                               maximumIterationHandling_ );
    }

private:

    //! Maximum number of iterations.
    const unsigned int maximumNumberOfIterations_;

    //! Flag indicating if run-time exception should be thrown.
    const MaximumIterationHandling maximumIterationHandling_;
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
class RootAbsoluteToleranceTerminationCondition : public TerminationCondition< ScalarType >
{
public:

    //! Default constructor, taking root absolute tolerance, maximum number of iterations and
    //! run-time exception flag.
    /*!
     * Default constructor, taking a root absolute tolerance (default=1.0e-15), a specified maximum
     * number of iterations (default=1000) and a flag indicating if a run-time exception should be
     * thrown if this number is exceeded (default=true).
     * \param anAbsoluteTolerance An absolute tolerance (default=1.0e-15).
     */
    RootAbsoluteToleranceTerminationCondition( const ScalarType absoluteTolerance = 1.0e-15 )
        : absoluteTolerance_( absoluteTolerance )
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
                    currentRootGuess, previousRootGuess, absoluteTolerance_ );
    }

private:

    //! Absolute tolerance.
    const ScalarType absoluteTolerance_;
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
class RootRelativeToleranceTerminationCondition : public TerminationCondition< ScalarType >
{
public:

    //! Default constructor, taking root relative tolerance, maximum number of iterations and
    //! run-time exception flag.
    /*!
     * Default constructor, taking a root relative tolerance (default=1.0e-12), a specified maximum
     * number of iterations (default=1000) and a flag indicating if a run-time exception should be
     * thrown if this number is exceeded (default=true).
     * \param relativeTolerance A relative tolerance (default=1.0e-12).
     */
    RootRelativeToleranceTerminationCondition( const ScalarType relativeTolerance = 1.0e-12 )
        : relativeTolerance_( relativeTolerance )
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
        return checkRootRelativeTolerance( currentRootGuess, previousRootGuess, relativeTolerance_ );
    }

private:

    //! Relative tolerance.
    const ScalarType relativeTolerance_;
};

template< typename ScalarType = double >
class RootFunctionTerminationCondition : public TerminationCondition< ScalarType >
{
public:

    RootFunctionTerminationCondition( const ScalarType rootFunctionTolerance = 1.0e-12 )
        : rootFunctionTolerance_( rootFunctionTolerance )
    { }

    bool checkTerminationCondition( const ScalarType currentRootGuess,
                                    const ScalarType previousRootGuess,
                                    const ScalarType currentRootFunctionValue,
                                    const ScalarType previousRootFunctionValue,
                                    const unsigned int numberOfIterations )
    {
        return checkRootFunctionValueCondition( currentRootFunctionValue, rootFunctionTolerance_ );
    }

private:

    const ScalarType rootFunctionTolerance_;

};

template< typename ScalarType = double >
class CombinedTerminationCondition : public TerminationCondition< ScalarType >
{
public:

    CombinedTerminationCondition(
            const std::vector< std::shared_ptr< TerminationCondition< ScalarType > > > terminationList )
        : terminationList_( terminationList )
    { }

    bool checkTerminationCondition( const ScalarType currentRootGuess,
                                    const ScalarType previousRootGuess,
                                    const ScalarType currentRootFunctionValue,
                                    const ScalarType previousRootFunctionValue,
                                    const unsigned int numberOfIterations )
    {
        bool terminate = false;
        for( unsigned int i = 0; i < terminationList_.size( ); i++ )
        {
            if( terminationList_.at( i )->checkTerminationCondition(
                        currentRootGuess, previousRootGuess,
                        currentRootFunctionValue, previousRootFunctionValue,
                        numberOfIterations ) )
            {
                terminate = true;
            }
        }
        return terminate;
    }

private:

    std::vector< std::shared_ptr< TerminationCondition< ScalarType > > > terminationList_;
};


template< typename DataType = double >
std::shared_ptr< TerminationCondition< DataType > > createTerminationCondition(
        const DataType relativeIndependentVariableTolerance = TUDAT_NAN,
        const DataType absoluteIndependentVariableTolerance = TUDAT_NAN,
        const DataType rootFunctionTolerance = TUDAT_NAN,
        const unsigned int maximumNumberOfIterations = 1000,
        const MaximumIterationHandling maximumIterationHandling = throw_exception )
{
    std::vector< std::shared_ptr< TerminationCondition< DataType > > > terminationConditionList;


    terminationConditionList.push_back(
                std::make_shared< MaximumIterationsTerminationCondition< DataType > >(
                    maximumNumberOfIterations, maximumIterationHandling ) );

    if( ( relativeIndependentVariableTolerance == relativeIndependentVariableTolerance ) )
    {
        terminationConditionList.push_back(
                    std::make_shared< RootRelativeToleranceTerminationCondition< DataType > >(
                        relativeIndependentVariableTolerance ) );
    }

    if( ( absoluteIndependentVariableTolerance == absoluteIndependentVariableTolerance ) )
    {
        terminationConditionList.push_back(
                    std::make_shared< RootAbsoluteToleranceTerminationCondition< DataType > >(
                        absoluteIndependentVariableTolerance ) );
    }

    if( ( rootFunctionTolerance == rootFunctionTolerance ) )
    {
        terminationConditionList.push_back(
                    std::make_shared< RootFunctionTerminationCondition< DataType > >(
                        rootFunctionTolerance ) );
    }

    if( terminationConditionList.size( ) == 1 )
    {
        return terminationConditionList.at( 0 );
    }
    else
    {
        return std::make_shared< CombinedTerminationCondition< DataType > >(
                    terminationConditionList );
    }
}

template< typename DataType = double >
std::function< bool( DataType, DataType, DataType, DataType, unsigned int ) > createTerminationConditionFunction(
        const DataType relativeIndependentVariableTolerance = TUDAT_NAN,
        const DataType absoluteIndependentVariableTolerance = TUDAT_NAN,
        const DataType rootFunctionTolerance = TUDAT_NAN,
        const unsigned int maximumNumberOfIterations = 1000,
        const MaximumIterationHandling maximumIterationHandling = throw_exception )
{
    return std::bind( &TerminationCondition< DataType >::checkTerminationCondition,
                      createTerminationCondition< DataType >(
                          relativeIndependentVariableTolerance, absoluteIndependentVariableTolerance, rootFunctionTolerance,
                          maximumNumberOfIterations, maximumIterationHandling ),
                      std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
                      std::placeholders::_4, std::placeholders::_5 );
}


} // namespace root_finders
} // namespace tudat

#endif // TUDAT_TERMINATION_CONDITIONS_H
