/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    NOTE: The code in this file, and the associated cpp file, is tested in the unitTestDependentVariableOutput.cpp file.
 */

#ifndef TUDAT_TRIMORIENTATION_H
#define TUDAT_TRIMORIENTATION_H

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/math/root_finders/rootFinder.h"

namespace tudat
{

namespace aerodynamics
{



//! Class to determine the trimmed angle-of-attack for a given set of aerodynamic coefficients.
/*!
 *  Class to determine the trimmed angle-of-attack for a given set of aerodynamic coefficients. The coefficient interface
 *  provided as input must be dependent on the angle of attack for this class to function.
 */
class TrimOrientationCalculator
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param coefficientInterface Object containing used to retrieve aerodynamic coefficients as function of independent
     * variables.
     * \param rootFinder Object to iteratively find the root of the equations C_m(alpha)=0, i.e. to determine the
     * angle of attack for which the pitch moment is zero.
     */
    TrimOrientationCalculator(
            const std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface,
            const std::shared_ptr< root_finders::RootFinder< double > > rootFinder = nullptr );

    //! Function to find the trimmed angle of attack for a given set of independent  variables
    /*!
     * Function to find the trimmed angle of attack for a given set of independent  variables. This function iteratively
     * changes the angle of attack in the untrimmedIndependentVariables vector (keeping any other variables the same) and
     * returns the value of the trimmed angle of attack. Note that this function will typically have some small numerical
     * error in the result, as a result of the error tolerances in the root finder.
     * \param untrimmedIndependentVariables Untrimmed list of independent variables (in order required as input for
     * coefficientInterface_)
     * \param untrimmedControlSurfaceIndependentVariables Untrimmed list of independent variables for control surfaces
     * with map key denoting the control surface name (in order required as input for coefficient interfaces)
     * \return Trimmed angle of attack.
     */
    double findTrimAngleOfAttack(
            const std::vector< double > untrimmedIndependentVariables,
            const std::map< std::string, std::vector< double > > untrimmedControlSurfaceIndependentVariables =
            std::map< std::string, std::vector< double > >( ) );

    //! Function to find the trimmed angle of attack for a given set of independent  variables
    /*!
     * Function to find the trimmed angle of attack for a given set of independent  variables. This function iteratively
     * changes the angle of attack in the untrimmedIndependentVariables vector (keeping any other variables the same) and
     * returns the value of the trimmed angle of attack. Note that this function will typically have some small numerical
     * error in the result, as a result of the error tolerances in the root finder.
     * \param untrimmedIndependentVariablesFunction Function returning untrimmed list of independent variables
     * (in order required as input for coefficientInterface_
     * \param untrimmedControlSurfaceIndependentVariablesFunction Function returning untrimmed list of independent variables
     * for control surfaces with map key denoting the control surface name (in order required as input for coefficient
     * interfaces)
     * \return Trimmed angle of attack.
     */
    double findTrimAngleOfAttackFromFunction(
            const std::function< std::vector< double >( ) > untrimmedIndependentVariablesFunction,
            const std::function< std::map< std::string, std::vector< double > >( ) >
            untrimmedControlSurfaceIndependentVariablesFunction )
    {
        return findTrimAngleOfAttack( untrimmedIndependentVariablesFunction( ),
                                      untrimmedControlSurfaceIndependentVariablesFunction( ) );
    }

private:

    //! Function to get the moment coefficient for a given angle of attack
    /*!
     * Function to get the moment coefficient for a given perturbed angle of attack, keeping all other independent  variables
     * constant. This function is used as input to the root finder to determine the trim point.
     * \param perturbedAngleOfAttack Angle of attack to use
     * \param unperturbedConditions Untrimmed list of independent variables (in order required as input for
     * coefficientInterface_
     * \param unperturbedControlSurfaceIndependentVariables Untrimmed list of independent variables for control surfaces
     * with map key denoting the control surface name (in order required as input for coefficient interfaces)
     * \return Moment coefficient at given independent variables.
     */
    double getPerturbedMomentCoefficient(
            const double perturbedAngleOfAttack,
            const std::vector< double >& unperturbedConditions,
            const std::map< std::string, std::vector< double > > unperturbedControlSurfaceIndependentVariables );

    //! Object containing used to retrieve aerodynamic coefficients as function of independent variables.
    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface_;

    //! Object to iteratively find the root of the equations C_m(alpha)=0, i.e. to determine the
    //!  angle of attack for which the pitch moment is zero.
    std::shared_ptr< root_finders::RootFinder< double > > rootFinder_;

    //! Index in independent variable list of coefficientInterface_ corresponding to the angle of attack.
    int variableIndex_;

    //! Index in list of each of the control surface interfaces corresponding to the angle of attack.
    std::map< std::string, int > controlSurfaceVariableIndex_;
};

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_TRIMORIENTATION_H
