/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#ifndef TUDAT_RUNGE_KUTTA_COEFFICIENTS_H
#define TUDAT_RUNGE_KUTTA_COEFFICIENTS_H

#include <memory>

#include <Eigen/Core>

#include <iostream>

namespace tudat
{
namespace numerical_integrators
{

// Struct that defines the coefficients of a Runge-Kutta integrator
/*
 * Struct that defines coefficients of a Runge-Kutta integrators
 */
struct RungeKuttaCoefficients
{
    // Enum of order estimates that can be integrated.
    enum OrderEstimateToIntegrate { lower, higher };

    // Main table of the Butcher tableau.
    Eigen::MatrixXd aCoefficients;

    // Bottom rows of the Butcher tableau.
    Eigen::MatrixXd bCoefficients;

    // First column of the Butcher tableau.
    Eigen::VectorXd cCoefficients;

    // Order of the higher order estimate.
    unsigned int higherOrder;

    // Order of the lower order estimate.
    unsigned int lowerOrder;

    // Order estimate to integrate.
    OrderEstimateToIntegrate orderEstimateToIntegrate;

    // Name of the coefficients.
    std::string name;

    // Default constructor.
    /*
     * Default constructor that initializes coefficients to 0.
     */
    RungeKuttaCoefficients( ) :
        aCoefficients( ),
        bCoefficients( ),
        cCoefficients( ),
        higherOrder( 0 ),
        lowerOrder( 0 ),
        orderEstimateToIntegrate( lower ),
        name( "Undefined" )        
    { }

    // Constructor.
    /*
     * Constructor that sets the coefficients.
     * \param aCoefficients_ Main table of the Butcher tableau.
     * \param bCoefficients_ Bottom rows of the Butcher tableau.
     * \param cCoefficients_ First column of the Butcher tableau.
     * \param higherOrder_ Order of the integrator.
     * \param lowerOrder_ Order of the embedded low-order integrator.
     * \param order Enum denoting whether to use the lower or higher order scheme for numerical
     * integration.
     */
    RungeKuttaCoefficients( const Eigen::MatrixXd& aCoefficients_,
                            const Eigen::MatrixXd& bCoefficients_,
                            const Eigen::MatrixXd& cCoefficients_,
                            const unsigned int higherOrder_,
                            const unsigned int lowerOrder_,
                            OrderEstimateToIntegrate order,
                            std::string name_ = "Undefined" ) :
        aCoefficients( aCoefficients_ ),
        bCoefficients( bCoefficients_ ),
        cCoefficients( cCoefficients_ ),
        higherOrder( higherOrder_ ),
        lowerOrder( lowerOrder_ ),
        orderEstimateToIntegrate( order ),
        name( name_ )
    { }

    // Enum of predefined coefficient sets.
    //! @get_docstring(CoefficientSets.__docstring__)
    enum CoefficientSets
    {
        undefinedCoefficientSet = -1,
        rungeKuttaFehlberg45,
        rungeKuttaFehlberg56,
        rungeKuttaFehlberg78,
        rungeKutta87DormandPrince
    };

    // Enum of predefined fixed-step coefficient sets.
    enum FixedStepCoefficientSets
    {
        undefinedFixedStepCoefficientSet = -1,
        forwardEuler,
        rungeKutta4
    };

    // Get coefficients for a specified coefficient set.
    /*
     * Returns coefficients for a specified coefficient set.
     * \param coefficientSet The set to get the coefficients for.
     * \return The requested coefficient set.
     */
    static const RungeKuttaCoefficients& get( CoefficientSets coefficientSet );

    // Get coefficients for a specified coefficient set for fixed-step integration.
    /*
     * Returns coefficients for a specified coefficient set for fixed-step integration.
     * \param coefficientSet The set to get the coefficients for.
     * \return The requested coefficient set.
     */
    static const RungeKuttaCoefficients& get( FixedStepCoefficientSets coefficientSet );


    void printButcherTableau( )
    {
        std::cout << "Butcher tableau of the " << name << " coefficients: " << std::endl;

        // Create a zero matrix of the same size as aCoefficients plus 1.
        Eigen::MatrixXd ButcherTable = Eigen::MatrixXd::Zero( aCoefficients.rows( ) + bCoefficients.rows( ), bCoefficients.cols( ) + 1 );
        
        // Set the table first column to the values of cCoefficients.
        ButcherTable.block( 0, 0, cCoefficients.rows( ), 1 ) = cCoefficients;

        // Set the table last row(s) to the values of bCoefficients.
        ButcherTable.block( aCoefficients.rows( ), 1, bCoefficients.rows( ), bCoefficients.cols( ) ) = bCoefficients;
        // ButcherTable.block( bCoefficients.cols( ), 1, cCoefficients.cols( ), cCoefficients.rows( ) ) = bCoefficients;

        // Set the rest of the table to the values of aCoefficients.
        ButcherTable.block( 0, 1, aCoefficients.rows( ), aCoefficients.cols( ) ) = aCoefficients;
        
        // Feed the full Butcher tableau into a stringstream (to make use of the precision/formatting from IOFormat implemented in Eigen).
        std::stringstream tableStream;
        tableStream << ButcherTable;
        int line_i = 0;
        int first_column_end = -1;
        std::string line;
        // Go trough each of the table lines.
        while(std::getline(tableStream,line,'\n'))
        {
            // If the index at which the first column ends is not known yet, find it.
            if (first_column_end == -1)
            {
                bool has_encountered_non_space = false;
                // Go trough each of the line characters.
                for (unsigned int i = 0; i < line.length( ); i++)
                {
                    // If the character is not a space, remember it.
                    if (line[i] != ' ')
                    {
                        has_encountered_non_space = true;
                    }
                    // If the character is a space and we have encountered a non-space character before...
                    if (line[i] == ' ' && has_encountered_non_space)
                    {
                        // Remember the index at which the first column ends.
                        first_column_end = i;
                        break;
                    }
                }
            }

            // Add a vertical bar after the end of the first column.
            line.insert(first_column_end + 1, "| ");

            // If we are in the last row(s)...
            if (line_i >= ButcherTable.rows( ) - bCoefficients.rows( ))
            {
                // Replace every character in the line before the first column end by a space (do not print 0 in bottom left tableau section).
                for (int i = 0; i < first_column_end; i++)
                {
                    line[i] = ' ';
                }                
            }

            // Print the line.
            std::cout << line << std::endl;

            // Print a dash line to show the separation between a and b coefficients.
            if(line_i == ButcherTable.rows( ) - bCoefficients.rows( ) - 1)
            {
                std::string dash_line(line.length( ) - 1, '-');
                dash_line.insert(first_column_end + 1, "|");
                std::cout << dash_line << std::endl;
            }
            line_i++;
        }
    }

};

// Typedef for shared-pointer to RungeKuttaCoefficients object.
typedef std::shared_ptr< RungeKuttaCoefficients > RungeKuttaCoefficientsPointer;

} // namespace numerical_integrators
} // namespace tudat

#endif // TUDAT_RUNGE_KUTTA_COEFFICIENTS_H
