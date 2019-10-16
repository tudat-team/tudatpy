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

#ifndef TUDAT_COMPOSITE_FUNCTION_HODOGRAPHIC_SHAPING_H
#define TUDAT_COMPOSITE_FUNCTION_HODOGRAPHIC_SHAPING_H

#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include <math.h>
#include <vector>
#include <Eigen/Dense>

namespace tudat
{
namespace shape_based_methods
{

class CompositeFunctionHodographicShaping
{
public:

    //! Constructor.
    CompositeFunctionHodographicShaping( std::vector< std::shared_ptr< BaseFunctionHodographicShaping > >& compositeFunctionComponents,
                       Eigen::VectorXd& compositeFunctionCoefficients ){

        compositeFunctionComponents_ = compositeFunctionComponents;
        compositeFunctionCoefficients_ = compositeFunctionCoefficients;

    }


    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CompositeFunctionHodographicShaping( ) { }

    //! Return coefficients.
    /*!
     * Returns coefficients of composite function.
     */
    Eigen::VectorXd getCompositeFunctionCoefficients()
    {
        return compositeFunctionCoefficients_;
    }

    //! Reset coefficients.
    /*!
     * Resets the coefficients of the composite function.
     */
    void resetCompositeFunctionCoefficients( Eigen::VectorXd compositeFunctionCoefficients );

    //! Returns number of composite function components.
    /*!
     * Returns number of composite function components.
     */
    int getNumberOfCompositeFunctionComponents(){

        return compositeFunctionComponents_.size();
    }

    double evaluateCompositeFunctionCurrentTime( const double independentVariable );

    double evaluateCompositeFunctionDerivativeCurrentTime( const double independentVariable );

    double evaluateCompositeFunctionIntegralCurrentTime( const double independentVariable );

    double getComponentFunctionCurrentValue( int componentIndex, double currentTime );

    double getComponentFunctionDerivativeCurrentTime( int componentIndex, double currentTime );

    double getComponentFunctionIntegralCurrentTime( int componentIndex, double currentTime );


protected:

private:

    std::vector<  std::shared_ptr< BaseFunctionHodographicShaping > > compositeFunctionComponents_;

    double frequency_;
    double scaleFactor_;

    int numberOfIntegrationSteps_;
    int numberOfComponents_;

    //! Vector of coefficients.
    /*!
     * VectorXd containing coefficients corresponding to components of the composite function.
     */
    Eigen::VectorXd compositeFunctionCoefficients_;

};


} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_COMPOSITE_FUNCTION_HODOGRAPHIC_SHAPING_H
