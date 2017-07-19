/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TIDALLOVENUMBERPARAMETER_H
#define TUDAT_TIDALLOVENUMBERPARAMETER_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"

using std::string;

namespace tudat
{

//! List of all orders in a single degree (starting at 2)
/*!
 *  List of all orders in a single degree (starting at 2), generated at compile time and used for tidal love number
 *  partial interface function calls over a full degree.
 */
static const std::vector< std::vector< int > > fullDegreeOrders = { { 0, 1, 2 }, { 0, 1, 2, 3 }, { 0, 1, 2, 3, 4 } };


namespace estimatable_parameters
{

template< typename ParameterScalar >
class TidalLoveNumber: public EstimatableParameter< ParameterScalar >
{
public:
    TidalLoveNumber(
            const boost::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariationModel,
            const std::string& associatedBody,
            const int degree,
            const std::vector< int > orders,
            const bool sumOrders,
            const EstimatebleParametersEnum loveNumberType,
            const bool useComplexComponents = 0 ):
        EstimatableParameter< ParameterScalar >( loveNumberType, associatedBody ),
        degree_( degree ),
        orders_( orders ),
        sumOrders_( sumOrders ),
        gravityFieldVariationModel_( gravityFieldVariationModel ),
        useComplexComponents_( useComplexComponents )
    { }

    virtual ~TidalLoveNumber( ){ }

    int getParameterSize( )
    {
        return parameterSize_;
    }

    int getDegree( )
    {
        return degree_;
    }

    bool useComplexComponents( )
    {
        return useComplexComponents_;
    }

    std::vector< std::string > getDeformingBodies( )
    {
        return gravityFieldVariationModel_->getDeformingBodies( );
    }

    std::vector< int > getOrders( )
    {
        return orders_;
    }

    bool getSumOrders( )
    {
        return sumOrders_;
    }

protected:

    int degree_;

    std::vector< int > orders_;

    bool sumOrders_;

    boost::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariationModel_;

    bool useComplexComponents_;

    int parameterSize_;

};

class FullDegreeTidalLoveNumber: public TidalLoveNumber< Eigen::VectorXd >
{
public:

    FullDegreeTidalLoveNumber(
            const boost::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariationModel,
            const std::string& associatedBody,
            const int degree,
            const bool useComplexComponents = 0 ):
        TidalLoveNumber< Eigen::VectorXd >(
            gravityFieldVariationModel, associatedBody, degree, fullDegreeOrders.at( degree - 2 ), 1, full_degree_tidal_love_number, useComplexComponents )
    {
        if( useComplexComponents_ )
        {
            parameterSize_ = 2;
        }
        else
        {
            parameterSize_ = 1;
        }
    }

    Eigen::VectorXd getParameterValue( );

    void setParameterValue( Eigen::VectorXd parameterValue );


private:
};

class SingleDegreeVariableTidalLoveNumber: public TidalLoveNumber< Eigen::VectorXd >
{
public:

    SingleDegreeVariableTidalLoveNumber(
            const boost::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariationModel,
            const std::string& associatedBody,
            const int degree,
            const std::vector< int > orders,
            const bool useComplexComponents = 0 ):
        TidalLoveNumber< Eigen::VectorXd >(
            gravityFieldVariationModel, associatedBody, degree, orders, 0, single_degree_variable_tidal_love_number, useComplexComponents )
    {
        if( useComplexComponents_ )
        {
            parameterSize_ = 2 * orders_.size( );
        }
        else
        {
            parameterSize_ = orders_.size( );
        }
    }

    Eigen::VectorXd getParameterValue( );

    void setParameterValue( Eigen::VectorXd parameterValue );


};

}

}

#endif // TUDAT_TIDALLOVENUMBERPARAMETER_H
