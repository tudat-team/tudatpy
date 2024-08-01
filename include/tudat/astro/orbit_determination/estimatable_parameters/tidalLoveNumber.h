/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"

using std::string;

namespace tudat
{

namespace estimatable_parameters
{

//! List of all orders in a single degree (starting at 2)
/*!
 *  List of all orders in a single degree (starting at 2), generated at compile time and used for tidal love number
 *  partial interface function calls over a full degree.
 */
static const std::vector< std::vector< int > > fullDegreeOrders =
{ { 0, 1, 2 }, { 0, 1, 2, 3 }, { 0, 1, 2, 3, 4 }, { 0, 1, 2, 3, 4, 5 }, { 0, 1, 2, 3, 4, 5, 6 } };


//!Pure virtual base class for estimating tidal Love number properties
template< typename ParameterScalar >
class TidalLoveNumber: public EstimatableParameter< ParameterScalar >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param gravityFieldVariationModel Tidal gravity field variation object of which estimated paraemeter is a property
     * \param associatedBody Deformed body
     * \param degree Degree of Love number that is to be estimated
     * \param orders List of orders at which Love numbers are to be estimated.
     * \param sumOrders True of the contributions of the various orders are to be summed (i.e. assumed to be constant for all
     * orders at given degree), or if they are handled separately
     * \param loveNumberType Type of Love number property that is to be estimated.
     * \param useComplexComponents True if the complex Love number is estimated, false if only the real part is considered
     */
    TidalLoveNumber(
            const std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariationModel,
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

    //! Destructor
    virtual ~TidalLoveNumber( ){ }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value
     */
    int getParameterSize( )
    {
        return parameterSize_;
    }

    //! Function to retrieve degree of Love number that is to be estimated
    /*!
     * Function to retrieve degree of Love number that is to be estimated
     * \return Degree of Love number that is to be estimated
     */
    int getDegree( )
    {
        return degree_;
    }

    //! Function to retrieve whether complex or true Love number(s) is(are) to be estimated
    /*!
     * Function to retrieve whether complex or true Love number(s) is(are) to be estimated
     * \return True if the complex Love number is estimated, false if only the real part is considered
     */
    bool useComplexComponents( )
    {
        return useComplexComponents_;
    }

    //! Function to retrieve list of bodies causing tidal deformation
    /*!
     * Function to retrieve list of bodies causing tidal deformation
     * \return :ist of bodies causing tidal deformation
     */
    std::vector< std::string > getDeformingBodies( )
    {
        return gravityFieldVariationModel_->getDeformingBodies( );
    }

    //! Function to retrieve list of orders at which Love numbers are to be estimated.
    /*!
     * Function to retrieve list of orders at which Love numbers are to be estimated.
     * \return List of orders at which Love numbers are to be estimated.
     */
    std::vector< int > getOrders( )
    {
        return orders_;
    }

    //! Function to retrieve whether contributions of the various orders are to be summed
    /*!
     * Function to retrieve contributions of the various orders are to be summed
     * \return True of the contributions of the various orders are to be summed, false if not
     */
    bool getSumOrders( )
    {
        return sumOrders_;
    }

    //! Function to retrieve parameter description.
    /*!
    * Function to retrieve parameter description.
    * \return Description tidal Love number.
    */
    std::string getParameterDescription( )
    {
        std::string parameterDescription =
                getParameterTypeString( this->parameterName_.first ) + "of " + this->parameterName_.second.first + " due to ";
        for ( unsigned int i = 0 ; i < getDeformingBodies( ).size( ) ; i++ )
        {
            parameterDescription += getDeformingBodies( )[ i ];
            if ( i != getDeformingBodies( ).size( ) - 1 )
            {
                parameterDescription += " & ";
            }
        }
        parameterDescription += " for degree " + std::to_string(degree_) + " and orders: ";
        for ( unsigned int i = 0 ; i < orders_.size( ) ; i++ )
        {
            if (i != orders_.size( )-1)
            {
                parameterDescription += std::to_string(orders_[i]) + ", ";
            }
            else
            {
                parameterDescription += std::to_string(orders_[i]) + ".";
            }
        }
        return parameterDescription;
    }

protected:

    //! Degree of Love number that is to be estimated
    int degree_;

    //! List of orders at which Love numbers are to be estimated.
    std::vector< int > orders_;

    //! True of the contributions of the various orders are to be summed
    /*!
     *  True of the contributions of the various orders are to be summed (i.e. assumed to be constant for all
     *  orders at given degree), or if they are handled separately
     */
    bool sumOrders_;

    //! Tidal gravity field variation object of which estimated paraemeter is a property
    std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariationModel_;

    //! True if the complex Love number is estimated, false if only the real part is considered
    bool useComplexComponents_;

    //! Size of the estimated parameter
    int parameterSize_;

};

//! Class for estimating the tidal Love number k_{n} at a single degree that is constant for all orders
class FullDegreeTidalLoveNumber: public TidalLoveNumber< Eigen::VectorXd >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param gravityFieldVariationModel Tidal gravity field variation object of which estimated paraemeter is a property
     * \param associatedBody Deformed body
     * \param degree Degree of Love number that is to be estimateds
     * \param useComplexComponents True if the complex Love number is estimated, false if only the real part is considered
     */
    FullDegreeTidalLoveNumber(
            const std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariationModel,
            const std::string& associatedBody,
            const int degree,
            const bool useComplexComponents = 0 ):
        TidalLoveNumber< Eigen::VectorXd >(
            gravityFieldVariationModel, associatedBody, degree, fullDegreeOrders.at( degree - 2 ), 1,
            full_degree_tidal_love_number, useComplexComponents )
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


    //! Get value of Love number k_{n}
    /*!
     *  Get value of Love number k_{n} (complex part as second vector entry if complex value is used)
     *  \return Value of Love number k_{n} (complex part as second vector entry if complex value is used)
     */
    Eigen::VectorXd getParameterValue( );

    //! Reset value of Love number k_{n}
    /*!
     *  Reset value of Love number k_{n} (complex part as second vector entry if complex value is used)
     *  \param parameterValue New value of Love number k_{n} (complex part as second vector entry if complex value is used)
     */
    void setParameterValue( Eigen::VectorXd parameterValue );

    //! Function to retrieve parameter description.
    /*!
    * Function to retrieve parameter description.
    * \return Description tidal Love number.
    */
    std::string getParameterDescription( )
    {
        std::string parameterDescription =
                getParameterTypeString( this->parameterName_.first ) + "of " + this->parameterName_.second.first + " due to ";
        for ( unsigned int i = 0 ; i < getDeformingBodies( ).size( ) ; i++ )
        {
            parameterDescription += getDeformingBodies( )[ i ];
            if ( i != getDeformingBodies( ).size( ) - 1 )
            {
                parameterDescription += " & ";
            }
        }
        parameterDescription += " for degree " + std::to_string(degree_) + ".";
        return parameterDescription;
    }

};


//! Class for estimating the tidal Love number k_{n} at a single degree that is constant for all orders
class ModeCoupledTidalLoveNumber: public EstimatableParameter< Eigen::VectorXd >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param gravityFieldVariationModel Tidal gravity field variation object of which estimated paraemeter is a property
     * \param associatedBody Deformed body
     * \param degree Degree of Love number that is to be estimateds
     * \param useComplexComponents True if the complex Love number is estimated, false if only the real part is considered
     */
    ModeCoupledTidalLoveNumber(
        const std::shared_ptr< gravitation::ModeCoupledSolidBodyTideGravityFieldVariations > gravityFieldVariationModel,
        const std::string& associatedBody,
        const std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > loveNumberIndices,
        const bool useComplexComponents = 0 );

    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd currentParameters = Eigen::VectorXd::Zero( parameterSize_ );
        std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > loveNumbers = gravityFieldVariationModel_->getLoveNumbers( );

        int counter = 0;
        for( auto it: loveNumberIndices_ )
        {
            for( unsigned int i = 0; i < it.second.size( ); i++ )
            {
                currentParameters( counter ) = loveNumbers.at( it.first ).at( it.second.at( i ) );
                counter++;
            }
        }
        return currentParameters;
    }


    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > loveNumbers = gravityFieldVariationModel_->getLoveNumbers( );

        int counter = 0;
        for( auto it: loveNumberIndices_ )
        {
            for( unsigned int i = 0; i < it.second.size( ); i++ )
            {
                gravityFieldVariationModel_->resetLoveNumber( it.first, it.second.at( i ), parameterValue( counter ) );
                counter++;
            }
        }
    }

    int getParameterSize( )
    {
        return parameterSize_;
    }

    int getMaximumForcingDegree( )
    {
        return maximumForcingDegree_;
    }

    std::vector< std::string > getDeformingBodies( )
    {
        return gravityFieldVariationModel_->getDeformingBodies( );
    }


    std::vector< std::pair< int, int > > getParameterForcingDegreeAndOrderIndices( )
    {
        return parameterForcingDegreeAndOrderIndices_;
    }

    std::map< int, std::vector< int > > getForcingOrdersPerDegree( )
    {
        return forcingOrdersPerDegree_;
    }

    std::vector< std::pair< int, int > > getResponseDegreeOrders( )
    {
        return responseDegreeOrders_;
    }

    std::vector< int > getResponseIndices( )
    {
        return responseIndices_;
    }

private:
    std::shared_ptr< gravitation::ModeCoupledSolidBodyTideGravityFieldVariations > gravityFieldVariationModel_;

    // Map with forcing degree/orders (key) and list of associated response degree orders (values)
    std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > loveNumberIndices_;

    // NUmber of separate Love numbers
    int parameterSize_;

    // Maximum degree of forcing
    int maximumForcingDegree_;

    // The list of forcing orders per forcing degree
    std::map< int, std::vector< int > > forcingOrdersPerDegree_;

    // For each independent Love number, the forcing degree, and the index in the associated vector of orders in forcingOrdersPerDegree_
    std::vector< std::pair< int, int > > parameterForcingDegreeAndOrderIndices_;

    // List of all response degree/orders
    std::vector< std::pair< int, int > > responseDegreeOrders_;

    // For each independent Love number, the response degree/order (as an index in responseDegreeOrders_);
    std::vector< int > responseIndices_;



};


//! Class for estimating the tidal Love numbers k_{n,m} at a single degree that may vary for different orders
class SingleDegreeVariableTidalLoveNumber: public TidalLoveNumber< Eigen::VectorXd >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param gravityFieldVariationModel Tidal gravity field variation object of which estimated paraemeter is a property
     * \param associatedBody Deformed body
     * \param degree Degree of Love number that is to be estimated
     * \param orders List of orders at which Love numbers are to be estimated.
     * \param useComplexComponents True if the complex Love number is estimated, false if only the real part is considered
     */
    SingleDegreeVariableTidalLoveNumber(
            const std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariationModel,
            const std::string& associatedBody,
            const int degree,
            const std::vector< int > orders,
            const bool useComplexComponents = 0 ):
        TidalLoveNumber< Eigen::VectorXd >(
            gravityFieldVariationModel, associatedBody, degree, orders, 0, single_degree_variable_tidal_love_number,
            useComplexComponents )
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

    //! Get values of Love numbers k_{n,m}
    /*!
     *  Get values of Love numbers k_{n,m}, concatenated by order (complex part as second entry of single k_{n,m} if
     * complex value is used)
     *  \return Value of Love number k_{n} (complex part as second vector entry if complex value is used)
     */
    Eigen::VectorXd getParameterValue( );

    //! Reset value of Love number k_{n,m}
    /*!
     *  Reset value of Love number k_{n,m}, concatenated by order (complex part as second entry of single k_{n,m} if complex
     *  value is used)
     *  \param parameterValue New values of Love numbers k_{n,m}
     */
    void setParameterValue( Eigen::VectorXd parameterValue );

    //! Function to retrieve parameter description.
    /*!
    * Function to retrieve parameter description.
    * \return Description tidal Love number.
    */
    std::string getParameterDescription( )
    {
        std::string parameterDescription =
                getParameterTypeString( this->parameterName_.first ) + "of " + this->parameterName_.second.first + " due to ";
        for ( unsigned int i = 0 ; i < getDeformingBodies( ).size( ) ; i++ )
        {
            parameterDescription += getDeformingBodies( )[ i ];
            if ( i != getDeformingBodies( ).size( ) - 1 )
            {
                parameterDescription += " & ";
            }
        }
        parameterDescription += " for degree " + std::to_string(degree_) + " and orders: ";
        for ( unsigned int i = 0 ; i < orders_.size( ) ; i++ )
        {
            if (i != orders_.size( )-1)
            {
                parameterDescription += std::to_string(orders_[i]) + ", ";
            }
            else
            {
                parameterDescription += std::to_string(orders_[i]) + ".";
            }

        }
        return parameterDescription;
    }
};

}

}

#endif // TUDAT_TIDALLOVENUMBERPARAMETER_H
