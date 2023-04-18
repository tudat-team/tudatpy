/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TABULATEDMEDIACORRECTION_H
#define TUDAT_TABULATEDMEDIACORRECTION_H

namespace tudat
{

namespace observation_models
{

class TabulatedMediaReferenceCorrection
{
public:

    TabulatedMediaReferenceCorrection( const double startTime,
                                       const double endTime ):
       startTime_( startTime ),
       endTime_( endTime )
    { }

    virtual double computeReferenceCorrection( const double time ) = 0;

protected:

    bool isTimeValid( const double time )
    {
        if ( time < startTime_ || time > endTime_ )
        {
            throw std::runtime_error(
                    "Error when computing tabulated media reference correction: selected time (" + std::to_string( time ) +
                    ") is outside validity interval (" + std::to_string( startTime_ ) + " to " + std::to_string( endTime_ ) +
                    ")." );
        }

        return true;
    }

    const double startTime_;
    const double endTime_;

private:

};

class ConstantReferenceCorrection: protected TabulatedMediaReferenceCorrection
{
public:

    ConstantReferenceCorrection( const double startTime,
                                 const double endTime,
                                 const double constantCorrection ):
        TabulatedMediaReferenceCorrection( startTime, endTime ),
        constantCorrection_( constantCorrection )
    { }

    double computeReferenceCorrection( const double time ) override
    {
        isTimeValid( time );
        
        return constantCorrection_;
    }

private:

    const double constantCorrection_;

};

class PowerSeriesReferenceCorrection: protected TabulatedMediaReferenceCorrection
{
public:

    PowerSeriesReferenceCorrection( const double startTime,
                                    const double endTime,
                                    const std::vector< double > coefficients ):
        TabulatedMediaReferenceCorrection( startTime, endTime ),
        coefficients_( coefficients )
    { }

    double computeReferenceCorrection( const double time ) override
    {
        isTimeValid( time );

        const double normalizedTime = 2.0 * ( ( time - startTime_ ) / ( time - endTime_ ) ) - 1.0;

        double correction = 0;
        for ( unsigned int i = 0; i < coefficients_.size( ); ++i )
        {
            correction += coefficients_.at( i ) * std::pow( normalizedTime, i );
        }

        return correction;
    }

private:

    const std::vector< double > coefficients_;

};

class FourierSeriesReferenceCorrection: protected TabulatedMediaReferenceCorrection
{
public:

    FourierSeriesReferenceCorrection( const double startTime,
                                      const double endTime,
                                      const std::vector< double > coefficients ):
        TabulatedMediaReferenceCorrection( startTime, endTime )
    {
        if ( coefficients.size( ) < 2 || coefficients.size( ) % 2 != 0 )
        {
            throw std::runtime_error(
                    "Error when computing Fourier series tabulated media reference correction: size of specified coefficients ("
                    + std::to_string( coefficients.size( ) ) + ") is invalid." );
        }

        period_ = coefficients.at( 0 );

        cosineCoefficients_.push_back( coefficients.at( 1 ) );
        sineCoefficients_.push_back( 0.0 );

        for ( unsigned int i = 2; i < coefficients.size( ); i = i + 2 )
        {
            cosineCoefficients_.push_back( coefficients.at( i ) );
            sineCoefficients_.push_back( coefficients.at( i + 1 ) );
        }
    }

    double computeReferenceCorrection( const double time ) override
    {
        isTimeValid( time );

        const double normalizedTime = 2.0 * mathematical_constants::PI * ( time - startTime_ ) / period_;

        double correction = 0;
        for ( unsigned int i = 0; i < sineCoefficients_.size( ); ++i )
        {
            correction += cosineCoefficients_.at( i ) * std::cos( i * normalizedTime );
            correction += sineCoefficients_.at( i ) * std::sin( i * normalizedTime );
        }

        return correction;
    }

private:

    double period_;

    std::vector< double > sineCoefficients_;
    std::vector< double > cosineCoefficients_;

};

} // namespace observation_models

} // namespace tudat

#endif //TUDATBUNDLE_TABULATEDMEDIACORRECTION_H
