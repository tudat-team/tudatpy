/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <vector>
#include <map>
#include <tuple>

#include "Tudat/Astrodynamics/Ephemerides/fullPlanetaryRotationModel.h"

namespace tudat
{
    
    namespace ephemerides
    {
        
        void PlanetaryOrientationAngleCalculator::updateCorrections( const double ephemerisTime )
        {
            currentEphemerisTime_ = ephemerisTime;
            double currentMeanMotion = bodyMeanMotion_;
            double currentMeanAnomaly = bodyMeanAnomalyAtEpoch_ + bodyMeanMotion_ * currentEphemerisTime_;
            
            currentAngleICorrection_ = 0.0;
            currentAnglePsiCorrection_ = 0.0;
            
            for( std::map< double, std::pair< double, double > >::iterator correctionIterator = meanMotionDirectNutationCorrections_.begin( );
                correctionIterator != meanMotionDirectNutationCorrections_.end( ); correctionIterator++ )
            {
                double a_m = correctionIterator->first * currentMeanMotion;

                double I_m = correctionIterator->second.first + coreFactor_ * a_m /
                        ( a_m * a_m - freeCoreNutationRate_ * freeCoreNutationRate_ ) *
                        ( a_m * correctionIterator->second.first + freeCoreNutationRate_ *
                          correctionIterator->second.second * std::sin( angleIAtEpoch_ ) );

                double Psi_m = correctionIterator->second.second + coreFactor_ * a_m /
                        ( a_m * a_m - freeCoreNutationRate_ * freeCoreNutationRate_ ) *
                        ( a_m * correctionIterator->second.second + freeCoreNutationRate_ *
                          correctionIterator->second.first / std::sin( angleIAtEpoch_ ) );

                currentAngleICorrection_ += I_m * std::cos(
                            correctionIterator->first * ( currentMeanMotion * currentEphemerisTime_ + bodyMeanAnomalyAtEpoch_ )  );
                currentAnglePsiCorrection_ += Psi_m * std::sin(
                            correctionIterator->first * ( currentMeanMotion * currentEphemerisTime_ + bodyMeanAnomalyAtEpoch_ )  );
            }
            
            double currentPhaseAngleCorrection;
            
            for( unsigned int i = 0; i < meanMotionTimeDependentPhaseNutationCorrections_.size( ); i++ )
            {
                currentPhaseAngleCorrection = phaseAngleCorrectionFunctions_[ i ]( currentEphemerisTime_ );
                for( std::map< double, std::pair< double, double > >::iterator correctionIterator =
                    meanMotionTimeDependentPhaseNutationCorrections_[ i ].begin( );
                    correctionIterator != meanMotionTimeDependentPhaseNutationCorrections_[ i ].end( ); correctionIterator++ )
                {
                    double a_m = correctionIterator->first * currentMeanMotion;

                    double I_m = correctionIterator->second.first + coreFactor_ * a_m /
                            ( a_m * a_m - freeCoreNutationRate_ * freeCoreNutationRate_ ) *
                            ( a_m * correctionIterator->second.first + freeCoreNutationRate_ *
                              correctionIterator->second.second * std::sin( angleIAtEpoch_ ) );

                    double Psi_m = correctionIterator->second.second + coreFactor_ * a_m /
                            ( a_m * a_m - freeCoreNutationRate_ * freeCoreNutationRate_ ) *
                            ( a_m * correctionIterator->second.second + freeCoreNutationRate_ *
                              correctionIterator->second.first / std::sin( angleIAtEpoch_ ) );


                    currentAngleICorrection_ += I_m * std::cos(
                                correctionIterator->first * ( currentMeanMotion * currentEphemerisTime_ + bodyMeanAnomalyAtEpoch_ ) +
                                currentPhaseAngleCorrection );
                    currentAnglePsiCorrection_ += Psi_m * std::sin(
                                correctionIterator->first * ( currentMeanMotion * currentEphemerisTime_ + bodyMeanAnomalyAtEpoch_ ) +
                                currentPhaseAngleCorrection );
                }
            }
            
            currentAnglePsi_ = anglePsiAtEpoch_ + anglePsiRateAtEpoch_ * currentEphemerisTime_ + currentAnglePsiCorrection_;

            currentAngleI_ = angleIAtEpoch_ + angleIRateAtEpoch_ * currentEphemerisTime_ + currentAngleICorrection_;
            
            currentAnglePhiCorrection_ = -currentAnglePsiCorrection_ * std::cos( currentAngleI_ );
            
            for( std::map< double, std::pair< double, double > >::iterator correctionIterator = rotationRateCorrections_.begin( );
                correctionIterator != rotationRateCorrections_.end( ); correctionIterator++ )
            {
                currentAnglePhiCorrection_ += correctionIterator->second.first * std::cos(
                            correctionIterator->first * currentMeanAnomaly );
                currentAnglePhiCorrection_ += correctionIterator->second.second * std::sin(
                            correctionIterator->first * currentMeanAnomaly );
            }

            currentAnglePhi_ = anglePhiAtEpoch_ + anglePhiRateAtEpoch_ * currentEphemerisTime_ + currentAnglePhiCorrection_;
        }
        
        
        void PlanetaryOrientationAngleCalculator::calculatePolarMotion( const double ephemerisTime )
        {
            currentEphemerisTime_ = ephemerisTime;
            double currentMeanAnomaly = bodyMeanAnomalyAtEpoch_ + bodyMeanMotion_ * currentEphemerisTime_;

            currentXPolarMotion_ = 0;
            for( std::map< double, std::pair< double, double > >::iterator polarMotionIterator = xPolarMotionCoefficients_.begin( );
                polarMotionIterator != xPolarMotionCoefficients_.end( ); polarMotionIterator++ )
            {
                    currentXPolarMotion_ += polarMotionIterator->second.first *
                            std::cos( polarMotionIterator->first * currentMeanAnomaly );
                    currentXPolarMotion_ += polarMotionIterator->second.second *
                            std::sin( polarMotionIterator->first * currentMeanAnomaly );
            }

            currentYPolarMotion_ = 0;
            for( std::map< double, std::pair< double, double > >::iterator polarMotionIterator = yPolarMotionCoefficients_.begin( );
                polarMotionIterator != yPolarMotionCoefficients_.end( ); polarMotionIterator++ )
            {
                    currentYPolarMotion_ += polarMotionIterator->second.first *
                            std::cos( polarMotionIterator->first * currentMeanAnomaly );
                    currentYPolarMotion_ += polarMotionIterator->second.second *
                            std::sin( polarMotionIterator->first * currentMeanAnomaly );
            }

        }

        void PlanetaryOrientationAngleCalculator::calculateCurrentMeanPhiAngleDerivative( const double ephemerisTime )
        {
            currentEphemerisTime_ = ephemerisTime;
            double currentMeanMotion = bodyMeanMotion_;
            double currentMeanAnomaly = bodyMeanAnomalyAtEpoch_ + bodyMeanMotion_ * currentEphemerisTime_;

            currentMeanPhiAngleDerivative_ = anglePhiRateAtEpoch_ ;

            for( std::map< double, std::pair< double, double > >::iterator correctionIterator = rotationRateCorrections_.begin( );
                correctionIterator != rotationRateCorrections_.end( ); correctionIterator++ )
            {
                currentMeanPhiAngleDerivative_ += - correctionIterator->first * currentMeanMotion * correctionIterator->second.first * std::sin(
                            correctionIterator->first * currentMeanAnomaly );
                currentMeanPhiAngleDerivative_ += correctionIterator->first * currentMeanMotion * correctionIterator->second.second * std::cos(
                            correctionIterator->first * currentMeanAnomaly );
            }
        }

        boost::shared_ptr< interpolators::CubicSplineInterpolator< double, Eigen::Vector3d > >
        createInterpolatorForPlanetaryRotationAngles( double intervalStart,
                                                     double intervalEnd,
                                                     double timeStep,
                                                     boost::shared_ptr< PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator )
        {
            using namespace interpolators;
            
            std::map< double, Eigen::Vector3d > orientationMap;
            
            double currentTime = intervalStart;
            while( currentTime < intervalEnd )
            {
                orientationMap[ currentTime ] = planetaryOrientationCalculator->updateAndGetRotationAngles( currentTime  );

                currentTime += timeStep;
            }
            
            boost::shared_ptr< CubicSplineInterpolator< double, Eigen::Vector3d > > interpolator =
            boost::make_shared< CubicSplineInterpolator< double, Eigen::Vector3d > >( orientationMap );
            return interpolator;
        }
        
        Eigen::Quaterniond PlanetaryRotationModel::getPolarMotionRotation( const double ephemerisTime )
        {
            Eigen::Vector2d currentPolarMotion = planetaryOrientationAnglesCalculator_->getPolarMotion( ephemerisTime );

            return Eigen::Quaterniond( Eigen::AngleAxisd( -currentPolarMotion.x( ), Eigen::Vector3d::UnitY( ) ) *
                                       Eigen::AngleAxisd( -currentPolarMotion.y( ), Eigen::Vector3d::UnitX( ) ) );
        }

        Eigen::Quaterniond PlanetaryRotationModel::getRotationFromBodyFixedToIntermediateInertialFrame( const double ephemerisTime )
        {
            Eigen::Vector3d currentAngleCorrections = planetaryOrientationAnglesCalculator_->updateAndGetRotationAngles( ephemerisTime );

            return Eigen::Quaterniond( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                                       Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) *
                                       Eigen::AngleAxisd( currentAngleCorrections.z( ), Eigen::Vector3d::UnitZ( ) ) );

        }
        
        Eigen::Quaterniond PlanetaryRotationModel::getRotationToBaseFrame( const double ephemerisTime )
        {
            return rotationFromMeanOrbitToIcrf_ * getRotationFromBodyFixedToIntermediateInertialFrame( ephemerisTime ) *
                    getPolarMotionRotation( ephemerisTime );
        }
        
        Eigen::Matrix3d PlanetaryRotationModel::getDerivativeOfRotationToBaseFrame( const double ephemerisTime )
        {
            Eigen::Vector3d currentAngleCorrections = planetaryOrientationAnglesCalculator_->updateAndGetRotationAngles( ephemerisTime );

            double currentPhiAngle = currentAngleCorrections.z( );

            double meanPhiAngleDerivative = planetaryOrientationAnglesCalculator_->getcurrentMeanPhiAngleDerivative( ephemerisTime );
              
            return meanPhiAngleDerivative * ( rotationFromMeanOrbitToIcrf_ ).toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentPhiAngle ), -std::cos( currentPhiAngle ), 0.0,
                      std::cos( currentPhiAngle ), -std::sin( currentPhiAngle ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    getPolarMotionRotation( ephemerisTime ).toRotationMatrix();

        }

    }
    
}
