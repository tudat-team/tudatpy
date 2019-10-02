/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef FULLPLANETARYROTATIONMODEL_H
#define FULLPLANETARYROTATIONMODEL_H

#include <vector>
#include <map>

#include <boost/function.hpp>

#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{
    
    namespace ephemerides
    {
        
        class PlanetaryOrientationAngleCalculator
        {
        public:
            PlanetaryOrientationAngleCalculator( )
            { }
            
            PlanetaryOrientationAngleCalculator(
                    const double anglePsiAtEpoch, const double anglePsiRateAtEpoch, const double angleIAtEpoch,
                    const double angleIRateAtEpoch, const double anglePhiAtEpoch, const double anglePhiRateAtEpoch,
                    const double coreFactor, const double freeCoreNutationRate,
                    const double bodyMeanMotion, const double bodyMeanAnomalyAtEpoch, const std::string baseFrame,

                    const std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections =
                    ( std::map< double, std::pair< double, double > >( ) ), // Konopliv table 5, 1st four terms

                    const std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections =
                    ( std::vector< std::map< double, std::pair< double, double > > > ( ) ), // Konopliv table 5, last six terms

                    const std::vector< std::function< double( const double ) > > phaseAngleCorrectionFunctions =
                    ( std::vector< std::function< double( const double ) > > ( ) ),

                    const std::map< double, std::pair< double, double > > rotationRateCorrections =
                    ( std::map< double, std::pair< double, double > >( ) ),  // Konopliv table 7, 1st four terms

                    const std::map< double, std::pair< double, double > > xPolarMotionCoefficients =
                    ( std::map< double, std::pair< double, double > >( ) ),

                    const std::map< double, std::pair< double, double > > yPolarMotionCoefficients =
                    ( std::map< double, std::pair< double, double > >( ) ) ):

                bodyMeanAnomalyAtEpoch_( bodyMeanAnomalyAtEpoch),
                bodyMeanMotion_( bodyMeanMotion ),
                anglePsiAtEpoch_( anglePsiAtEpoch ),
                anglePsiRateAtEpoch_( anglePsiRateAtEpoch ),
                angleIAtEpoch_( angleIAtEpoch ),
                angleIRateAtEpoch_( angleIRateAtEpoch ),
                anglePhiAtEpoch_( anglePhiAtEpoch ),
                anglePhiRateAtEpoch_( anglePhiRateAtEpoch ),
                coreFactor_ ( coreFactor ),
                freeCoreNutationRate_( freeCoreNutationRate ),
                baseFrame_( baseFrame ),
                meanMotionDirectNutationCorrections_( meanMotionDirectNutationCorrections ),
                meanMotionTimeDependentPhaseNutationCorrections_( meanMotionTimeDependentPhaseNutationCorrections ),
                phaseAngleCorrectionFunctions_( phaseAngleCorrectionFunctions ),
                rotationRateCorrections_( rotationRateCorrections ),
                xPolarMotionCoefficients_( xPolarMotionCoefficients ),
                yPolarMotionCoefficients_( yPolarMotionCoefficients ){ }
            
            Eigen::Vector3d updateAndGetRotationAngles( const double ephemerisTime )
            {
                //if( std::fabs( currentEphemerisTime_ - ephemerisTime ) > 1.0E-8 )
                {
                    updateCorrections( ephemerisTime );
                }
                return ( Eigen::Vector3d( ) << currentAnglePsi_, currentAngleI_, currentAnglePhi_ ).finished( );
            }
            
            Eigen::Quaterniond getRotationMatrixFromMeanOfDateEquatorToInertialPlanetCenteredAtEpoch( const double ephemerisTime )
            {
                return Eigen::Quaterniond(
                            Eigen::AngleAxisd( anglePsiAtEpoch_ + ephemerisTime * anglePsiRateAtEpoch_, Eigen::Vector3d::UnitZ( ) )  *
                            Eigen::AngleAxisd( angleIAtEpoch_ + ephemerisTime * angleIRateAtEpoch_, Eigen::Vector3d::UnitX( ) ) );
            }
            
            double getCurrentAnglePsi( )
            {
                return currentAnglePsi_;
            }
            
            double getCurrentAngleI( )
            {
                return currentAngleI_;
            }
            
            double getCurrentAnglePhi( )
            {
                return currentAnglePhi_;
            }
            
            double getMeanPhiAngleDerivative( )
            {
                return anglePhiRateAtEpoch_;
            }

            double getCorefactor()
            {
                return coreFactor_;
            }

            double getFreeCoreNutationRate()
            {
                return freeCoreNutationRate_;
            }
            
            std::string getBaseFrame( )
            {
                return baseFrame_;
            }
            
            double getAnglePsiRateAtEpoch( )
            {
                return anglePsiRateAtEpoch_;
            }
            
            void setAnglePsiRateAtEpoch( const double anglePsiRateAtEpoch )
            {
                anglePsiRateAtEpoch_ = anglePsiRateAtEpoch;
                currentEphemerisTime_ = -1.0E100;
            }

            Eigen::Vector2d getPolarMotion( const double ephemerisTime )
            {
                    calculatePolarMotion( ephemerisTime );                  
                    return ( Eigen::Vector2d( ) << currentXPolarMotion_ , currentYPolarMotion_ ).finished( );
            }

            double getcurrentMeanPhiAngleDerivative( const double ephemerisTime )
            {
                calculateCurrentMeanPhiAngleDerivative( ephemerisTime );
                return currentMeanPhiAngleDerivative_;
            }

            std::map< double, std::pair< double, double > > getRotationRateCorrections( )
            {
                return rotationRateCorrections_;
            }

            std::map< double, std::pair< double, double > > getMeanMotionDirectNutationCorrections( )
            {
                return meanMotionDirectNutationCorrections_;
            }

            std::vector< std::map< double, std::pair< double, double > > > getMeanMotionTimeDependentPhaseNutationCorrections( )
            {
                return meanMotionTimeDependentPhaseNutationCorrections_;
            }

            std::vector< std::function< double( const double ) > > getphaseAngleCorrectionFunctions ( )
            {
                return phaseAngleCorrectionFunctions_;
            }

            double getBodyMeanMotion( )
            {
                return bodyMeanMotion_;
            }

            double getBodyMeanAnomalyAtEpoch( )
            {
                return bodyMeanAnomalyAtEpoch_;
            }

            double getAngleIAtEpoch( )
            {
                return angleIAtEpoch_;
            }

            void resetRotationRateCorrections( std::map< double, std::pair< double, double > > rotationRateCorrections ){
                rotationRateCorrections_ = rotationRateCorrections;
            }

            std::map< double, std::pair< double, double > > getXpolarMotionCoefficients( )
            {
                return xPolarMotionCoefficients_;
            }

            std::map< double, std::pair< double, double > > getYpolarMotionCoefficients( )
            {
                return yPolarMotionCoefficients_;
            }

            void resetXpolarMotionCoefficients( std::map< double, std::pair< double, double > > xPolarMotionCoefficients ){
                xPolarMotionCoefficients_ = xPolarMotionCoefficients;
            }

            void resetYpolarMotionCoefficients( std::map< double, std::pair< double, double > > yPolarMotionCoefficients ){
                yPolarMotionCoefficients_ = yPolarMotionCoefficients;
            }

            void resetCoreFactor( const double coreFactor ){
                coreFactor_ = coreFactor;
            }

            void resetFreeCoreNutationRate( const double freeCoreNutationRate ){
                freeCoreNutationRate_ = freeCoreNutationRate;
            }



        private:
            
            void updateCorrections( const double ephemerisTime );
            void calculatePolarMotion ( const double ephemerisTime );
            void calculateCurrentMeanPhiAngleDerivative( const double ephemerisTime );
            
            double currentEphemerisTime_;
            
            double bodyMeanAnomalyAtEpoch_;
            double bodyMeanMotion_;
            
            double currentAnglePsiCorrection_;
            double currentAngleICorrection_;
            double currentAnglePhiCorrection_;
            
            double currentAnglePsi_;
            double currentAngleI_;
            double currentAnglePhi_;
            double currentXPolarMotion_;
            double currentYPolarMotion_;
            double currentMeanPhiAngleDerivative_;
            
            double anglePsiAtEpoch_;
            double anglePsiRateAtEpoch_;
            double angleIAtEpoch_;
            double angleIRateAtEpoch_;
            double anglePhiAtEpoch_;
            double anglePhiRateAtEpoch_;
            double coreFactor_;
            double freeCoreNutationRate_;
            
            std::string baseFrame_;
            
            std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections_;
            
            std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections_;
            
            std::vector< std::function< double( const double ) > > phaseAngleCorrectionFunctions_;
            
            std::map< double, std::pair< double, double > > rotationRateCorrections_;

            std::map< double, std::pair< double, double > > xPolarMotionCoefficients_;

            std::map< double, std::pair< double, double > > yPolarMotionCoefficients_;

        };
        
        std::shared_ptr< interpolators::CubicSplineInterpolator< double, Eigen::Vector3d > >
        createInterpolatorForPlanetaryRotationAngles( double intervalStart,
                                                     double intervalEnd,
                                                     double timeStep,
                                                     std::shared_ptr< PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator );
        
        class PlanetaryRotationModel: public RotationalEphemeris
        {
        public:
            PlanetaryRotationModel( const double angleN,
                                   const double angleJ,
                                   const std::shared_ptr< PlanetaryOrientationAngleCalculator > planetaryOrientationAnglesCalculator,
                                   const std::string& baseFrameOrientation = "",
                                   const std::string& targetFrameOrientation = "" ):

            RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
            planetaryOrientationAnglesCalculator_( planetaryOrientationAnglesCalculator )

            {
                rotationFromMeanOrbitToIcrf_ = spice_interface::computeRotationQuaternionBetweenFrames( "J2000", "ECLIPJ2000", 0.0 ) *
                Eigen::AngleAxisd( angleN, Eigen::Vector3d::UnitZ( ) ) *
                Eigen::AngleAxisd( angleJ, Eigen::Vector3d::UnitX( ) );
            }
            
            Eigen::Quaterniond getPolarMotionRotation( const double ephemerisTime );

            Eigen::Quaterniond getRotationFromBodyFixedToIntermediateInertialFrame( const double ephemerisTime );
            
            Eigen::Quaterniond getRotationToBaseFrame( const double ephemerisTime );
            
            Eigen::Quaterniond getRotationToTargetFrame( const double ephemerisTime )
            {
                return getRotationToBaseFrame( ephemerisTime ).inverse( );
            }
            
            Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double ephemerisTime );
            
            Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( const double secondsSinceEpoch )
            {
                return getDerivativeOfRotationToBaseFrame( secondsSinceEpoch ).
                transpose( );
            }
            
            Eigen::Quaterniond getRotationFromMeanOrbitToIcrf( )
            {
                return rotationFromMeanOrbitToIcrf_;
            }

            std::shared_ptr< PlanetaryOrientationAngleCalculator > getPlanetaryOrientationAngleCalculator()
            {
                return planetaryOrientationAnglesCalculator_;
            }

        private:
            
            std::shared_ptr< PlanetaryOrientationAngleCalculator > planetaryOrientationAnglesCalculator_;
            
            Eigen::Quaterniond rotationFromMeanOrbitToIcrf_;
                        
        };
        
    }
    
}

#endif // FULLPLANETARYROTATIONMODEL_H
