/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_KEYS_H
#define TUDAT_JSONINTERFACE_KEYS_H

#include <string>
#include <vector>
#include <iostream>

namespace tudat
{

namespace json_interface
{

//! Special keys (used internally by json_interface, can't be used in JSON files).
struct SpecialKeys
{
    static const std::string root;
    static const std::string up;
    static const std::string rootObject;
    static const std::string accessHistory;
    static const std::string keyPath;

    static const std::vector< std::string > all;
};

//! Keys recognised by json_interface.
struct Keys
{
    static const std::string simulationType;
    static const std::string startEpoch;
    static const std::string endEpoch;
    static const std::string globalFrameOrigin;
    static const std::string globalFrameOrientation;

    static const std::string spice;
    struct Spice
    {
        static const std::string kernels;
        static const std::string preloadKernels;
        static const std::string preloadOffsets;
        static const std::string preloadOffset;
        static const std::string startPreloadOffset;
        static const std::string endPreloadOffset;
    };

    static const std::string bodies;
    struct Body
    {
        static const std::string useDefaultSettings;
        static const std::string mass;
        static const std::string referenceArea;

        static const std::string aerodynamics;
        struct Aerodynamics
        {
            static const std::string type;
            static const std::string referenceArea;
            static const std::string dragCoefficient;
            static const std::string forceCoefficients;
            static const std::string momentCoefficients;
            static const std::string areCoefficientsInAerodynamicFrame;
            static const std::string areCoefficientsInNegativeAxisDirection;
        };

        static const std::string atmosphere;
        struct Atmosphere
        {
            static const std::string type;
            static const std::string densityScaleHeight;
            static const std::string constantTemperature;
            static const std::string densityAtZeroAltitude;
            static const std::string specificGasConstant;
            static const std::string file;
            static const std::string spaceWeatherFile;
        };

        static const std::string ephemeris;
        struct Ephemeris
        {
            static const std::string type;
            static const std::string frameOrigin;
            static const std::string frameOrientation;
            static const std::string makeMultiArc;
            static const std::string correctForStellarAbberation;
            static const std::string correctForLightTimeAbberation;
            static const std::string convergeLighTimeAbberation;
            static const std::string initialTime;
            static const std::string finalTime;
            static const std::string timeStep;
            static const std::string interpolator;
            static const std::string useLongDoubleStates;
            static const std::string bodyIdentifier;
            static const std::string useCircularCoplanarApproximation;
            static const std::string constantState;
            // static const std::string customStateFunction;
            static const std::string initialStateInKeplerianElements;
            static const std::string epochOfInitialState;
            static const std::string centralBodyGravitationalParameter;
            static const std::string rootFinderAbsoluteTolerance;
            static const std::string rootFinderMaximumNumberOfIterations;
            static const std::string bodyStateHistory;
        };

        static const std::string gravityField;
        struct GravityField
        {
            static const std::string type;
            static const std::string gravitationalParameter;
            static const std::string referenceRadius;
            static const std::string cosineCoefficients;
            static const std::string sineCoefficients;
            static const std::string associatedReferenceFrame;
            static const std::string file;
            static const std::string maximumDegree;
            static const std::string maximumOrder;
            static const std::string gravitationalParameterIndex;
            static const std::string referenceRadiusIndex;
        };

        static const std::string rotationModel;
        struct RotationModel
        {
            static const std::string type;
            static const std::string originalFrame;
            static const std::string targetFrame;
            static const std::string initialOrientation;
            static const std::string initialTime;
            static const std::string rotationRate;
        };

        static const std::string shapeModel;
        struct ShapeModel
        {
            static const std::string type;
            static const std::string radius;
            static const std::string equatorialRadius;
            static const std::string flattening;
        };

        static const std::string radiationPressure;
        struct RadiationPressure
        {
            static const std::string type;
            static const std::string referenceArea;
            static const std::string radiationPressureCoefficient;
            static const std::string ocultingBodies;
        };

        static const std::string gravityFieldVariations;
        struct GravityFieldVariation
        {
            static const std::string bodyDeformationType;
            static const std::string modelInterpolation;
            static const std::string deformingBodies;
            static const std::string loveNumbers;
            static const std::string referenceRadius;
            static const std::string cosineCoefficientCorrections;
            static const std::string sineCoefficientCorrections;
            static const std::string minimumDegree;
            static const std::string minimumOrder;
        };
    };


    static const std::string propagation;
    struct Propagator
    {
        static const std::string integratedStateType;
        static const std::string initialStates;

        static const std::string termination;
        struct Termination
        {
            static const std::string conditions;
            static const std::string stopIfSingleConditionMet;
            static const std::string variable;
            static const std::string limitValue;
            static const std::string useAsLowerLimit;
        };

        static const std::string output;
        struct Output
        {
            static const std::string variables;
            static const std::string printVariableTypes;
            static const std::string onlyLastStep;
            static const std::string saveFrequency;
            static const std::string numericalPrecision;
        };

        struct DependentVariable
        {
            static const std::string name;
            static const std::string associatedBody;
            static const std::string secondaryBody;
            static const std::string accelerationType;
            static const std::string bodyUndergoingAcceleration;
            static const std::string bodyExertingAcceleration;
            static const std::string torqueType;
            static const std::string bodyUndergoingTorque;
            static const std::string bodyExertingTorque;
            static const std::string baseFrame;
            static const std::string targetFrame;
            static const std::string angle;
        };

        static const std::string printInterval;

        // Hybrid
        static const std::string propagators;

        // Translational
        static const std::string type;
        static const std::string centralBodies;
        static const std::string bodiesToPropagate;
        // static const std::string initialStateTypes;

        static const std::string accelerations;
        struct Acceleration
        {
            static const std::string type;
            static const std::string maximumDegree;
            static const std::string maximumOrder;
            static const std::string maximumDegreeOfBodyExertingAcceleration;
            static const std::string maximumOrderOfBodyExertingAcceleration;
            static const std::string maximumDegreeOfBodyUndergoingAcceleration;
            static const std::string maximumOrderOfBodyUndergoingAcceleration;
            static const std::string maximumDegreeOfCentralBody;
            static const std::string maximumOrderOfCentralBody;
            static const std::string calculateSchwarzschildCorrection;
            static const std::string calculateLenseThirringCorrection;
            static const std::string calculateDeSitterCorrection;
            static const std::string primaryBody;
            static const std::string centralBodyAngularMomentum;
            static const std::string constantAcceleration;
            static const std::string sineAcceleration;
            static const std::string cosineAcceleration;

            static const std::string direction;
            struct ThrustDirection
            {
                static const std::string type;
                static const std::string relativeBody;
                static const std::string colinearWithVelocity;
                static const std::string towardsRelativeBody;
            };

            static const std::string magnitude;
            struct ThrustMagnitude
            {
                static const std::string type;
                static const std::string originID;
                static const std::string constantMagnitude;
                static const std::string specificImpulse;
                static const std::string bodyFixedDirection;
                static const std::string useAllEngines;
            };

            // static const std::string thrustFrame;
            // static const std::string centralBody;

            /*
            static const std::string fullThrustInterpolationInterface;
            struct FullThrustInterpolationInterface
            {

            };
            */
        };

        static const std::string massRates;
        struct MassRate
        {
            // static const std::string type;
        };

        static const std::string torques;
        struct Torque
        {
            static const std::string type;
        };

    };

    static const std::string integrator;
    struct Integrator
    {
        static const std::string type;
        static const std::string initialTime;
        static const std::string initialTimeStep;
        static const std::string saveFrequency;
        static const std::string rungeKuttaCoefficientSet;
        static const std::string minimumStepSize;
        static const std::string maximumStepSize;
        static const std::string relativeErrorTolerance;
        static const std::string absoluteErrorTolerance;
        static const std::string safetyFactorForNextStepSize;
        static const std::string maximumFactorIncreaseForNextStepSize;
        static const std::string minimumFactorDecreaseForNextStepSize;
    };

    static const std::string interpolator;
    struct Interpolator
    {
        static const std::string type;
        static const std::string lookupScheme;
        static const std::string useLongDoubleTimeStep;
        static const std::string order;
        static const std::string boundaryHandling;
    };

    struct ModelInterpolation
    {
        static const std::string initialTime;
        static const std::string finalTime;
        static const std::string timeStep;
        static const std::string interpolator;
    };

    static const std::string validation;
    struct Validation
    {
        static const std::string usingDefaultValueForMissingKey;
        static const std::string unusedKey;
    };
};

// FIXME: what about arrays?
//! Class for specifying a key path (key.subkey.subsubkey ...) used to access data from `json` objects.
class KeyPath : public std::vector< std::string >
{
public:
    //! Inherit constructors.
    using vector< std::string >::vector;

    //! Constructor with a single string key.
    /*!
     * Constructor with a single string key.
     * \param key The key to be accessed.
     */
    KeyPath( const std::string& key ) : KeyPath( std::vector< std::string >( { key } ) ) { }

    //! Constructor with a single char key.
    /*!
     * Constructor with a single char key.
     * \param key The key to be accessed.
     */
    //! Constructor with a single char key.
    KeyPath( const char* key ) : KeyPath( std::string( key ) ) { }

    //! Constructor with an element index.
    /*!
     * Constructor with an element index.
     * \param key The index of the element to be accessed.
     */
    //! Constructor with a single char key.
    KeyPath( unsigned int vectorIndex ) : KeyPath( std::to_string( vectorIndex ) ) { }

    /*
    //! -DOC
    KeyPath& operator+ ( const KeyPath& subkeys ) const
    {
        KeyPath* compoundKeyPath = new KeyPath( *this );
        for ( std::string subkey : subkeys )
        {
            compoundKeyPath->push_back( subkey );
        }
        return *compoundKeyPath;
    }

    //! -DOC
    KeyPath& operator+ ( const unsigned int vectorIndex ) const
    {
        KeyPath* compoundKeyPath = new KeyPath( *this );
        compoundKeyPath->push_back( std::to_string( vectorIndex ) );
        return *compoundKeyPath;
    }
    */

    //! -DOC
    bool isAbsolute( ) const
    {
        if ( size( ) == 0 )
        {
            return false;
        }
        return *begin( ) == SpecialKeys::root;
    }

    //! -DOC
    KeyPath canonical( const KeyPath& basePath = SpecialKeys::root ) const;
};

//! String representation for `KeyPath`, as key.subkey.vectorIndex.subsubkey ...
inline std::ostream& operator<< ( std::ostream& stringRepresentation, KeyPath const& keyPath )
{
    for ( unsigned int i = 0; i < keyPath.size( ); ++i )
    {
        if ( keyPath.at( i ) != SpecialKeys::root )
        {
            stringRepresentation << keyPath.at( i );
            if ( i < keyPath.size() - 1 )
            {
                stringRepresentation << ".";
            }
        }
    }
    return stringRepresentation;
}

inline KeyPath operator / ( KeyPath path1, const KeyPath& path2 )
{
    for ( std::string subkey : path2 )
    {
        path1.push_back( subkey );
    }
    return path1;
}

inline KeyPath operator / ( const KeyPath& path, const std::string& str )
{
    return path / KeyPath( str );
}

inline KeyPath operator / ( const std::string& str, const KeyPath& path )
{
    return KeyPath( str ) / path;
}

inline KeyPath operator / ( const std::string& str1, const std::string& str2 )
{
    return KeyPath( str1 ) / KeyPath( str2 );
}

inline KeyPath operator / ( const KeyPath& path, const unsigned int vectorIndex )
{
    return path / KeyPath( vectorIndex );
}

inline KeyPath operator / ( const unsigned int vectorIndex, const KeyPath& path )
{
    return KeyPath( vectorIndex ) / path;
}

inline KeyPath operator / ( const std::string& str, const unsigned int vectorIndex )
{
    return KeyPath( str ) / vectorIndex;
}

inline KeyPath operator / ( const unsigned int vectorIndex, const std::string& str )
{
    return vectorIndex / KeyPath( str );
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_KEYS_H
