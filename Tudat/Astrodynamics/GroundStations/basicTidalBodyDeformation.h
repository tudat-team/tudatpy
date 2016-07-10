#ifndef BASICTIDALBODYDEFORMATION_H
#define BASICTIDALBODYDEFORMATION_H

#include <vector>
#include <map>
#include <assert.h>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/GroundStations/groundStation.h"
#include "Tudat/Astrodynamics/GroundStations/bodyDeformationModel.h"

namespace tudat
{

namespace site_displacements
{

//! Calculate displacement at given site due to degree 2 tide of single external body, band-independent Love and Shida numbers.
/*!
 *  Calculate displacement at given site due to degree 2 tide of single external body, band-independent Love and Shida numbers.
 *  This calculation is part of 'Step 1' of solid earth tide site displacement calcualtion,
 *  IERS 2010 Conventions, Section 7.1.1.
 *  \param gravitationalParameterRatio. ratio of gravtational parameters of body causing deformation and body being deformed.
 *  \param stationPositionUnitVector Unit vector of position vetor of site at which displacement is to be calculated, in
 *  body-fixed, body-centered frame.
 *  \param relativeBodyState Position of body causing deformation in frame fixed to and centered at body being deformed.
 *  \param bodyEquatorialRadius Equatorial radius of body being deformed.
 *  \param degreeTwoLoveNumber Degree two displacement Love number h_{2}.
 *  \param degreeTwoShidaNumber Degree two displacement Shida number l_{2}.
 *  \return Displacement at given site due to degree two tide from body at given relative state.
 */
Eigen::Vector3d calculateDegreeTwoInPhaseDisplacement( const double gravitationalParameterRatio,
                                                       const Eigen::Vector3d& stationPositionUnitVector,
                                                       const Eigen::Vector3d& relativeBodyState,
                                                       const double bodyEquatorialRadius,
                                                       const double degreeTwoLoveNumber,
                                                       const double degreeTwoshidaNumber );

//! Calculate displacement due to degree 3 tide, band-independet love and shida numbers.
//! Calculate displacement at given site due to degree 3 tide of single external body, band-independent Love and Shida numbers.
/*!
 *  Calculate displacement at given site due to degree 3 tide of single external body, band-independent Love and Shida numbers.
 *  This calculation is part of 'Step 1' of solid earth tide site displacement calcualtion,
 *  IERS 2010 Conventions, Section 7.1.1.
 *  \param gravitationalParameterRatio. ratio of gravtational parameters of body causing deformation and body being deformed.
 *  \param stationPositionUnitVector Unit vector of position vetor of site at which displacement is to be calculated, in
 *  body-fixed, body-centered frame.
 *  \param relativeBodyState Position of body causing deformation in frame fixed to and centered at body being deformed.
 *  \param bodyEquatorialRadius Equatorial radius of body being deformed.
 *  \param degreeThreeLoveNumber Degree three displacement Love number h_{3}.
 *  \param degreeThreeShidaNumber Degree three displacement Shida number l_{3}.
 *  \return Displacement at given site due to degree three tide from body at given relative state.
 */
Eigen::Vector3d calculateDegreeThreeInPhaseDisplacement( const double gravitationalParameterRatio,
                                                         const Eigen::Vector3d& stationPositionUnitVector,
                                                         const Eigen::Vector3d& relativeBodyState,
                                                         const double bodyEquatorialRadius,
                                                         const double degreeThreeLoveNumber,
                                                         const double degreeThreeShidaNumber );


class BasicTidalBodyDeformation: public BodyDeformationModel
{
public:
    BasicTidalBodyDeformation(
            const boost::function< basic_mathematics::Vector6d( const double ) > deformedBodyEphemeris,
            const std::vector< boost::function< basic_mathematics::Vector6d( const double ) > >& deformingBodyEphemerides,
            const boost::function< Eigen::Quaterniond( const double ) > deformedBodyRotation,
            const std::vector< int >& maximumTidalDegrees,
            const boost::function< double( ) > gravitionalParameterOfDeformedBody,
            const std::vector< boost::function< double( ) > >& gravitionalParametersOfDeformingBodies,
            const boost::function< double( ) > deformedBodyEquatorialRadius,
            const std::vector< double >& nominalDisplacementLoveNumbers,
            const std::vector< double >& nominalDisplacementShidaNumbers );

    virtual ~BasicTidalBodyDeformation( ){ }

    virtual Eigen::Vector3d calculateSiteDisplacement( const double ephemerisTime,
                                                       const Eigen::Vector3d& nominalSiteUnitVector );

    virtual Eigen::Vector3d calculateSiteDisplacement(
            const double time,
            const boost::shared_ptr< ground_stations::NominalGroundStationState > siteState )
    {
        return calculateSiteDisplacement( time, siteState->getNominalCartesianPosition( ).normalized( ) );
    }


    double getDisplacementLoveNumber( const int index )
    {
        return nominalDisplacementLoveNumbers_[ index ];
    }

    void setDisplacementLoveNumber( const int index, const double value )
    {
        nominalDisplacementLoveNumbers_[ index ] = value;
    }

    double getDisplacementShidaNumber( const int index )
    {
        return nominalDisplacementShidaNumbers_[ index ];
    }

    void setDisplacementShidaNumber( const int index, const double value )
    {
        nominalDisplacementShidaNumbers_[ index ] = value;
    }

    int getNumberOfDeformingBodies( )
    {
        return deformingBodyEphemerides_.size( );
    }

    boost::function< basic_mathematics::Vector6d( const double ) > getDeformingBodyStateFunction( int bodyIndex )
    {
        return deformingBodyEphemerides_[ bodyIndex ];
    }

    boost::function< Eigen::Quaterniond( const double ) > getDeformedBodyRotationFunction( )
    {
        return deformedBodyRotation_;
    }

    int getMaximumTidalDegree( const int bodyIndex )
    {
        return maximumTidalDegrees_[ bodyIndex ];
    }

    boost::function< double( ) > getDeformedBodyGravitationalParameterFunction( )
    {
        return gravitionalParameterOfDeformedBody_;
    }

    boost::function< double( ) > getDeformingBodyGravitationalParameterFunction( int bodyIndex )
    {
        return gravitionalParametersOfDeformingBodies_[ bodyIndex ];
    }

    boost::function< double( ) > getDeformedBodyEquatorialRadiusFunction( )
    {
        return deformedBodyEquatorialRadius_;
    }

    double getDisplacementLoveNumbers( const int index )
    {
        return nominalDisplacementLoveNumbers_[ index ];
    }

    double getDisplacementShidaNumbers( const int index )
    {
        return nominalDisplacementShidaNumbers_[ index ];
    }

protected:
    boost::function< basic_mathematics::Vector6d( const double ) > deformedBodyEphemeris_;

    std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > deformingBodyEphemerides_;

    boost::function< Eigen::Quaterniond( const double ) > deformedBodyRotation_;

    int numberOfBodies_;

    std::vector< int > maximumTidalDegrees_;

    boost::function< double( ) > gravitionalParameterOfDeformedBody_;

    std::vector< boost::function< double( ) > > gravitionalParametersOfDeformingBodies_;

    boost::function< double( ) > deformedBodyEquatorialRadius_;

    std::vector< double > nominalDisplacementLoveNumbers_;

    std::vector< double > nominalDisplacementShidaNumbers_;

};


}

}

#endif // BASICTIDALBODYDEFORMATION_H
