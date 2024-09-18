/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *   References
 *     Blake, W.B. Missile Datcom User's Manual - 1997 Fortran 90 Version, AFRL-VA-WP-TR-1998-3009
 *         Air Force Research Laboratory, 1998.
 *
 */

#ifndef TUDAT_MISSILE_DATCOM_DATA_H
#define TUDAT_MISSILE_DATCOM_DATA_H

#include <string>
#include <vector>

#include <memory>

namespace tudat
{
namespace input_output
{

//! Class which converts and contains the Missile Datcom output aerodynamic data
/*!
 * Class which converts and contains the Missile Datcom output aerodynamic data
 * It is initialized with the data from the MissileDatcomReader. It first converts
 * the vector obtained from the for004.dat file to usable data. Subsequently, coefficients
 * can be retrieved per type.
 */
class MissileDatcomData
{
public:

    //! Enum of static coefficients.
    /*!
     * Enum of static coefficients available via Missile Datcom output. See documentation
     * for details.
     */
    enum StaticCoefficientNames
    {
        cn = 0,
        cm = 1,
        ca = 2,
        cy = 3,
        cln = 4,
        cll = 5,
        cna = 6,
        cma = 7,
        cyb = 8,
        cnb = 9,
        clb = 10
    };

    //! Enum of dynamic coefficients.
    /*!
     * Enum of dynamic coefficients available via Missile Datcom output. See documentation
     * for details.
     */
    enum DynamicCoefficientNames
    {
        cnq = 0,
        cmq = 1,
        caq = 2,
        cyq = 3,
        clnq = 4,
        cllq = 5,
        cnr = 6,
        cmr = 7,
        car = 8,
        cyr = 9,
        clnr = 10,
        cllr = 11,
        cnp = 12,
        cmp = 13,
        cap = 14,
        cyp = 15,
        clnp = 16,
        cllp = 17,
        cnad = 18,
        cmad = 19
    };

    //! Constructor that reads and processes Missile Datcom output.
    /*!
     * Constructor. Reads Missile Datcom file and sets class member variables from file contents.
     * \param fileNameAndPath Name of the file containg missile datcom output, including the
     *         relative path wrt Tudat root path.
     */
    MissileDatcomData( const std::string& fileNameAndPath );

    //! Access the static coefficient database.
    /*!
     * Accesses the static coefficient database.
     * \param machIndex The machnumber index at which the coefficient is accessed.
     * \param angleOfAttackIndex The angle of attack index at which the coefficient is accessed.
     * \param coefficientIndex The index of the coefficient.
     * \return The value of the coefficient.
     */
    double getStaticCoefficient( int machIndex, int angleOfAttackIndex,
                                 StaticCoefficientNames coefficientIndex );

    //! Access the dynamic coefficient database.
    /*!
     * Accesses the dynamic coefficient database.
     * \param machIndex The machnumber index at which the coefficient is accessed.
     * \param angleOfAttackIndex The angle of attack index at which the coefficient is accessed.
     * \param coefficientIndex The index of the coefficient.
     * \return Dynamic coefficient value.
     */
    double getDynamicCoefficient( int machIndex, int angleOfAttackIndex,
                                  DynamicCoefficientNames coefficientIndex );

    //! Retrieve the angle of attacks.
    /*!
     * Retrieve the vector with the different angle of attacks.
     * \return The vector with the different angles of attack.
     */
    std::vector< double > getAngleOfAttacks( ) { return angleOfAttack_; }

    //! Retrieve the Mach numbers.
    /*!
     * Retrieve the vector with the different angle of attacks.
     * \return The vector with the different angles of attack.
     */
    std::vector< double > getMachNumbers( ) {return machNumber_; }

    //! Retrieve the side slip angle.
    /*!
     * Retrieve the vector with the different side slip angles.
     * \return The vector with the different side slip angles.
     */
    double getSideSlipAngle( ) { return sideslipAngle_; }

    //! Retrieve the Reynolds Numbers.
    /*!
     * Retrieve the vector with the different Reynolds Numbers.
     * \return The vector with the different Reynolds Numbers.
     */
    std::vector< double > getReynoldsNumbers( ) { return reynoldsNumbers_; }

    //! Write the database to space-separated files.
    /*!
     * Writes the stored data to a file. Per angle of attack a txt file is generated.
     * \param fileNameBase Base for the output filename.
     * \param basePrecision Number of digits of precision of the base of the output floating-point
     *         numbers in scientific notation.
     * \param exponentWidth Number of digits used to represent the exponent of the output
     *         floating-point numbers in scientific notation.
     */
    void writeAllCoefficientsToFiles(
            const std::string& fileNameBase,
            const int basePrecision = std::numeric_limits< double >::digits10,
            const int exponentWidth = 2 );

    //! Write the force and moment coefficients to tudat aerodynamic coefficients input files.
    /*!
     * \param fileNameBase Base for the output filename.
     * \param basePrecision Number of digits of precision of the base of the output floating-point
     *         numbers in scientific notation.
     * \param exponentWidth Number of digits used to represent the exponent of the output
     *         floating-point numbers in scientific notation.
     */
    void writeForceAndMomentCoefficientsToFiles(
            const std::string& fileNameBase,
            const int basePrecision = std::numeric_limits< double >::digits10,
            const int exponentWidth = 2 );

private:    

    //! Convert the MissileDatcomData.
    /*!
     * Converts the MissileDatcomData to usable data by categorizing the output
     * into various the coefficients and flight conditions.
     * \param datcomData Vector of data read from Missile Datcom output file that is to be 
     *          converted to usable data.
     */
    void convertDatcomData( const std::vector< double >& datcomData );

    //! Sorted data array with the static coefficients.
    /*!
     * The first index is the Mach number, the second index is the angle of attack.
     * Third index are the coefficients.
     * \see StaticCoefficientNames.
     */    
    double staticCoefficients_[ 20 ][ 20 ][ 11 ];

    //! Sorted data array with the dynamic coefficients.
    /*!
     * The first index is the Machnumber, the second index is the angle of attack.
     * the third index are the coefficients
     * \see DynamicCoefficientNames
     */
    double dynamicCoefficients_[ 20 ][ 20 ][ 20 ];

    //! The number of angles of attack angles.
    int numberOfAnglesOfAttack_;

    //! Data vector with the different angles of attack.
    std::vector< double > angleOfAttack_;

    //! The sideslip angle of the datcom data file.
    double sideslipAngle_;

    //! The roll angle of the datcom data file.
    double rollAngle_;

    //! The number of mach numbers.
    int numberOfMachNumbers_;

    //! Datavector containing the different mach numbers.
    std::vector< double > machNumber_;

    //! Datavector containing the different altitudes used.
    std::vector< double > altitudes_;

    //! Datavector containing the different Reynold numbers.
    std::vector< double > reynoldsNumbers_;

    //! Datavector containing the different free-steam velocities.
    std::vector< double > freeStreamVelocities_;

    //! Datavector containing the different free-steam temperatures.
    std::vector< double > freeStreamStaticTemperatures_;

    //! Datavector containing the different free-steam pressures.
    std::vector< double > freeStreamStaticPressure_;
};

//! Typedef for shared-pointer to MissileDatcomData object.
typedef std::shared_ptr< MissileDatcomData > MissileDatcomDataPointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_MISSILE_DATCOM_DATA_H
