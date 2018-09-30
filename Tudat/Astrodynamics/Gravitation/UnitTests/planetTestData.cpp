/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Dekens, E. Orbit Analysis of a Low Frequency Array for Radio Astronomy, MSc thesis, Delft
 *          University of Technlogy, Delft, The Netherlands, 2012.
 *      MathWorks. gravityzonal, MATLAB 2012b, 2012.
 *      Melman, J. Propagate software, J.C.P.Melman@tudelft.nl, 2012.
 *      Ronse, A. A parametric study of space debris impact footprints, MSc thesis, Delft
 *          University of Technlogy, Delft, The Netherlands, in preparation.
 *
 *    Notes
 *      This file is used by the unitTestGravitationalAcceleration.cpp file. The test data in
 *      included in this file was obtained using the gravityzonal() function in MATLAB
 *      (Mathworks, 2012), and code written as part of MSc thesis work by (Dekens, 2012).
 *
 */

#include "Tudat/Astrodynamics/Gravitation/UnitTests/planetTestData.h"

namespace tudat
{
namespace unit_tests
{

using Eigen::Vector3d;

//! Get planet test data generated using MATLAB (Mathworks, 2012).
std::vector< PlanetTestData > getPlanetMatlabTestData( )
{
    return { getMercuryMatlabTestData( ), getVenusMatlabTestData( ),
                getEarthMatlabTestData( ), getMoonMatlabTestData( ), getMarsMatlabTestData( ),
                getJupiterMatlabTestData( ), getSaturnMatlabTestData( ),
                getUranusMatlabTestData( ), getNeptuneMatlabTestData( ) };
}

//! Add (raw) expected accelerations to planet data.
void addExpectedAccelerations( PlanetTestData& planetData, const Eigen::MatrixXd& rawData )
{
    // Add expected accelerations.
    unsigned int rawDataCounter = 0;
    for ( unsigned int k = 0; k < planetData.body1Positions.size( ); k++ )
    {
        for ( unsigned int j = 0; j < planetData.body2Positions.size( ); j++ )
        {
            planetData.expectedAcceleration[ k ][ j ] [ central ] = rawData.row( rawDataCounter );
            rawDataCounter++;

            for ( PlanetTestData::KeyIntValueDoubleMap::iterator coefficientIterator
                  = planetData.zonalCoefficients.begin( );
                  coefficientIterator != planetData.zonalCoefficients.end( );
                  coefficientIterator++ )
            {
                planetData.expectedAcceleration[ k ][ j ] [ coefficientIterator->first ]
                        = rawData.row( rawDataCounter );
                rawDataCounter++;
            }
        }
    }
}

//! Get Mercury test data generated using MATLAB (Mathworks, 2012).
PlanetTestData getMercuryMatlabTestData( )
{
    // Add Mercury data.
    PlanetTestData mercuryData(
                "Mercury",
                2.2032e13,
                2439.0e3,
    { { 2, 0.00006 } },
    { Vector3d( 0.0, 0.0, 0.0 ), Vector3d( 101.1, -253.6, 301.9 ) },
    { Vector3d( 1513.3e3, -7412.67e3, 3012.1e3 ), Vector3d( -5413.3e3, 1292.68e3, -234.6e3 ) } );

    // Set matrix with raw data generated using the gravityzonal() function in MATLAB
    // (Mathworks, 2012).
    Eigen::Matrix< double, 8, 3 > rawData;
    rawData << -6.174552714649318e-02,     3.024510782481964e-01,    -1.228994266291893e-01,
            -6.174568462599339e-02,     3.024518496375884e-01,    -1.229017246366501e-01,
            6.899950362122290e-01,    -1.647687701422098e-01,     2.990280152501966e-02,
            6.900068357055310e-01,    -1.647715878262475e-01,     2.990434476502257e-02,
            -6.174962116268525e-02,     3.024809921551458e-01,    -1.229034674116504e-01,
            -6.174977870471156e-02,     3.024817638759681e-01,    -1.229057657744596e-01,
            6.899448181058168e-01,    -1.647860228710743e-01,     2.993854430409946e-02,
            6.899566157555552e-01,    -1.647888406149411e-01,     2.994008928301272e-02;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( mercuryData, rawData );

    // Return complete benchmark data set for Mercury.
    return mercuryData;
}

//! Get Venus test data generated using MATLAB (Mathworks, 2012).
PlanetTestData getVenusMatlabTestData( )
{
    // Add Venus data.
    PlanetTestData venusData(
                "Venus",
                3.257e14,
                6052.0e3,
    { { 2, 0.000027 } },
    { Vector3d( 0.0, 0.0, 0.0 ), Vector3d( 101.1, -253.6, 301.9 ) },
    { Vector3d( 6485.67e3, -1015.35e3, -4556.31e3 ), Vector3d( -9145.01e3, -3211.55e3, 2032.1e3 ) } );

    // Set matrix with raw data generated using the gravityzonal() function in MATLAB
    // (Mathworks, 2012).
    Eigen::Matrix< double, 8, 3 > rawData;
    rawData << -4.139826709649422e+00,     6.481015916077354e-01,     2.908309216386705e+00,
            -4.139766549333593e+00,     6.480921733399730e-01,     2.908402075753050e+00,
            3.066666265459139e+00,     1.076953665970327e+00,    -6.814396614152983e-01,
            3.066702883991467e+00,     1.076966525688085e+00,    -6.814684119433942e-01,
            -4.139672237863435e+00,     6.479256408539760e-01,     2.908438730972106e+00,
            -4.139612062211696e+00,     6.479162223914803e-01,     2.908531580164727e+00,
            3.066747381498370e+00,     1.076885203277098e+00,    -6.813489125253342e-01,
            3.066784004176275e+00,     1.076898063292949e+00,    -6.813776601420389e-01;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( venusData, rawData );

    // Return complete benchmark data set for Venus.
    return venusData;
}

//! Get Earth test data generated using MATLAB (Mathworks, 2012).
PlanetTestData getEarthMatlabTestData( )
{
    // Add Earth data.
    PlanetTestData earthData(
                "Earth",
                3.986004415e14,
                6378.1363e3,
    { { 2, 0.0010826269 }, { 3, -0.0000025323 }, { 4, -0.0000016204 } },
    { Vector3d( 0.0, 0.0, 0.0 ), Vector3d( 101.1, -253.6, 301.9 ) },
    { Vector3d( 6678.1e3, -4567.12e3, 7601.21e3 ), Vector3d( -8841.1e3, 1234.56e3, -9851.1e3 ) } );

    // Set matrix with raw data generated using the gravityzonal() function in MATLAB
    // (Mathworks, 2012).
    Eigen::Matrix< double, 16, 3 > rawData;
    rawData << -1.945785152656619e+00,     1.330712970216244e+00,    -2.214749937890272e+00,
            -1.944382972538238e+00,     1.329754026076105e+00,    -2.215528487941962e+00,
            -1.944383423600883e+00,     1.329754334555647e+00,    -2.215525868927271e+00,
            -1.944382813097204e+00,     1.329753917035160e+00,    -2.215524898366989e+00,
            1.499913732311049e+00,    -2.094460527945536e-01,     1.671262644735312e+00,
            1.498935088958935e+00,    -2.093093962770631e-01,     1.671421637452590e+00,
            1.498934433363110e+00,    -2.093093047304929e-01,     1.671422558346843e+00,
            1.498934347484291e+00,    -2.093092927384834e-01,     1.671422160014216e+00,
            -1.945951248850231e+00,     1.330772811959008e+00,    -2.214884553154215e+00,
            -1.944548885502504e+00,     1.329813781244983e+00,    -2.215663237316093e+00,
            -1.944549336580869e+00,     1.329814089722815e+00,    -2.215660617849919e+00,
            -1.944548725921873e+00,     1.329813672112975e+00,    -2.215659647116151e+00,
            1.499824437124720e+00,    -2.094742096267836e-01,     1.671195253047338e+00,
            1.498845876819640e+00,    -2.093375381995119e-01,     1.671354207858707e+00,
            1.498845221263043e+00,    -2.093374466406623e-01,     1.671355128623619e+00,
            1.498845135414254e+00,    -2.093374346505208e-01,     1.671354730344184e+00;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( earthData, rawData );

    // Return complete benchmark data set for Earth.
    return earthData;
}

//! Get Moon test data generated using MATLAB (Mathworks, 2012).
PlanetTestData getMoonMatlabTestData( )
{
    // Add Moon data.
    PlanetTestData moonData(
                "Moon",
                4902.799e9,
                1738.0e3,
    { { 2, 0.0002027 } },
    { Vector3d( 0.0, 0.0, 0.0 ), Vector3d( 101.1, -253.6, 301.9 ) },
    { Vector3d( 1843.3e3, -2498.1e3, -1432.1e3 ), Vector3d( -2490.1e3, -1641.6e3, -2009.5e3 ) } );

    // Set matrix with raw data generated using the gravityzonal() function in MATLAB
    // (Mathworks, 2012).
    Eigen::Matrix< double, 8, 3 > rawData;
    rawData << -2.261333728893904e-01,     3.064632880241883e-01,     1.756879527558705e-01,
            -2.261355534987627e-01,     3.064662432567999e-01,     1.757172547885159e-01,
            2.624731121408693e-01,     1.730355651943501e-01,     2.118146736464708e-01,
            2.624626541920230e-01,     1.730286707849584e-01,     2.118363165778046e-01,
            -2.261434565485390e-01,     3.064626497028208e-01,     1.757424642737660e-01,
            -2.261456297927033e-01,     3.064655948164013e-01,     1.757717714296486e-01,
            2.624568485483230e-01,     1.729910903662098e-01,     2.118247690439265e-01,
            2.624463852100570e-01,     1.729841937494688e-01,     2.118464060771237e-01;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( moonData, rawData );

    // Return complete benchmark data set for the Moon.
    return moonData;
}

//! Get Mars test data generated using MATLAB (Mathworks, 2012).
PlanetTestData getMarsMatlabTestData( )
{
    // Add Earth data.
    PlanetTestData marsData(
                "Mars",
                4.305e13,
                3397.2e3,
    { { 2, 0.001964 }, { 3, 0.000036 } },
    { Vector3d( 0.0, 0.0, 0.0 ), Vector3d( 101.1, -253.6, 301.9 ) },
    { Vector3d( -4121.1e3, -4567.12e3, 2654.9e3 ), Vector3d( 7512.1e3, 934.16e3, -7813.6e3 ) } );

    // Set matrix with raw data generated using the gravityzonal() function in MATLAB
    // (Mathworks, 2012).
    Eigen::Matrix< double, 12, 3 > rawData;
    rawData << 5.898668659584946e-01,     6.537072045949770e-01,    -3.800047420429515e-01,
            5.899628855534998e-01,     6.538136162357381e-01,    -3.806422241944675e-01,
            5.899680981754114e-01,     6.538193930113040e-01,    -3.806441314048611e-01,
            -2.511571466344103e-01,    -3.123240639767851e-02,     2.612374011185459e-01,
            -2.510432118553166e-01,    -3.121823814735727e-02,     2.612689821271689e-01,
            -2.510435138099840e-01,    -3.121827569664071e-02,     2.612683517549593e-01,
            5.899421694375344e-01,     6.537383171006586e-01,    -3.800007143788147e-01,
            5.900382635595487e-01,     6.538448028086878e-01,    -3.806382694569552e-01,
            5.900434774491907e-01,     6.538505805266012e-01,    -3.806401759292082e-01,
            -2.511420761315145e-01,    -3.123943102141862e-02,     2.612353345597516e-01,
            -2.510281431200038e-01,    -3.122525895392210e-02,     2.612669053884275e-01,
            -2.510284451320143e-01,    -3.122529652103756e-02,     2.612662750941693e-01;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( marsData, rawData );

    // Return complete benchmark data set for Mars.
    return marsData;
}

//! Get Jupiter test data generated using MATLAB (Mathworks, 2012).
PlanetTestData getJupiterMatlabTestData( )
{
    // Add Jupiter data.
    PlanetTestData jupiterData(
                "Jupiter",
                1.268e17,
                71492.0e3,
    { { 2, 0.01475 }, { 4, -0.00058 } },
    { Vector3d( 0.0, 0.0, 0.0 ), Vector3d( 101.1, -253.6, 301.9 ) },
    { Vector3d( -80123.4e3, 68534.7e3, -75612.5e3 ), Vector3d( 10450.6e3, -100104.5e3, 76091.3e3 ) } );

    // Set matrix with raw data generated using the gravityzonal() function in MATLAB
    // (Mathworks, 2012).
    Eigen::Matrix< double, 12, 3 > rawData;
    rawData << 4.651543414488976e+00,    -3.978764411507470e+00,     4.389664273209171e+00,
            4.629729160110997e+00,    -3.960105275980040e+00,     4.428053668322301e+00,
            4.629107792469540e+00,    -3.959573780250991e+00,     4.427832609529299e+00,
            -6.596983486340383e-01,     6.319137019964030e+00,    -4.803294064973992e+00,
            -6.558632035006745e-01,     6.282400824338628e+00,    -4.843607829959959e+00,
            -6.557660387997968e-01,     6.281470100380290e+00,    -4.843226427552548e+00,
            4.651509238622939e+00,    -3.978744880869845e+00,     4.389644009153466e+00,
            4.629695150539367e+00,    -3.960085842084608e+00,     4.428032894251573e+00,
            4.629073794689641e+00,    -3.959554355145752e+00,     4.427811835972138e+00,
            -6.596981096580464e-01,     6.319179854544868e+00,    -4.803319735152931e+00,
            -6.558629568242166e-01,     6.282443322831597e+00,    -4.843634073021384e+00,
            -6.557657908929294e-01,     6.281512580441177e+00,    -4.843252666565804e+00;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( jupiterData, rawData );

    // Return complete benchmark data set for Jupiter.
    return jupiterData;
}

//! Get Saturn test data generated using MATLAB (Mathworks, 2012).
PlanetTestData getSaturnMatlabTestData( )
{
    // Add Saturn data.
    PlanetTestData saturnData(
                "Saturn",
                3.794e16,
                60268.0e3,
    { { 2, 0.01645 }, { 4, -0.001 } },
    { Vector3d( 0.0, 0.0, 0.0 ), Vector3d( 101.1, -253.6, 301.9 ) },
    { Vector3d( 61823.7e3, 78919.2e3, -75612.5e3 ), Vector3d( -100405.3e3, 61001.2e3, -10030.8e3 ) } );

    // Set matrix with raw data generated using the gravityzonal() function in MATLAB
    // (Mathworks, 2012).
    Eigen::Matrix< double, 12, 3 > rawData;
    rawData << -1.184681053060722e+00,    -1.512269258596779e+00,     1.448905454130921e+00,
            -1.179206635267083e+00,    -1.505281053867206e+00,     1.458681567612995e+00,
            -1.179051590178298e+00,    -1.505083135684198e+00,     1.458580716980696e+00,
            2.323762616303018e+00,    -1.411801051434772e+00,     2.321510722204138e-01,
            2.338200643640199e+00,    -1.420572869189420e+00,     2.365865899134226e-01,
            2.338468227936784e+00,    -1.420735439922169e+00,     2.367301556630723e-01,
            -1.184670868171594e+00,    -1.512263589882015e+00,     1.448901152077055e+00,
            -1.179196482024189e+00,    -1.505275391665871e+00,     1.458677141226794e+00,
            -1.179041439883195e+00,    -1.505077476286146e+00,     1.458576290413381e+00,
            2.323750590783748e+00,    -1.411798193012421e+00,     2.321566241553212e-01,
            2.338188453505635e+00,    -1.420569938389298e+00,     2.365922266073728e-01,
            2.338456032559031e+00,    -1.420732506448698e+00,     2.367357943293829e-01;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( saturnData, rawData );

    // Return complete benchmark data set for Saturn.
    return saturnData;
}

//! Get Uranus test data generated using MATLAB (Mathworks, 2012).
PlanetTestData getUranusMatlabTestData( )
{
    // Add Uranus data.
    PlanetTestData uranusData(
                "Uranus",
                5.794e15,
                25559.0e3,
    { { 2, 0.012 } },
    { Vector3d( 0.0, 0.0, 0.0 ), Vector3d( 101.1, -253.6, 301.9 ) },
    { Vector3d( 95736e3, -94505.5e3, 12012.8e3 ), Vector3d( -26531.8e3, -25888.9e3, -32406.3e3 ) } );

    // Set matrix with raw data generated using the gravityzonal() function in MATLAB
    // (Mathworks, 2012).
    Eigen::Matrix< double, 8, 3 > rawData;
    rawData << -2.251549858651967e-01,     2.222610566211598e-01,    -2.825208713755991e-02,
            -2.252943868321627e-01,     2.223986658599372e-01,    -2.830600336561831e-02,
            1.287819389684475e+00,     1.256613851966410e+00,     1.572960051256680e+00,
            1.280536941285538e+00,     1.249507866757896e+00,     1.579323760413968e+00,
            -2.251561282890238e-01,     2.222618226486847e-01,    -2.825155029985575e-02,
            -2.252955307980970e-01,     2.223994331857989e-01,    -2.830546575701126e-02,
            1.287814893883585e+00,     1.256592367457817e+00,     1.572963219999945e+00,
            1.280532320141713e+00,     1.249486356630524e+00,     1.579326683554127e+00;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( uranusData, rawData );

    // Return complete benchmark data set for Uranus.
    return uranusData;
}

//! Get Neptune test data generated using MATLAB (Mathworks, 2012).
PlanetTestData getNeptuneMatlabTestData( )
{
    // Add Neptune data.
    PlanetTestData neptuneData(
                "Neptune",
                6.809e15,
                24764.0e3,
    { { 2, 0.004 } },
    { Vector3d( 0.0, 0.0, 0.0 ), Vector3d( 101.1, -253.6, 301.9 ) },
    { Vector3d( -34708.9e3, 12709.3e3, -56032.8e3 ), Vector3d( 26531.8e3, 25888.9e3, -25555.6e3 ) } );

    // Set matrix with raw data generated using the gravityzonal() function in MATLAB
    // (Mathworks, 2012).
    Eigen::Matrix< double, 8, 3 > rawData;
    rawData <<  7.813589654806269e-01,    -2.861089086655852e-01,     1.261397815574186e+00,
            7.797740562566650e-01,    -2.855285650995230e-01,     1.260899315132827e+00,
            -1.979182275573878e+00,    -1.931224116498111e+00,     1.906361067159251e+00,
            -1.976988235594378e+00,    -1.929083240959124e+00,     1.911167951328249e+00,
            7.813489390381045e-01,    -2.861101128320262e-01,     1.261384751305303e+00,
            7.797640661670834e-01,    -2.855297726878257e-01,     1.260886260257521e+00,
            -1.979140764528808e+00,    -1.931209887539759e+00,     1.906350867889173e+00,
            -1.976946725170032e+00,    -1.929068983477118e+00,     1.911157603431198e+00;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( neptuneData, rawData );

    // Return complete benchmark data set for Neptune.
    return neptuneData;
}

//! Get Earth test data generated using code from (Melman, 2012).
PlanetTestData getEarthMelmanTestData( )
{
    // Add Earth data.
    PlanetTestData earthData(
                "Earth",
                3.98600441e14,
                6378140.0,
    { { 2, 1.082629989051944e-3 }, { 3, -2.532153068197590e-6 }, { 4, -1.610987610000000e-6 } },
    { Vector3d( 0.0, 0.0, 0.0 ) },
    { Vector3d( -6.54e6, 6.58e6, -6.47e6 ), Vector3d( -6.57e6, -6.53e6, 6.55e6 ) } );

    // Set matrix with raw data generated code from (Dekens, 2012).
    Eigen::Matrix< double, 8, 3 > rawData;
    rawData <<  0.0,                      0.0,                      0.0,
            -5.91805238122957e-04,     5.9542484202585e-04,      1.25534258583695e-03,
            8.29983229092564e-07,    -8.35059579117595e-07,     2.17101021302807e-06,
            -7.33293586126873e-07,     7.37778562188811e-07,    -2.10467844084083e-07,
            0.0,                      0.0,                      0.0,
            -6.1368874587967e-04,     -6.09952436924543e-04,    -1.22366970451254e-03,
            -7.76518489162587e-07,    -7.71790827128112e-07,     2.1675897371097e-06,
            -7.21586841476186e-07,    -7.17193618697031e-07,     2.39786309956306e-07;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( earthData, rawData );

    // Return complete benchmark data set for Earth.
    return earthData;
}

//! Get Earth test data generated using code from (Ronse, 2012).
/*!
 * Returns Earth test data, generated using code written by (Ronse, 2012). This is used as
 * benchmark testing of Tudat code.
 * \return Structure containing Earth test data.
 */
PlanetTestData getEarthRonseTestData( )
{
    // Add Earth data.
    PlanetTestData earthData(
                "Earth",
                3.98600441e14,
                6378140.0,
    { { 2, 1.082629989051944e-3, }, { 3, -2.532153068197590e-6 }, { 4, -1.610987610000000e-6 } },
    { Vector3d( 0.0, 0.0, 0.0 ) },
    { Vector3d( -6.54e6, 6.58e6, -6.47e6 ), Vector3d( -6.57e6, -6.53e6, 6.55e6 ), Vector3d( -6.39e6, -6.40e6, -6.43e6 ) } );

    // Set matrix with raw data generated code from (Dekens, 2012).
    Eigen::Matrix< double, 12, 3 > rawData;
    rawData <<  1.8016172222914324,       -1.8126362878711966,        1.7823338575268453,
            -5.918052381229873e-04,     5.954248420258802e-04,     1.255342585837012e-03,
            8.299832290925637e-07,    -8.350595791175946e-07,     2.171010213028073e-06,
            -7.332935861268730e-07,     7.377785621888111e-07,    -2.104678440840827e-07,
            1.7934666829550308,        1.7825475555093382,       -1.7880071192321845,
            -6.136887458797010e-04,    -6.099524369245734e-04,    -1.223669704512601e-03,
            -7.765184891625865e-07,    -7.717908271281111e-07,     2.167589737109704e-06,
            -7.215868414761851e-07,    -7.171936186970302e-07,     2.397863099563051e-07,
            1.8640423715657017,        1.8669594957778546,        1.8757108684143133,
            -6.788586975941177e-04,    -6.799210742726687e-04,     1.329534713685728e-03,
            8.435002705885402e-07,     8.448203023109010e-07,     2.433609648747570e-06,
            -8.193069668755583e-07,    -8.205891374027501e-07,    -2.887862018589577e-07;

    // Add raw data to list of expected accelerations.
    addExpectedAccelerations( earthData, rawData );

    // Return complete benchmark data set for Earth.
    return earthData;
}

} // namespace unit_tests
} // namespace tudat
