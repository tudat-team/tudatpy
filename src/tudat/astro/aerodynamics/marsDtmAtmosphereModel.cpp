
#include <iostream>
#include <fstream>
#include <vector>
#include "tudat/astro/aerodynamics/marsDtmAtmosphereModel.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/io/basicInputOutput.h"
#include <tuple>

namespace tudat
{
namespace aerodynamics
{

    void MarsDtmAtmosphereModel::computeProperties(
            const double altitude, const double longitude,
            const double latitude, const double time )
    {
        // Compute the hash key??
        basic_astrodynamics::DateTime currentDateTime_ = basic_astrodynamics::DateTime::fromTime( time );
        currentF107_ =  f107Function_( time );
        if (currentF107_ > 100.0) {
            currentF107_ = 100.0;
        }
        //std::cout << "F107: " << currentF107_ << std::endl;
        currentDensity_ = getTotalDensity(
                altitude, latitude, longitude,
                currentDateTime_.getMinute(), currentDateTime_.getHour(), currentDateTime_.getDay(), currentDateTime_.getMonth(), currentDateTime_.getYear() );
        //std::cout << "current Density: " << currentDensity_ << std::endl;
        currentGeopotentialAltitude_ = computeGeopotentialAltitude( altitude );
        currentTemperature_138 = computeCurrentTemperature(
                latitude, longitude, currentDateTime_.getMinute(), currentDateTime_.getHour(), currentDateTime_.getDay(), currentDateTime_.getMonth(), currentDateTime_.getYear(), 3 );
        currentTemperature_inf = computeCurrentTemperature(
                latitude, longitude, currentDateTime_.getMinute(), currentDateTime_.getHour(), currentDateTime_.getDay(), currentDateTime_.getMonth(), currentDateTime_.getYear(), 1 );
        currentdTemperature_ = computeCurrentTemperature(
                latitude, longitude, currentDateTime_.getMinute(), currentDateTime_.getHour(), currentDateTime_.getDay(), currentDateTime_.getMonth(), currentDateTime_.getYear(), 5 );
        sigma = currentdTemperature_/(currentTemperature_inf-currentTemperature_138);
        currentTemperature_z = currentTemperature_inf - (currentTemperature_inf - currentTemperature_138)* exp(-sigma*currentGeopotentialAltitude_);
    }
    // function to load the coefficients from the file
    std::vector<std::vector<double>> loadCoefficients(const std::string& filename){
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }
        // define the number of rows and columns in the file
        const int rows = 70;
        //const int rows = 74;
        const int cols = 25;
        std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (!(file >> matrix[i][j])) {
                    throw std::runtime_error(
                    "Error reading value at row " + std::to_string(i) +
                    ", column " + std::to_string(j)
                );
                }
            }
        }
        // Check for extra data
        double extra;
        if (file >> extra) {
            throw std::runtime_error("File " + filename + " contains more data than expected (more than 70×25 values).");
        }
        file.close();
        return matrix;
    }

    // function to provide the nearest date in a list of dates to a target date (i
    marsDate tudat::aerodynamics::marsDate::findNearestDate(const marsDate& targetDate){
        if (targetDate.year < marsYears[0].year){
            std::cerr << "The target date is before the first date of the Martian years: 1955, 4, 11. Day of year can't be calculated. PROGRAM TERMINATED!!" << std::endl;
            exit(EXIT_FAILURE);
        }
            marsDate nearestDate = marsYears[0];
            double minDifference = marsDate::dateDifference(targetDate, nearestDate);
            for (const marsDate& date : marsYears){
                if (date.year < targetDate.year || (date.year == targetDate.year && date.month < targetDate.month) || (date.year == targetDate.year && date.month == targetDate.month && date.day < targetDate.day)){
                    double difference = marsDate::dateDifference(targetDate, date);
                    if (difference < minDifference){
                        minDifference = difference;
                        nearestDate = date;
                    }
                }

            }
            return nearestDate;

    }

    // function to compute the day of martian year
    double tudat::aerodynamics::marsDate::marsDayofYear(const marsDate& targetDate){
        marsDate nearestDate = findNearestDate(targetDate);
        double difference = marsDate::dateDifference(targetDate, nearestDate);
        // convert the difference in days to seconds, and devide them by the number of seconds in a martian year to get the day of year
        double dayOfYear = difference*24*3600 / 88668.0; // 88668 seconds is the length of a martian day (24.63 h * 3600 s/h)
        return dayOfYear;
    }

    std::tuple<double, double> computeSolarLongitude(const double longitude, const int day ,const int month,const int year, const int hours, const int minutes){

        //Day since J2000 epoch in Terrestrial Time (or TAI)
        double timeElaspsedJD = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                year, month, day, hours, minutes, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000);
        //std::cout << "timeElaspsedJD: " << timeElaspsedJD << std::endl;
        //Mars orbital parameters:
        //Mars Mean Anomaly
        double meanAnomaly = 19.3870 + 0.52402073 * timeElaspsedJD; //degrees
        //std::cout << "MeanAnomaly: " << MeanAnomaly << std::endl;
        //Angle of Fiction Mean Sun
        double omegaFMS = 270.3863 + 0.52403840 * timeElaspsedJD;
        //std::cout << "OmegaFMS: " << OmegaFMS << std::endl;
        //Perturbers
        double PBS =
                0.0071 * cos(unit_conversions::convertDegreesToRadians(0.985626 * timeElaspsedJD / 2.2353 + 49.409))
                + 0.0057 * cos(unit_conversions::convertDegreesToRadians(0.985626 * timeElaspsedJD / 2.7543 + 168.173))
                + 0.0039 * cos(unit_conversions::convertDegreesToRadians(0.985626 * timeElaspsedJD / 1.1177 + 191.837))
                + 0.0037 * cos(unit_conversions::convertDegreesToRadians(0.985626 * timeElaspsedJD / 15.7866 + 21.736))
                + 0.0021 * cos(unit_conversions::convertDegreesToRadians(0.985626 * timeElaspsedJD / 2.1354 + 15.704))
                + 0.0020 * cos(unit_conversions::convertDegreesToRadians(0.985626 * timeElaspsedJD / 2.4694 + 95.528))
                + 0.0018 * cos(unit_conversions::convertDegreesToRadians(0.985626 * timeElaspsedJD / 32.8493 + 49.095));
        //std::cout << "PBS: " << PBS << std::endl;
        //Equation of Center
        double EoC = (10.691 + 3.0E-7 * timeElaspsedJD) * sin(unit_conversions::convertDegreesToRadians(meanAnomaly))
                     + 0.623 * sin(unit_conversions::convertDegreesToRadians(2 * meanAnomaly))
                     + 0.050 * sin(unit_conversions::convertDegreesToRadians(3 * meanAnomaly))
                     + 0.005 * sin(unit_conversions::convertDegreesToRadians(4 * meanAnomaly))
                     + 0.0005 * sin(unit_conversions::convertDegreesToRadians(5 * meanAnomaly)) +
                     PBS;//Equation of Center
        //std::cout << "EoC: " << EoC << std::endl;
        //Areocentric Solar Longitude
        double Ls1 = omegaFMS + EoC;
        double Ls = basic_mathematics::computeModulo(Ls1, 360.0);
        return std::make_tuple(Ls, EoC);
    }



    MarsDtmAtmosphereModel::MarsDtmAtmosphereModel(const std::function< double( const double ) > f107Function ) :
    polarRadius_( 3378.0e3 ), // polar radius of Mars
    filename_ ( paths::getAtmosphereTablesPath( ) + "/dtm_mars.dat" ), // file name of the coefficients
    f107Function_( f107Function ),
    alpha_( {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.38, -0.40, 0.0}), // thermal diffusion coefficients
    mmass_( {44.01, 16.00, 28.0, 40.0, 28.0, 32.0, 4.0, 1.0, 2.0} )// molar mass of the species
    {
        //std::cout<<"MarsDtmAtmosphereModel constructor"<<std::endl;
        coefficients_ = loadCoefficients( filename_ ); // load the coefficients from the file
        currentLegendrePolynomials_.resize( 7 );
        currentLegendrePolynomials_[ 0 ] = 0.0;
        Ls_ = 0.0;
        //std::cout << "Loaded Matrix values:" << std::endl;
        //std::cout << coefficients_[0][3] << "\t";
    }

    // function to compute the local true solar time
    // equations from: https://www.giss.nasa.gov/tools/mars24/help/algorithm.html
    double MarsDtmAtmosphereModel::computeLocalSolarTime(const double longitude, const int day ,const int month,const int year, const int hours, const int minutes)
    {

        auto outResults = computeSolarLongitude(longitude, day, month, year, hours, minutes);
        double EoC = std::get<1>(outResults);
        Ls_ = std::get<0>(outResults);
        //std::cout << "Ls: " << Ls_ << std::endl;
        //Equation of Time
        double EoT = 2.861 * sin( unit_conversions::convertDegreesToRadians(2*Ls_))
                     - 0.071 * sin( unit_conversions::convertDegreesToRadians(4*Ls_))
                     + 0.002 * sin( unit_conversions::convertDegreesToRadians(6*Ls_)) - EoC; // Equation of Time
        //std::cout << "EoT: " << EoT << std::endl;
        //Mean Solar Time at Mars's prime merdian, i.e., Airy Mean Time
        double MST =  24*((basic_astrodynamics::convertCalendarDateToJulianDay(year,month,day,hours,minutes,0.0) - 2451549.5)/1.0274912517 + 44796.0 - 0.0009626); //double MST
        //std::cout << "MST1: " << MST << std::endl;
        MST = basic_mathematics::computeModulo( MST,24.0);// Modulo MST
        //std::cout << "MST: " << MST << std::endl;
        //Local Mean Solar Time
        double LMST =  basic_mathematics::computeModulo( MST-longitude*(24.0/360.0),24.0);
        //Local True Solar Time
        double LTST = basic_mathematics::computeModulo(LMST + EoT*(24.0/360.0),24.0);
        //std::cout << "LTST: " << LTST << std::endl;
        //Subsolar Longitude
        double subSolarLongitude = MST*(360.0/24.0) + EoT + 180.0;
        // return LTST in seconds.
        return LTST*3600.0;
    }
    // Function to define the dust storm parameters
    void MarsDtmAtmosphereModel::defineDustPars(const int rows )
    {
        double tau = 0.0; // dust opacity
        double lsdeb = 0.0; // dust storm onset
        double lsfin = 0.0; // dust storm end
        double taubg = 0.0; // background dust opacity

        if (rows >= 74){
            tau = coefficients_[1][rows-3];
            lsdeb = coefficients_[1][rows-2];
            lsfin = coefficients_[1][rows-1];
            taubg = coefficients_[1][rows];

            if (taubg<0.15){
                taubg = 0.0;
            }

        }

        taus = {tau, lsdeb, lsfin, taubg};


    }

    // function to compute the dust storm
    double MarsDtmAtmosphereModel::computeDustStorm(const double Ls){
        MarsDtmAtmosphereModel::defineDustPars(70);

        double dells = 15.0; //delta Ls of Stewart
        double xlsmaxd = taus[1] + 6.0; //tau max reached after 6 degree, and xlsdev is storm onset
        double dectau = 12.0; // decay of tau = max to normal in 12 degree of Ls
        double xlsmaxf = taus[2] - dectau; // start os storm decay (tau)
        double deuxse = 2.0/2.78128;
        double correc = pow(deuxse,dectau);
        double amr = 0.0;
        double da41 = 0.0;
        //Dust storm
        if (Ls >= taus[1] && Ls >= 180){
            if (Ls < xlsmaxd){
                double q = exp(-(Ls - taus[1])/dells);
                double da41 = q*(1.0-pow(q,4))*taus[0];
            }
            if (Ls >= xlsmaxd && Ls < xlsmaxf){
                double q = exp(-(xlsmaxd - Ls)/dells);
                double da41 = q*(1.0-pow(q,4))*taus[0];
            }
            if (Ls >= xlsmaxf && Ls < taus[2]){
                double qmax = exp(-(xlsmaxd - Ls)/dells);
                double taumax = qmax*(1.0-pow(qmax,4))*taus[0];
                double xmult = taumax/deuxse;
                double ipow = (Ls - xlsmaxf) + 1.0;
                double denom = pow(deuxse,ipow);
                double da41 = xmult*(denom-correc);
            }

        }
        return da41;

    }

    void MarsDtmAtmosphereModel::updateLegendrePolynomials( const double latitude )
    {
        for( int i = 1; i <= 6; i++ )
        {
            currentLegendrePolynomials_[ i ] = basic_mathematics::computeLegendrePolynomialExplicit( i, 0, sin( latitude ) );
        }
    }

    // Function to compute the spherical harmonic expansions defined as G(l) as in the paper (Bruinsma and Lemoine, 2002)
    double MarsDtmAtmosphereModel::computeGl(const double latitude, const double longitude, const int minutes_, const int hours_, const int day_ ,const int month_,const int year_, const int indexg)
    {
        updateLegendrePolynomials( latitude );
        using basic_mathematics::computeLegendrePolynomialExplicit;
        marsDate date2 = marsDate(year_, month_, day_, hours_ , minutes_, 0.0);
        double doy = date2.marsDayofYear(date2); //day of year
        //std::cout << "doy: " << doy << std::endl;
        double t = computeLocalSolarTime(longitude, day_, month_, year_, hours_, minutes_); //seconds
        if (currentF107_ > 100.0) {
            currentF107_ = 100.0;
        }
        //std::cout<<"F107: "<<currentF107_<<std::endl;
        double F = currentF107_ - 65.0;
        // Non-periodic term
        double NP = (coefficients_[1][indexg] * F) + coefficients_[30][indexg] * currentLegendrePolynomials_[1]
                + coefficients_[31][indexg] * currentLegendrePolynomials_[ 2 ]
                + coefficients_[32][indexg] * currentLegendrePolynomials_[ 3 ]
                + coefficients_[33][indexg] * currentLegendrePolynomials_[ 4 ]
                + coefficients_[34][indexg] * currentLegendrePolynomials_[ 5 ];
        // Annual-Periodic terms
        double PA = (coefficients_[2][indexg]+ coefficients_[4][indexg] * currentLegendrePolynomials_[2] + coefficients_[6][indexg]*F)* cos(Omega * doy)
              + (coefficients_[3][indexg]+ coefficients_[5][indexg] * currentLegendrePolynomials_[2] + coefficients_[7][indexg]*F)* sin(Omega * doy)
             + (coefficients_[8][indexg] * currentLegendrePolynomials_[1] + coefficients_[10][indexg] * currentLegendrePolynomials_[3] + coefficients_[12][indexg] *currentLegendrePolynomials_[5] ) * cos(Omega * doy)
             + (coefficients_[9][indexg] * currentLegendrePolynomials_[1] + coefficients_[11][indexg] * currentLegendrePolynomials_[3] + coefficients_[13][indexg] * currentLegendrePolynomials_[5]) * sin(Omega * doy);
        //Semi-Annual-Periodic terms
        double PSA = (coefficients_[14][indexg] +  coefficients_[16][indexg] * currentLegendrePolynomials_[2]  + coefficients_[18][indexg]*F) * cos(2.0*Omega * doy)
              + (coefficients_[15][indexg]+ coefficients_[17][indexg] * currentLegendrePolynomials_[2]  + coefficients_[19][indexg]*F)* sin(2.0*Omega * doy)
              + (coefficients_[20][indexg] * currentLegendrePolynomials_[1] + coefficients_[22][indexg] *  currentLegendrePolynomials_[3]) * cos(2.0*Omega * doy)
              + (coefficients_[21][indexg] * currentLegendrePolynomials_[1] + coefficients_[23][indexg] *  currentLegendrePolynomials_[3]) * sin(2.0*Omega * doy);
        //Diurnal terms
        double PD = (coefficients_[24][indexg]*computeLegendrePolynomialExplicit(1,1,sin(latitude)) + coefficients_[26][indexg]*computeLegendrePolynomialExplicit(2,1,sin(latitude)) + coefficients_[28][indexg]*computeLegendrePolynomialExplicit(3,1,sin(latitude)))*cos(omega*t)
                  + (coefficients_[25][indexg]*computeLegendrePolynomialExplicit(1,1,sin(latitude)) + coefficients_[27][indexg]*computeLegendrePolynomialExplicit(2,1,sin(latitude)) + coefficients_[29][indexg]*computeLegendrePolynomialExplicit(3,1,sin(latitude)))*sin(omega*t);
        double Gl = NP + PA + PSA + PD;
        return Gl;
    }

    //Function to compute the spherical harmonic expansions defined as G(l) as in the subroutine provided by Bruinsma
    double MarsDtmAtmosphereModel::computeGl_Subr(const double latitude, const double longitude, const int minutes_, const int hours_, const int day_ ,const int month_,const int year_, const int indexg)
    {
        //Get the zonal Legendre polynomials for the current latitude
        updateLegendrePolynomials( latitude );
        using basic_mathematics::computeLegendrePolynomialExplicit;
        //Get the dust storm parameters
        MarsDtmAtmosphereModel::defineDustPars(70);
        //Get the day of year
        marsDate date2 = marsDate(year_, month_, day_, hours_ , minutes_, 0.0);
        double doy = date2.marsDayofYear(date2);
        //std::cout << "doy: " << doy << std::endl;
        //Get the local solar time
        double t = computeLocalSolarTime(longitude, day_, month_, year_, hours_, minutes_); //hours
        //Get the local solar time in radians + pi.
        double hl0 = omega*t+ mathematical_constants::PI;
        //std::cout << "hl: " << omega*t << std::endl;
        double cos2h = cos(hl0)* cos(hl0) - sin(hl0)*sin(hl0);
        double sin2h = 2.0*sin(hl0)*cos(hl0);
        //flux terms:
        double ff0 = 0.0;
        if (currentF107_ > 100.0) {
            currentF107_ = 100.0;
        }
        double F = currentF107_ - 65.0;
        //std::cout<<"currentF107 "<<currentF107_<<std::endl;
        double F2 = F*F;
        //std::cout << "F: " << F << std::endl;
        double f0 = coefficients_[4][indexg]*F + coefficients_[39][indexg] * F2;
        double f1f = 1.0 + f0*ff0; // coupling terms
        //flux + latitude terms
        f0 = f0+ coefficients_[1][indexg]*currentLegendrePolynomials_[2] + coefficients_[2][indexg]*currentLegendrePolynomials_[3]
                + coefficients_[3][indexg]*currentLegendrePolynomials_[4] + coefficients_[58][indexg]*currentLegendrePolynomials_[6]
                + coefficients_[59][indexg]*currentLegendrePolynomials_[1] + coefficients_[60][indexg]*currentLegendrePolynomials_[5];
        // symmetrical and seasonal annual terms
        double PA = coefficients_[5][indexg]* cos(Omega*doy) + coefficients_[6][indexg]* sin(Omega*doy)
                + coefficients_[7][indexg]* cos(Omega*doy) * F + coefficients_[8][indexg]* sin(Omega*doy) * F
                + coefficients_[9][indexg]* cos(Omega*doy) * currentLegendrePolynomials_[2] + coefficients_[10][indexg]* sin(Omega*doy) * currentLegendrePolynomials_[2]
                + coefficients_[11][indexg] * cos(Omega*doy) * currentLegendrePolynomials_[1] + coefficients_[12][indexg] * sin(Omega*doy) * currentLegendrePolynomials_[1]
                + coefficients_[17][indexg] * cos(Omega*doy) * currentLegendrePolynomials_[3] + coefficients_[18][indexg] * sin(Omega*doy) * currentLegendrePolynomials_[3]
                + coefficients_[61][indexg] * cos(Omega*doy) * currentLegendrePolynomials_[1] *F + coefficients_[62][indexg] * sin(Omega*doy) * currentLegendrePolynomials_[1] *F
                + coefficients_[63][indexg] * cos(Omega*doy) * currentLegendrePolynomials_[5] + coefficients_[64][indexg] * sin(Omega*doy) * currentLegendrePolynomials_[5];
        //symmetrical and seasonal semi-annual terms
        double PSA = coefficients_[13][indexg] * cos(2.0*Omega*doy) + coefficients_[14][indexg] * sin(2.0*Omega*doy)
                + coefficients_[15][indexg] * cos(2.0*Omega*doy) * F + coefficients_[16][indexg] * sin(2.0*Omega*doy) * F
                + coefficients_[19][indexg] * cos(2.0*Omega*doy) * currentLegendrePolynomials_[1] + coefficients_[20][indexg] * sin(2.0*Omega*doy) * currentLegendrePolynomials_[1]
                + coefficients_[65][indexg] * cos(2.0*Omega*doy) * currentLegendrePolynomials_[2] + coefficients_[66][indexg] * sin(2.0*Omega*doy) * currentLegendrePolynomials_[2]
                + coefficients_[67][indexg] * cos(2.0*Omega*doy) * currentLegendrePolynomials_[3] + coefficients_[68][indexg] * sin(2.0*Omega*doy) * currentLegendrePolynomials_[3];
        // diurnal terms
        double PD = coefficients_[21][indexg] * computeLegendrePolynomialExplicit(1,1,sin(latitude)) * cos(hl0) + coefficients_[22][indexg] * computeLegendrePolynomialExplicit(1,1,sin(latitude)) * sin(hl0)
                + coefficients_[23][indexg] * computeLegendrePolynomialExplicit(2,1,sin(latitude)) * cos(hl0) + coefficients_[24][indexg] * computeLegendrePolynomialExplicit(2,1,sin(latitude)) * sin(hl0)
                + coefficients_[25][indexg] * computeLegendrePolynomialExplicit(3,1,sin(latitude)) * cos(hl0) + coefficients_[26][indexg] * computeLegendrePolynomialExplicit(3,1,sin(latitude)) * sin(hl0)
                + coefficients_[27][indexg] * computeLegendrePolynomialExplicit(1,1,sin(latitude)) * cos(hl0) * F + coefficients_[28][indexg] * computeLegendrePolynomialExplicit(1,1,sin(latitude)) * sin(hl0) * F
                + coefficients_[30][indexg] * computeLegendrePolynomialExplicit(2,2,sin(latitude)) * (cos2h) + coefficients_[31][indexg] * computeLegendrePolynomialExplicit(2,2,sin(latitude)) * sin2h
                + coefficients_[32][indexg] * computeLegendrePolynomialExplicit(3,2,sin(latitude)) * (cos2h) + coefficients_[33][indexg] * computeLegendrePolynomialExplicit(3,2,sin(latitude)) * sin2h
                + coefficients_[34][indexg] * computeLegendrePolynomialExplicit(4,2,sin(latitude)) * (cos2h) + coefficients_[35][indexg] * computeLegendrePolynomialExplicit(4,2,sin(latitude)) * sin2h
                + coefficients_[36][indexg] * computeLegendrePolynomialExplicit(2,2,sin(latitude)) * (cos2h) * F + coefficients_[37][indexg] * computeLegendrePolynomialExplicit(2,2,sin(latitude)) * sin2h * F;
        //dust terms
        double da41 = computeDustStorm( Ls_ );
        double fpds = coefficients_[40][indexg]*da41 + coefficients_[69][indexg]*taus[3];
        double Gl = fpds + f0 + PA + PSA + PD;
        return Gl;
    }

    // Function to compute the geopotential altitude
    double MarsDtmAtmosphereModel::computeGeopotentialAltitude( const double altitude )
    {
        return ((altitude - 138.0E3)*(polarRadius_ + 138.0E3 )/(polarRadius_+altitude))/1000.0; //km //removed the 80 from Bruinsma's formula, because in Bates, 1959, it is not there, results seem to be more realistic withoout the 80
    }

    // Function to compute the current temperature (exospheric, 138 km, and partial temperatures)
    double MarsDtmAtmosphereModel::computeCurrentTemperature( const double latitude, const double longitude, const int minutes_, const int hours_, const int day_ ,const int month_,const int year_, const int indexg)
    {
        double T0 = coefficients_[0][indexg];
        double Ti;
        double Gl = computeGl_Subr(latitude, longitude, minutes_, hours_, day_ , month_, year_, indexg);
        //std::cout << "Gl: " << Gl << std::endl;
        //std::cout << "T0: " << T0 << std::endl;
        Ti = T0*(1.0 + Gl);
        return Ti;
    }
    // Function to compute gamma parameter
    double MarsDtmAtmosphereModel::computeGamma( const double latitude, const double longitude, const int minutes_, const int hours_, const int day_ ,const int month_,const int year_, const int indexm)
    {
        double universalGasConstant = 8.314;//physical_constants::MOLAR_GAS_CONSTANT; //J/mol/K //kg m^2 / s^2 / K / mol
        //double g138 =celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/((138.0E3+3376.78E3)*(138.0E3+3376.78E3)); //m/s2
        //Get gravity at 138 km
        double g138  = 3.727/pow((1+(138.0/3376.78)),2); //m/s2
        currentTemperature_138 = computeCurrentTemperature( latitude, longitude, minutes_, hours_, day_ , month_, year_, 3); //K
        currentTemperature_inf = computeCurrentTemperature( latitude, longitude, minutes_, hours_, day_ , month_, year_, 1); //K
        currentdTemperature_ = computeCurrentTemperature( latitude, longitude, minutes_, hours_, day_ , month_, year_, 5); //K
        double tmt = currentTemperature_inf - currentTemperature_138;
        if (tmt <= 0.0) {
            currentTemperature_138 = 0.99 * currentTemperature_inf;
        }
        sigma = currentdTemperature_/(currentTemperature_inf-currentTemperature_138);
        //std::cout << "sigma: " << sigma << std::endl;
        double gamma = mmass_[indexm] * g138 / ( universalGasConstant * sigma * currentTemperature_inf); //km^-1
        /*
        std::cout << "gamma: " << gamma << std::endl;

        std::cout << "currentTemperature_138: " << currentTemperature_138 << std::endl;
        std::cout << "currentTemperature_inf: " << currentTemperature_inf << std::endl;
        std::cout << "currentdTemperature_: " << currentdTemperature_ << std::endl;
         */
        return gamma;
    }

    // Function to compute the height distribution function
    double MarsDtmAtmosphereModel::heightDistributionFunction(const double altitude, const double latitude,
                                                              const double longitude, const int minutes_, const int hours_, const int day_ ,const int month_,const int year_,const int indexm)
    {
        currentGeopotentialAltitude_ = computeGeopotentialAltitude( altitude ); //km
        //std::cout << "currentGeopotentialAltitude_: " << currentGeopotentialAltitude_ << std::endl;
        double gamma= computeGamma( latitude, longitude, minutes_, hours_, day_ , month_, year_, indexm); //km^-1
        currentTemperature_z = currentTemperature_inf - (currentTemperature_inf - currentTemperature_138)* exp(-sigma*currentGeopotentialAltitude_);
       // std::cout << "T138: " << currentTemperature_138 << std::endl;
       // std::cout << "Tinf: " << currentTemperature_inf << std::endl;
       // std::cout << "exp(-sigma*currentGeopotentialAltitude_): " << exp(-sigma*currentGeopotentialAltitude_) << std::endl;
       // std::cout << "Tz: " << currentTemperature_z << std::endl;

        double fi = pow((currentTemperature_138/currentTemperature_z),(1+alpha_[indexm]+gamma))* exp(-sigma*gamma*currentGeopotentialAltitude_);
        //std::cout << "fi: " << fi << std::endl;
        return fi;
    }

    // Function to compute the total density
    double MarsDtmAtmosphereModel::getTotalDensity( const double altitude, const double latitude,
                                                    const double longitude, const int minutes_, const int hours_, const int day_ ,const int month_,const int year_)
    {
        //std::cout<<"getting total density"<<std::endl;
        //std::cout<<"F107: "<<currentF107_<<std::endl;
        double avogadroConstant_ = 6.022E23;//physical_constants::AVOGADRO_CONSTANT;
        double rho0;
        double rho = 0;
        int indexm = 0;
        for (int col = 7; col < 25; col += 2) {
            rho0 = coefficients_[0][col]*mmass_[indexm]/avogadroConstant_; //g/cm3
            double fi = heightDistributionFunction(altitude, latitude, longitude, minutes_, hours_, day_ , month_, year_, indexm);
            double gl = computeGl_Subr(latitude, longitude, minutes_, hours_, day_ , month_, year_, col);

            rho += rho0*fi*exp(gl);
            indexm +=1;
        }
        return rho*1000;
    }


} // namespace aerodynamics

} // namespace tudat
