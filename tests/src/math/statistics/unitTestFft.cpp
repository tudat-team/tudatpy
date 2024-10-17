#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/math/statistics/fastFourierTransform.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
    namespace unit_tests
    {

        BOOST_AUTO_TEST_SUITE( test_fft_interface )

            using namespace tudat::fftw_interface;
            BOOST_AUTO_TEST_CASE( testFft )
            {
                std::vector< double > timeDomainData;
                int numberOfDataPoints = 1E7;
                timeDomainData.resize( numberOfDataPoints );
                for( unsigned int i = 0; i < timeDomainData.size(); i++ )
                {
                    timeDomainData[ i ] =  std::sin( 5.0 * 2.0 * mathematical_constants::PI * static_cast< double >( i ) / static_cast< double >( numberOfDataPoints ) ) +
                                           std::cos( 17.0 * 2.0 * mathematical_constants::PI * static_cast< double >( i ) / static_cast< double >( numberOfDataPoints ) ) + 19.0+ static_cast< double >( i );
                }

                std::vector< std::complex< double > > frequencyDomainData = performFftOfRealData( timeDomainData );

                std::vector< double > backTransformedTimeDomainData = performInverseFftToRealData( frequencyDomainData );

                for( unsigned int i = 0; i < frequencyDomainData.size( ); i++ )
                {
                    BOOST_CHECK_SMALL( ( backTransformedTimeDomainData[ i ] - timeDomainData[ i ] ) / timeDomainData[ i ], 2.0E-10 );
                }
            }
        BOOST_AUTO_TEST_SUITE_END( )

    }

}


///* Start reading here */

//#define NUM_POINTS 256


///* Never mind this bit */

//#include <stdio.h>
//#include <math.h"
//#include "Mathematics/Statistics/fastFourierTransform.h"

//#define REAL 0
//#define IMAG 1

//void acquire_from_somewhere(fftw_complex* signal) {
//    /* Generate two sine waves of different frequencies and
//     * amplitudes.
//     */

//    int i;
//    for (i = 0; i < NUM_POINTS; ++i) {
//        double theta = (double)i / (double)NUM_POINTS * M_PI;

//        signal[i][REAL] = 1.0 * cos(32.0 * theta) +
//                          0.5 * cos(16.0 * theta);

//        signal[i][IMAG] = 1.0 * sin(32.0 * theta) +
//                          0.5 * sin(16.0 * theta);
//    }
//}

//void do_something_with(fftw_complex* result) {
//    int i;
//    for (i = 0; i < NUM_POINTS; ++i) {
//        double mag = sqrt(result[i][REAL] * result[i][REAL] +
//                          result[i][IMAG] * result[i][IMAG]);

//        printf("%g\n", mag);
//    }
//}


///* Resume reading here */

//int main() {
//    fftw_complex signal[NUM_POINTS];
//    fftw_complex result[NUM_POINTS];

//    fftw_plan plan = fftw_plan_dft_1d(NUM_POINTS,
//                                      signal,
//                                      result,
//                                      FFTW_FORWARD,
//                                      FFTW_ESTIMATE);

//    acquire_from_somewhere(signal);
//    fftw_execute(plan);

//    fftw_destroy_plan(plan);
//    do_something_with(result);

//    return 0;
//}
