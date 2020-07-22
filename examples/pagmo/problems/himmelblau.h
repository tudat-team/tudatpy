#ifndef TUDAT_EXAMPLE_PAGMO_HIMMELBLAU_H
#define TUDAT_EXAMPLE_PAGMO_HIMMELBLAU_H

#include<math.h>

struct HimmelblauFunction {

    // Empty constructor
    // Without an empty constructor the problem is not accepted
    // as a multithreading type
    HimmelblauFunction( ){ }

    //Actual constructor allowing the user to define the boundaries
    HimmelblauFunction( const double x_min, const double x_max, const double y_min,
            const double y_max ) :
        x_min_( x_min ), x_max_( x_max ), y_min_( y_min ), y_max_( y_max )
    { }

    // Mandatory, computes the fitness, i.e. the Himmelblau's function
    std::vector< double > fitness( const std::vector< double > &x ) const
    {

        std::vector< double > return_value;

        return_value.push_back( pow( x[0]*x[0] + x[1] - 11.0, 2.0 )
                + pow( x[0] + x[1]*x[1] - 7.0, 2.0 ) );

        return return_value;

    }

    // Mandatory, returns the box-bounds
    std::pair< std::vector< double >, std::vector< double > > get_bounds( ) const
    {

        std::pair< std::vector< double >, std::vector< double > > box_bounds;

        box_bounds.first.push_back( x_min_ );
        box_bounds.first.push_back( y_min_ );

        box_bounds.second.push_back( x_max_ );
        box_bounds.second.push_back( y_max_ );

        return box_bounds;

    }

private:

    //Storage members
    double x_min_;
    double x_max_;
    double y_min_;
    double y_max_;

};

#endif
