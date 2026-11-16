#ifndef TUDAT_CONSTANT_MODEL_SCALING_H
#define TUDAT_CONSTANT_MODEL_SCALING_H

namespace tudat
{
namespace simulation_setup
{

//! A numerical scaling value with no callback or references to another model.
//! Retaining this type inside std::function lets worker copying recognize known constants
//! without evaluating an arbitrary function or relying on a user safety assertion.
template< typename Value >
class ConstantModelScaling final
{
public:
    explicit ConstantModelScaling( const Value& value ): value_( value ) {}

    Value operator( )( const double ) const
    {
        return value_;
    }

private:
    const Value value_;
};

}  // namespace simulation_setup
}  // namespace tudat

#endif
