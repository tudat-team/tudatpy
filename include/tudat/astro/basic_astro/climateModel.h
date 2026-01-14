#ifndef TUDAT_CLIMATEMODEL_H
#define TUDAT_CLIMATEMODEL_H

namespace tudat
{

namespace environment
{

class ClimateModelSettings {

    public:

    ClimateModelSettings() = default;

    virtual ~ClimateModelSettings() = default;

};

class ClimateModel {

    public:

    ClimateModel() = default;

    virtual ~ClimateModel() = default;

    virtual void update( 
        double altitude,
        double longitude,
        double latitude,
        double time ) = 0;
        
};

}

}

#endif