#include <iostream>
#include <cstring>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

extern "C" {
#include <cspice/SpiceUsr.h>
#include <cspice/SpiceZfc.h>
}

int main() {
/* From https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcov_c.html */
#define  FILSIZ         256
#define  LNSIZE         81
#define  MAXCOV         100000
#define  WINSIZ         ( 2 * MAXCOV )
#define  TIMLEN         51

    SPICEDOUBLE_CELL        (cover, WINSIZ);

    SpiceBoolean found;

    SpiceChar file[FILSIZ];
    SpiceChar idch[LNSIZE];
    SpiceChar meta[FILSIZ];
    SpiceChar source[FILSIZ];
    SpiceChar timstr[TIMLEN];
    SpiceChar type[LNSIZE];

    SpiceDouble b;
    SpiceDouble e;

    SpiceInt count;
    SpiceInt handle;
    SpiceInt i;
    SpiceInt idcode;
    SpiceInt niv;


//    furnsh_c("/home/dominik/.conda/pkgs/tudat-resources-1.1.2-hc8dc577_6/resource/spice_kernels/tudat_merged_spk_kernel.bsp");
    furnsh_c("/home/dominik/.conda/pkgs/tudat-resources-1.1.2-hc8dc577_6/resource/spice_kernels/naif0012.tls");
//    furnsh_c("/home/dominik/dev/tudat-bundle/spice/lro/data/spk/lrorg_2009169_2010001_v01.bsp");
//    furnsh_c("/home/dominik/dev/tudat-bundle/spice/lro/data/spk/lrorg_2010001_2010091_v01.bsp");

    std::string path = "/home/dominik/dev/tudat-bundle/spice/lro/data/spk";
    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path), {})) {
        if (entry.path().extension() == ".bsp")
        {
            furnsh_c(entry.path().string().c_str());
        }
    }

    idcode = -85;

    ktotal_c("SPK", &count);

    for (i = 0; i < count; i++) {
        kdata_c(i, "SPK", FILSIZ, LNSIZE, FILSIZ,
                file, type, source, &handle, &found);

        spkcov_c(file, idcode, &cover);
    }

    niv = wncard_c(&cover);

    printf("\nCoverage for object %d\n", (int) idcode);

    for (i = 0; i < niv; i++) {
        wnfetd_c(&cover, i, &b, &e);

        timout_c(b,
                 "YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
                 TIMLEN,
                 timstr);

        printf("\n"
               "Interval:  %d\n"
               "Start:     %s\n",
               (int) i,
               timstr);

        timout_c(e,
                 "YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
                 TIMLEN,
                 timstr);
        printf("Stop:      %s\n", timstr);

    }
    return (0);
}