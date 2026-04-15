# Set the source files.
set(metric_partials_SOURCES
  "metricPartial.cpp"
  "schwarzschildMetricPartial.cpp"
  "solarSystemMetricPartials.cpp"
)

# Set the header files.
set(metric_partials_HEADERS
  "metricPartial.h"
  "christoffelSymbolPartial.h"
  "schwarzschildMetricPartial.h"
  "schwarzschildMetricChristoffelSymbolPartials.h"
  "numericalSolarSystemMetricPartial.h"
  "solarSystemMetricPartials.h"
)

TUDAT_ADD_LIBRARY("metric_partials"
        "${metric_partials_SOURCES}"
        "${metric_partials_HEADERS}"
        )
