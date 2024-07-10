// /*    Copyright (c) 2010-2019, Delft University of Technology
//  *    All rigths reserved
//  *
//  *    This file is part of the Tudat. Redistribution and use in source and
//  *    binary forms, with or without modification, are permitted exclusively
//  *    under the terms of the Modified BSD license. You should have received
//  *    a copy of the license with this file. If not, please or visit:
//  *    http://tudat.tudelft.nl/LICENSE.
//  *
//  */

// #include <tudat/astro/aerodynamics/rarefiedFlowAerodynamicCoefficientInterface.h>
// #include <tudat/astro/aerodynamics/rarefiedFlowInteractionModel.h>
// #include <tudat/astro/system_models/vehicleExteriorPanels.h>
// // #include "tudat/astro/ephemerides/rotationalEphemeris.h"



// #include <Eigen/Core>
// #include <iostream>
// #include <iostream>
// #include <fstream>
// #include <vector>


// int main( )
// {
    
//     // print hello world
//     std::cout << "Hello World!" << std::endl;


//     tudat::aerodynamics::RarefiedFlowInteractionModel InteractionModel;    
    
//     double Gj = InteractionModel.get_Gj(1.523);

//     // print Gj
//     std::cout << "Gj: " << Gj << std::endl;


//     std::string mesh_data_path = "/tudat-bundle/build/zz-test/mesh_data";
//     std::vector<double> mesh_panel_normals;
//     std::vector<double> mesh_panel_position_vectors;
//     std::vector<double> mesh_panel_surface_areas;


//     // create vehicle with exterior panels
//     std::ifstream mesh_panel_normals_input(mesh_data_path + "/free_panel_unit_normals3.dat");
//     double dval;
//     while (mesh_panel_normals_input >> dval) {
//         mesh_panel_normals.push_back(dval);
//     }
//     mesh_panel_normals_input.close();

//     std::ifstream mesh_panel_position_vectors_input(mesh_data_path + "/free_panel_position_vectors3.dat");
//     while (mesh_panel_position_vectors_input >> dval) {
//         mesh_panel_position_vectors.push_back(dval);
//     }
//     mesh_panel_position_vectors_input.close();

//     std::ifstream mesh_panel_surface_areas_input(mesh_data_path + "/free_panel_surface_areas3.dat");
//     while (mesh_panel_surface_areas_input >> dval) {
//         mesh_panel_surface_areas.push_back(dval);
//     }
//     mesh_panel_surface_areas_input.close();

//     // for (int i = 0; i < mesh_panel_surface_areas.size(); i++) {
//     //     std::cout << mesh_panel_surface_areas[i] << std::endl;
//     // }


//     // Initialize vehicle class


    





//     // initialize vehicleExteriorPanels
//     std::map< std::string, std::vector< std::shared_ptr< tudat::system_models::VehicleExteriorPanel > > > vehicleExteriorPanels;

//     // initialize single oart of vehicle
//     std::vector< std::shared_ptr< tudat::system_models::VehicleExteriorPanel > > panels;

//     for (int i_panel = 0; i_panel < mesh_panel_surface_areas.size(); i_panel ++) {

//         Eigen::Vector3d normal_vector(mesh_panel_normals[i_panel * 3], mesh_panel_normals[i_panel * 3 + 1], mesh_panel_normals[i_panel * 3 + 2]);

//         Eigen::Vector3d position_vector(mesh_panel_position_vectors[i_panel * 3], mesh_panel_position_vectors[i_panel * 3 + 1], mesh_panel_position_vectors[i_panel * 3 + 2]);

//         // Create a shared pointer using std::make_shared
//         auto panel = std::make_shared<tudat::system_models::VehicleExteriorPanel>(
//             mesh_panel_surface_areas[i_panel],  // Assuming each panel corresponds to three values in the vectors
//             normal_vector,
//             position_vector
//         );

//         // Add the shared pointer to the vector
//         panels.push_back(panel);

//     }
    
//     vehicleExteriorPanels["vehicle"] = panels;

//     for (auto it = vehicleExteriorPanels.begin(); it != vehicleExteriorPanels.end(); ++it) {
//         const std::string& vehiclePartName = it->first;

//         std::cout << "Vehicle part name: " << vehiclePartName << std::endl;
//     }

//     // // initialize vehiclePartOrientation
//     // std::map< std::string, std::shared_ptr< tudat::ephemerides::RotationalEphemeris > > vehiclePartOrientation;
//     // vehiclePartOrientation["vehicle"] = std::make_shared< tudat::ephemerides::RotationalEphemeris >();

//     tudat::aerodynamics::RarefiedFlowAerodynamicCoefficientInterface coeff_interface(
//         vehicleExteriorPanels,
//         1.0, 1.0, Eigen::Vector3d(0.0, 0.0, 0.0)
//     );

//     // print coefficients 

//     std::cout << "Coefficients: " << coeff_interface.getCurrentAerodynamicCoefficients() << std::endl;

//     // upadte coefficients

//     std::vector<double> independentVariables = {0.0, 0.0, 7000.0, 395.0, 2.33420255e+07, 7.69265012e+10, 2.67610829e+11, 2.86047330e+10, 1.18813533e+09, 2.63177678e+06, 2.10337778e+06, 4.06079823e-33};
//     //                                          aoa side    V       T       He             O               N2               O2              Ar              H               N                O+
//     coeff_interface.updateCurrentCoefficients(independentVariables, 0.0);

//     std::cout << "Coefficients: " << coeff_interface.getCurrentAerodynamicCoefficients() << std::endl;

//     // print goodbye
//     std::cout << "Goodbye World!" << std::endl;


//     return 0;

// }