#include <cstddef>
#include <iostream>
#include <linux/limits.h>
#include <vector>
#include <tuple>
#include <cmath>
#include <mutex>
#include <complex>
#include <string>
#include <utility>
//#include<dlfcn.h>

#include "../ClipperUtils.hpp" //configuration of clipper (and basic function like offcet)
#include "../PolylineCollection.hpp" //defintion of polyline object like struct Chaining 
#include "../Surface.hpp"  // No idea
//#include "density_function.h"
//#include "triangle.hpp"

#include "FillPolynomial.hpp"
#include "polynome/polynome.hpp"

std::mutex console_mutex_polynomial;

constexpr double pi = 3.14159265358979323846;

std::string json_path = "/home/makipoea/Documents/prepa/tipe/programme/a-playground-for-oscar/si_avec_parsimonie/analyse_contraintes/polynome.json";

const int scale_from_mm = 1e6;


namespace Slic3r {

std::vector<double> sim_to_slicer(double p_x, double p_y, double center_sim_x, double center_sim_y, double center_slicr_x, double center_slicr_y){
    return {(p_x-center_sim_x)*scale_from_mm+center_slicr_x, (p_y-center_sim_y)*scale_from_mm+center_slicr_y};
}

std::vector<double> slicer_to_sim(double p_x, double p_y, double center_sim_x, double center_sim_y, double center_slicr_x, double center_slicr_y){
    return {((p_x-center_slicr_x)/scale_from_mm+center_sim_x), ((p_y-center_slicr_y)/scale_from_mm+center_sim_y)};
}

void
FillPolynomial::_fill_surface_single(
    unsigned int                    thickness_layers,
    const direction_t               &direction,
    ExPolygon                       &expolygon,
    Polylines*                      polylines_out)
{
    std::lock_guard<std::mutex> lock(console_mutex_polynomial);

    //std::cout<<"compute polynome"<< this->z << "/ " <<density <<std::endl;
    Polylines path_out;

    BoundingBox b_box_slicr = expolygon.bounding_box();

    Tensor2D b_box_sim = open_bbox_file(json_path);

    double center_sim_x = (b_box_sim[0][0]+b_box_sim[2][0])/2;
    double center_sim_y = (b_box_sim[0][1]+b_box_sim[2][1])/2;

    Point center_slicr = b_box_slicr.polygon().centroid();

    
    //std::cout << "ratio_x " << b_box_slicr.polygon()[1].x - b_box_slicr.polygon()[0].x << "  " << b_box_sim[1][0]-b_box_sim[0][0] << std::endl;

    //std::cout << "ratio_y " << b_box_slicr.max.y - b_box_slicr.min.y << "  " << b_box_sim[3][1]-b_box_sim[0][1] << std::endl;


    //std::cout << "center_sim = " << center_sim_x << "  " << center_sim_y << std::endl;

    //std::cout << "center_slicr = " << center_slicr.x << "  " << center_slicr.y << std::endl;
    
    //for (int i=0; i<4; i++){
    //   std::cout << b_box_sim[i][0] << "   " << b_box_sim[i][1] << std::endl; 
       // std::cout<< "slicer " <<b_box_slicr.polygon()[i].x << "   " << b_box_slicr.polygon()[i].y  << std::endl;
       // std::cout<< "sim " <<b_box_sim[i][0] << "   " << b_box_sim[i][1]  << std::endl;
       // std::vector<double> ps = slicer_to_sim(b_box_slicr.polygon()[i].x, b_box_slicr.polygon()[i].y, center_sim_x, center_sim_y, center_slicr.x, center_slicr.y);
       // std::cout<< "slicer -> sim " << ps[0] << "   " << ps[1]  << std::endl;
       // std::vector<double> ps2 = sim_to_slicer(ps[0], ps[1], center_sim_x, center_sim_y, center_slicr.x,  center_slicr.y);
       // std::cout<< "slicer -> sim -> slicr " << ps2[0] << "   " << ps2[1]  << std::endl;
    //}
    
    Tensor3D polynome_totale = open_polynome_file(json_path);
    Tensor2D polynome_densite = slice_polynome_on_z_axis(polynome_totale, this->z);
    //polynome_densite = {{0.1}};
    /*
    polynome_densite = {{0.1, 0, 0, 0},
                        {0  , 0, 0, 0}, 
                        {0  , 0, 0, 0}, 
                        {0  , 0, 0, 0}};
    */
                       
    //polynome_densite = {{0.1}};
    
    double eta_x = 0.9;

    Tensor2D polynome_densite_x = multiply_by_eta(polynome_densite, eta_x);

    Tensor2D polynome_densite_y = multiply_by_eta(polynome_densite, 1-eta_x);

    int n = polynome_densite.size()-1;

    int nb_subdivision = 500;
    double largeur_x = b_box_sim[1][0] - b_box_sim[0][0];
    double largeur_y = b_box_sim[3][1] - b_box_sim[0][1];
    std::vector<double> l_x = {};
    std::vector<double> l_y = {};

    Tensor2D poly_phi_x = {};
    Tensor2D poly_phi_y = {};
    std::vector<double> a_x(nb_subdivision, b_box_sim[0][1]);
    std::vector<double> b_x(nb_subdivision, b_box_sim[3][1]);

    std::vector<double> a_y(nb_subdivision, b_box_sim[0][0]);
    std::vector<double> b_y(nb_subdivision, b_box_sim[1][0]);


    for (int i=0; i<nb_subdivision; i++){
        double x = b_box_sim[0][0]+(largeur_x)*i/(nb_subdivision-1);
        l_x.push_back(x);
        poly_phi_x.push_back(evaluate_on_x_axis(polynome_densite_x,  x));

        double y = b_box_sim[0][1]+(largeur_y)*i/(nb_subdivision-1);
        l_y.push_back(y);
        poly_phi_y.push_back(evaluate_on_y_axis(polynome_densite_y,  y));
    }
    //print_vector(a_y);
    //print_vector(b_y);

    std::cout << "xxxxxxx" << std::endl;
    Tensor2D l_points_x = compute_point_integration(poly_phi_x, density, a_x, b_x);
    std::cout << "yyyyyyy" << std::endl;
    Tensor2D l_points_y = compute_point_integration(poly_phi_y, density, a_y, b_y);

    //print_matrix(poly_phi_y);
    //print_vector(a_y);
    //print_vector(b_y);
    //print_matrix(l_points_y);

    //print_vector(l_y);

    auto l_polyline_x = extract_blocks_with_xy(l_points_x, l_x, Axis_maki::X);
    auto l_polyline_y = extract_blocks_with_xy(l_points_y, l_y, Axis_maki::Y);


    //std::cout << " xxxxxxxxxxxxxxxxxx" <<  std::endl;


    //print_l_point(l_polyline_x, 2, " | ");

    //std::cout << " y y y y y yy y" << std::endl;

    //print_l_point(l_polyline_y, 2, " | ");
        
    for (const auto& curve : l_polyline_x) {
        Polyline poly;
        for (const auto& [x_sim, y_sim] : curve) {
            std::vector<double> p_slicer = sim_to_slicer(x_sim, y_sim, center_sim_x, center_sim_y, center_slicr.x, center_slicr.y);
            poly.points.emplace_back(p_slicer[0], p_slicer[1]);
        }
        if (poly.points.size() >= 2) {
            path_out.emplace_back(std::move(poly));
        }
    }

    // Remplissage selon Y
    for (const auto& curve : l_polyline_y) {
        Polyline poly;
        for (const auto& [x_sim, y_sim] : curve) {
            std::vector<double> p_slicer = sim_to_slicer(x_sim, y_sim, center_sim_x, center_sim_y, center_slicr.x, center_slicr.y);
            poly.points.emplace_back(p_slicer[0], p_slicer[1]);
        }
        if (poly.points.size() >= 2) {
            path_out.emplace_back(std::move(poly));
        }
    }


    *polylines_out = intersection_pl(path_out, expolygon);

}

} // namespace Slic3r
