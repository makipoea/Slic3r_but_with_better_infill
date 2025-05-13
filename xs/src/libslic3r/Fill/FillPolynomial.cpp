#include <cstddef>
#include <iostream>
#include <linux/limits.h>
#include <vector>
#include <tuple>
#include <cmath>
#include <mutex>
#include <complex>
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


const int scale_from_mm = 10e6;


namespace Slic3r {

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

    BoundingBox b_box = expolygon.bounding_box();

    double max_norme = 0;
    double min_norme = INFINITY;

    for (int i=0; i<4; i++){
        double norme = std::sqrt(pow(b_box.polygon()[i].x, 2)+pow(b_box.polygon()[i].y, 2));
        if (norme > max_norme){
            max_norme = norme;
        }
        if (norme < min_norme){
            min_norme = norme;
        }
    }

    //std::cout << min_norme << "   " << max_norme << std::endl; 

    Tensor3D polynome_totale = open_polynome_file("/home/makipoea/Documents/prepa/tipe/programme/a-playground-for-oscar/si_avec_parsimonie/analyse_contraintes/polynome.json");
    Tensor2D polynome_densite = slice_polynome_on_z_axis(polynome_totale, this->z);

    //polynome_densite = {{0.9}};

    int n = polynome_densite.size()-1;

    std::vector<double> angle_equi;

    for (int i = 0; i<=n; i++){
        angle_equi.push_back(i*pi/(n+1));
    }

    Tensor2D polynome_phi = decomposition_polynomial(polynome_densite, angle_equi);

    std::vector<double> a_vect;
    std::vector<double> b_vect; 

    for (int i=0; i<=n; i++){
        double a = INFINITY;
        double b = -INFINITY;
        for (int j=0; j<4; j++){
            
            double theta = angle_equi[i];
            double v = std::cos(theta)*b_box.polygon()[j].x + std::sin(theta)*b_box.polygon()[j].y;
            if (v > b){
                b = v;
            }
            if (v < a){
                a = v;
            }
        }
        
        a_vect.push_back(a/scale_from_mm);
        b_vect.push_back(b/scale_from_mm);    
    }

    Tensor2D point_tensor = compute_point_integration(polynome_phi, density , a_vect, b_vect);
    
    for (int i=0; i<=n; i++){
        for (int j=0; j<point_tensor[i].size(); j++){
            double theta = angle_equi[i];
            double x = point_tensor[i][j];
            double height = 200;
            Point p1 = Point(x*scale_from_mm, height*scale_from_mm).rotated(theta);
            Point p2 = Point(x*scale_from_mm, -height*scale_from_mm).rotated(theta);

            path_out.push_back(Line(p1, p2));  
        }
    }
    

    *polylines_out = intersection_pl(path_out, expolygon);

}

} // namespace Slic3r
