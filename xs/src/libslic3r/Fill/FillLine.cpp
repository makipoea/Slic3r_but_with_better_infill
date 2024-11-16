#include <cstddef>
#include <iostream>
#include <linux/limits.h>
#include <vector>
#include <tuple>
#include <cmath>
#include <mutex>
//#include<dlfcn.h>

#include "../ClipperUtils.hpp" //configuration of clipper (and basic function like offcet)
#include "../PolylineCollection.hpp" //defintion of polyline object like struct Chaining 
#include "../Surface.hpp"  // No idea

//#include "density_function.h"
//#include "triangle.hpp"

#include "FillLine.hpp"

std::mutex console_mutex_ligne;

namespace Slic3r {

void
FillLine::_fill_surface_single(
    unsigned int                    thickness_layers,
    const direction_t               &direction,
    ExPolygon                       &expolygon,
    Polylines*                      polylines_out)
{
    std::lock_guard<std::mutex> lock(console_mutex_ligne);

    std::cout<<"compute line" << std::endl;

    const int nb_ligne = 10;

    if (expolygon.contour.points.empty()) return ; // si le polygone est vide on ne fait rien 

    Polylines path_out;

    BoundingBox b_box = expolygon.bounding_box();

    long hauteur = b_box.polygon()[3].y - b_box.polygon()[0].y;
    long largeur = b_box.polygon()[1].x - b_box.polygon()[0].x;
    //cout << "largeur = " << unscale(largeur) << endl;
    //cout << "hauteur  = " << unscale(hauteur) << endl;
    long x_min = b_box.polygon()[0].x;
    long x_max = b_box.polygon()[1].x;


    for (int i = 1; i <= nb_ligne; i++){
        coord_t y = b_box.polygon()[0].y + (i * hauteur) / (nb_ligne+1);
        std::vector<Pointf> points;

        

        path_out.push_back(Line(Point(x_min, y), Point(x_max, y)));
    }
    *polylines_out = intersection_pl(path_out, expolygon);

}

} // namespace Slic3r