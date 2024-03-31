#include <iostream>

#include "../ClipperUtils.hpp" //configuration of clipper (and basic function like offcet)
#include "../PolylineCollection.hpp" //defintion of polyline object like struct Chaining 
#include "../Surface.hpp"  // No idea

#include "Fillsin.hpp"

namespace Slic3r {

void
Fillsin::_fill_surface_single(
    unsigned int                    thickness_layers,
    const direction_t               &direction,
    ExPolygon                       &expolygon,
    Polylines*                      polylines_out)
{
    
    std::cout << "la fonction Fillsin a été appelé" << std::endl;
}


}