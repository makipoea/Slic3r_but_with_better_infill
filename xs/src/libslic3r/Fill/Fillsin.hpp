#ifndef slic3r_Fillsin_hpp_
#define slic3r_Fillsin_hpp_

#include <map>

#include "../libslic3r.h"
#include "Fill.hpp"


namespace Slic3r {

class Fillsin : public Fill

{
public:
    virtual Fill* clone() const {return new Fillsin(*this); };
    virtual ~Fillsin() {}

protected:
	virtual void _fill_surface_single(
	    unsigned int                     thickness_layers,
	    const direction_t               &direction, 
	    ExPolygon                       &expolygon, 
	    Polylines*                      polylines_out);
};

} //namespace Slic3r

#endif //slic3r_Fillsin_hpp_