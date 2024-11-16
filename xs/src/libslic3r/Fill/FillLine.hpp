#ifndef slic3r_FillLine_hpp_
#define slic3r_FillLine_hpp_

#include "../libslic3r.h"

#include "Fill.hpp"

namespace Slic3r {

class FillLine : public Fill
{
public:

    FillLine(){}
    virtual Fill* clone() const { return new FillLine(*this); };
    virtual ~FillLine() {}

protected:

    virtual void _fill_surface_single(
        unsigned int                     thickness_layers,
        const std::pair<float, Point>   &direction, 
        ExPolygon                       &expolygon, 
        Polylines                       *polylines_out);
};

} // namespace Slic3r

#endif // slic3r_sin_hpp_