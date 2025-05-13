#ifndef slic3r_FillPolynomial_hpp_
#define slic3r_FillPolynomial_hpp_

#include "../libslic3r.h"

#include "Fill.hpp"

namespace Slic3r {

class FillPolynomial : public Fill
{
public:

    FillPolynomial(){}
    virtual Fill* clone() const { return new FillPolynomial(*this); };
    virtual ~FillPolynomial() {}

protected:

    virtual void _fill_surface_single(
        unsigned int                     thickness_layers,
        const std::pair<float, Point>   &direction, 
        ExPolygon                       &expolygon, 
        Polylines                       *polylines_out);
};

} // namespace Slic3r

#endif // slic3r_sin_hpp_
