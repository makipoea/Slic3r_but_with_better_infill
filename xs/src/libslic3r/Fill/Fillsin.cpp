#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <mutex>

#include "../ClipperUtils.hpp" //configuration of clipper (and basic function like offcet)
#include "../PolylineCollection.hpp" //defintion of polyline object like struct Chaining 
#include "../Surface.hpp"  // No idea

#include "Fillsin.hpp"

const double PI = 3.14159265358979323846;

std::mutex console_mutex_sin;

// sinusoidal

namespace Slic3r {

void
Fillsin::_fill_surface_single(
    unsigned int                    thickness_layers,
    const direction_t               &direction,
    ExPolygon                       &expolygon,
    Polylines*                      polylines_out)
{
    std::lock_guard<std::mutex> lock(console_mutex_sin);

    std::cout << "---------------------------------------------" << std::endl;

    if (expolygon.contour.points.empty()) return ; // si le polygone est vide on ne fait rien 

    size_t min_x_index_point = 0;
    size_t max_x_index_point = 0;

    coord_t min_x = expolygon.contour[0].x;
    coord_t max_x = expolygon.contour[0].x;

    //std::cout << "nombre de point " << expolygon.contour.points.size() << std::endl;

    for (unsigned int idx_point = 0; idx_point < expolygon.contour.points.size(); ++idx_point){
        //std::cout <<  "x = " << expolygon.contour[idx_point].x << " y = " << expolygon.contour[idx_point].y << std::endl;
        Point p_tmp = expolygon.contour[idx_point];
        
        if (p_tmp.x < min_x){ 
            min_x_index_point = idx_point; 
            min_x = p_tmp.x;
        } 

        if (p_tmp.x > max_x){
            max_x_index_point = idx_point; 
            max_x = p_tmp.x;
        }
    }

    //std::cout << "max_x = " << max_x << "minx = " << min_x << std::endl;
    //std::cout << "max_x_index point " << max_x_index_point << "min_x index point " << min_x_index_point << std::endl;

    float angle_of_rotation = std::atan(static_cast<double>(expolygon.contour[max_x_index_point].y-expolygon.contour[min_x_index_point].y)/(max_x-min_x));// on tourne le polygone puis

    std::cout<< "angle of rotation = " << angle_of_rotation <<std::endl;

    expolygon.rotate(-angle_of_rotation, expolygon.contour[min_x_index_point]);

    //std::cout<< " x = "<<expolygon.contour[0].x << " / " << expolygon.contour[1].x << " / " << expolygon.contour[2].x << " / " << expolygon.contour[3].x << " / " << expolygon.contour[4].x << "/" << expolygon.contour[5].x << std::endl;

    //std::cout<< "y = " <<expolygon.contour[0].y << " / " << expolygon.contour[1].y << " / " << expolygon.contour[2].y << " / " << expolygon.contour[3].y << " / " << expolygon.contour[4].y << " / " <<  expolygon.contour[5].y << std::endl;

    Polyline path_out;

    size_t nb_point = expolygon.contour.points.size();
    size_t indice_haut = min_x_index_point;
    size_t indice_bas = min_x_index_point;
    coord_t min_x_rotate = expolygon.contour[min_x_index_point].x;
    coord_t max_x_rotate = expolygon.contour[max_x_index_point].x;

    unsigned int nb_period = static_cast<unsigned int>(100.0*density);
    unsigned int nb_subdivision = 1000;

    std::cout<<" min_x_ rotate = " << min_x_rotate << " max_x_rotate = " << max_x_rotate<<std::endl;


/*
    for (size_t i = 1; i <= nb_subdivision; i++) {
        coord_t x = min_x_rotate + i * (max_x_rotate - min_x_rotate) / nb_subdivision;

        while (expolygon.contour[(indice_haut - 1 + nb_point) % nb_point].x < x) {
            indice_haut = (indice_haut - 1 + nb_point) % nb_point;
        }

        while (expolygon.contour[(indice_bas + 1) % nb_point].x < x) {
            indice_bas = (indice_bas + 1) % nb_point;
        }

        float sin_x = std::sin(2 * PI * nb_period / (max_x_rotate - min_x_rotate) * x);
        
        //  / ! \ look at coord_t // mm factor 

        //float sin_x = std::sin(PI/(max_x_rotate-min_x_rotate)+ density*(PI/thickness_layers-PI/(max_x_rotate-min_x_rotate)));

        float t1 = static_cast<double>(x - expolygon.contour[indice_haut].x) / (expolygon.contour[(indice_haut - 1 + nb_point) % nb_point].x - expolygon.contour[indice_haut].x);
        float f1_x = expolygon.contour[indice_haut].y + t1 * (expolygon.contour[(indice_haut - 1 + nb_point) % nb_point].y - expolygon.contour[indice_haut].y);

        //float t1 = (expolygon.contour[(indice_haut - 1 + nb_point) % nb_point].x - x)/(expolygon.contour[(indice_haut - 1 + nb_point) % nb_point].x - expolygon.contour[indice_haut].x);
        //float f1_x = expolygon.contour[indice_haut].y*t1 + (1-t1)*expolygon.contour[(indice_haut - 1 + nb_point) % nb_point].y;

        float t2 = static_cast<double>(x - expolygon.contour[indice_bas].x) / (expolygon.contour[(indice_bas + 1) % nb_point].x - expolygon.contour[indice_bas].x);
        float f2_x = expolygon.contour[indice_bas].y + t2 * (expolygon.contour[(indice_bas + 1) % nb_point].y - expolygon.contour[indice_bas].y);

        //std::cout << "x = " << static_cast<coord_t>(x) << " y = " << static_cast<coord_t>(std::floor((f1_x + f2_x) / 2 + (f1_x - f2_x) / 2 * sin_x)) << std::endl;
        path_out.append(Point(static_cast<coord_t>(x), static_cast<coord_t>((f1_x + f2_x) / 2 + (f1_x - f2_x) / 2 * sin_x)));
        //std::cout << "x = " << x << " indice _haut = " << indice_haut << " indice bas" << indice_bas << " t1 = " << t1 << " f1_x = "<< f1_x << std::endl;
        //std::cout << "x = " << x << " precedent " << expolygon.contour[indice_haut].x << "suivant = " << expolygon.contour[(indice_haut - 1 + nb_point) % nb_point].x << std::endl; 
        //path_out.append(Point(static_cast<coord_t>(x), static_cast<coord_t>(f1_x)));
    }

    path_out.rotate(angle_of_rotation, expolygon.contour[min_x_index_point]);
    *polylines_out=path_out;
*/

    for (size_t i = 1; i <= nb_subdivision; i++) {
    coord_t x = min_x_rotate + i * (max_x_rotate - min_x_rotate) / nb_subdivision;

    while (expolygon.contour[(indice_haut - 1 + nb_point) % nb_point].x < x) {
        indice_haut = (indice_haut - 1 + nb_point) % nb_point;
    }

    while (expolygon.contour[(indice_bas + 1) % nb_point].x < x) {
        indice_bas = (indice_bas + 1) % nb_point;
    }

    // Calculate frequency based on x
    float normalized_x = (x - min_x_rotate) / (max_x_rotate - min_x_rotate); // Normalize x to [0, 1]
    float frequency = 2.0f + 8.0f * std::abs(0.5f - normalized_x); // Adjust frequency (higher in the middle, lower at edges)

    float sin_x = std::sin(frequency * PI * nb_period * normalized_x); // Use frequency to control step size

    float t1 = static_cast<double>(x - expolygon.contour[indice_haut].x) / 
                (expolygon.contour[(indice_haut - 1 + nb_point) % nb_point].x - expolygon.contour[indice_haut].x);
    float f1_x = expolygon.contour[indice_haut].y + t1 * 
                  (expolygon.contour[(indice_haut - 1 + nb_point) % nb_point].y - expolygon.contour[indice_haut].y);

    float t2 = static_cast<double>(x - expolygon.contour[indice_bas].x) / 
                (expolygon.contour[(indice_bas + 1) % nb_point].x - expolygon.contour[indice_bas].x);
    float f2_x = expolygon.contour[indice_bas].y + t2 * 
                  (expolygon.contour[(indice_bas + 1) % nb_point].y - expolygon.contour[indice_bas].y);

    path_out.append(Point(static_cast<coord_t>(x), static_cast<coord_t>((f1_x + f2_x) / 2 + (f1_x - f2_x) / 2 * sin_x)));
    }

    path_out.rotate(angle_of_rotation, expolygon.contour[min_x_index_point]);
    *polylines_out = path_out;
    }

}
