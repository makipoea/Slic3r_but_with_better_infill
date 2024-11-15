#include "density_function.h"
#include <cmath>
#include<math.h>
#include <iostream>

int density(long x, long y) {
    std::cout << "x: " << x << ", y: " << y << std::endl;
    //x = (x-50)/50;
    //y = (y-50)/50;
    //std::cout << "density = " << static_cast<float>(5*exp(-x*x-y*y)) << std::endl;
    return 5;//static_cast<int>(5*exp(-x*x-y*y)); // parce que  
}