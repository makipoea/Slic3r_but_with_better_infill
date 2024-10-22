#ifndef slic3r_Fillpavage_hpp_
#define slic3r_Fillpavage_hpp_


#include <map>
#include<vector>
#include <dlfcn.h> // sur le long termes (demain) pour pouvoir charger une fonction de densité externe 

#include "../libslic3r.h"
#include "Fill.hpp"

#include <fstream>
#include <functional>
#include <cmath>
#include <mutex> // implementer imperativement un témoin
#include "Point.hpp"
//#include"triangle.hpp" // il faudra trouver comment compiler une librairie externe si ca va vite devnir ingérable 
 
using namespace std;

namespace Slic3r {

class Triangle {
    public:
        int depth;
        bool is_feuille;
		bool is_root = false;
        function<int(Point)> density;
        Point sommet[3];
        Point center;
        int max_depth;
        vector<Triangle*> children = {nullptr, nullptr, nullptr, nullptr};
        Triangle* feuilles[4];

        Triangle(int d, Point p0, Point p1, Point p2, function<int(Point)> densFunc = [](Point) { return 2; }, bool is_f = false, int max_dete = 10);
        ~Triangle();

        void display();
        void compute_feuille();
        void update();
};


class Tree {
    public:
        Triangle* root;
        function<int(Point)> density;
        int max_depth=5;

        Tree(function<int(Point)> densFunc = [](Point) { return 2; }, int max_dete = 10);
        ~Tree();

        void saveToFile(const string& filename);
        void createRoot(Point p0, Point p1, Point p2);
		void exportPolylines(Polylines* poly, Triangle* triangle);

    private:
        void saveTriangle(ofstream& outFile, Triangle* triangle);
};


class Fillpavage : public Fill

{
public:
    virtual Fill* clone() const {return new Fillpavage(*this); };
    virtual ~Fillpavage() {}

protected:
	virtual void _fill_surface_single(
	    unsigned int                     thickness_layers,
	    const direction_t               &direction, 
	    ExPolygon                       &expolygon, 
	    Polylines*                      polylines_out);
};
/*
int densityFunction_pavage(Point p){
	return 2; 
}
*/
} //namespace Slic3r



#endif //slic3r_Fillpavage_hpp_