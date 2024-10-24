#include <cstddef>
#include <iostream>
#include <linux/limits.h>
#include <vector>
#include <tuple>
#include <cmath>
#include <mutex>
#include<dlfcn.h>

#include "../ClipperUtils.hpp" //configuration of clipper (and basic function like offcet)
#include "../PolylineCollection.hpp" //defintion of polyline object like struct Chaining 
#include "../Surface.hpp"  // No idea
//#include "triangle.hpp"

#include "Fillpavage.hpp"

std::mutex console_mutex_pavage;

const double PI = 3.14159265358979323846;

namespace Slic3r {
    //cout.setf(std::ios::unitbuf);

Triangle::Triangle(int d, Point p0, Point p1, Point p2, function<int(Point)> densFunc, bool is_f, int max_dete)
    : depth(d), is_feuille(is_f), sommet{p0, p1, p2}, density(densFunc), max_depth(max_dete) {
    //children.resize(4, nullptr);
    try {
    center = Point(static_cast<coord_t>((p0.x+p1.x+p2.x)/3), static_cast<coord_t>((p0.y+p1.y+p2.y)/3));
    }
    catch (...){
        cout << "le centre a échouer" << endl;
    }
    try {
        compute_feuille();
    }
    catch (...){
        cout << "compute a echouer" <<endl;
    }
}


Triangle::~Triangle() {
    cout << "destruction du triangle " << endl;
    for (int i = 0; i < 4; ++i) {
        if (children[i] != nullptr) {
            delete children[i];
            children[i] = nullptr;
        }
        //if (feuilles[i] != nullptr) {
        //    delete feuilles[i];
        //    feuilles[i] = nullptr;
        //}
    }
}

void Triangle::display() {
    for (int i = 0; i < 3; ++i) {
        cout << "(" << sommet[i].x << " , " << sommet[i].y << ")  ";
    }
}

void Triangle::compute_feuille() {
    if (!is_feuille) {
        Point P01 = Point(static_cast<coord_t>((sommet[0].x+ sommet[1].x)/2), static_cast<coord_t>((sommet[0].y+ sommet[1].y)/2)); //(sommet[0] + sommet[1]) / 2;
        Point P12 = Point(static_cast<coord_t>((sommet[1].x+ sommet[2].x)/2), static_cast<coord_t>((sommet[1].y+ sommet[2].y)/2)); //(sommet[1] + sommet[2]) / 2;
        Point P20 = Point(static_cast<coord_t>((sommet[2].x+ sommet[0].x)/2), static_cast<coord_t>((sommet[2].y+ sommet[0].y)/2)); //(sommet[2] + sommet[0]) / 2;
        feuilles[0] = new Triangle(depth + 1, sommet[0], P01, P20, density, true, max_depth);
        feuilles[1] = new Triangle(depth + 1, P20, P12, sommet[2], density, true, max_depth);
        feuilles[2] = new Triangle(depth + 1, P20, P01, P12, density, true, max_depth);
        feuilles[3] = new Triangle(depth + 1, P01, sommet[1], P12, density, true, max_depth);
    }
}

void Triangle::update() {
    compute_feuille();
    for (int i = 0; i < 4; i++) {
        if (feuilles[i] != nullptr) {
            if (depth < min(density(feuilles[i]->center), max_depth)) {
                children[i] = feuilles[i];
                children[i]->is_feuille = false;
                children[i]->depth = depth + 1;
                children[i]->update();
            }
        } else {
            cout << "is_feuille " << is_feuille << " depth = " << depth << " Erreur: feuille[" << i << "] est nullptr" << endl;
        }
    }
}

// Implémentation des méthodes de la classe Tree
Tree::Tree(function<int(Point)> densFunc, int max_dete) : density(densFunc), max_depth(max_dete) {
    root = nullptr;
}

Tree::~Tree() {
    cout << "Destruction de l'arbre" << endl;
    if (root != nullptr){
        delete root;
    }
}

void Tree::saveToFile(const string& filename) {
    cout << "Sauvegarde des fichiers" << endl;
    ofstream outFile(filename);

    if (!outFile.is_open()) {
        cout << "Erreur d'ouverture du fichier : " << filename << endl;
        return;
    }

    saveTriangle(outFile, root);
    outFile.close();
    cout << "Données sauvegardées dans le fichier : " << filename << endl;
}
//Triangle(int d, Point p0, Point p1, Point p2, function<int(Point)> densFunc = [](Point) { return 0; }, bool is_f = false, int max_dete = 10);
void Tree::createRoot(Point p0, Point p1, Point p2) {
    root = new Triangle(0, p0, p1, p2, density, false, max_depth);
    root->is_root= true;
}

void Tree::saveTriangle(ofstream& outFile, Triangle* triangle) {
    if (triangle == nullptr) {
        return;
    }

    outFile << "Triangle Depth " << triangle->depth << " : ";
    for (int i = 0; i < 3; ++i) {
        outFile << "(" << triangle->sommet[i].x << ", " << triangle->sommet[i].y << ") ";
    }
    outFile << "\n";

    for (Triangle* child : triangle->children) {
        saveTriangle(outFile, child);
    }
}

void Tree::exportPolylines(Polylines* poly, Triangle* triangle){
    if (triangle == nullptr || poly == nullptr) return;
    if (not triangle->is_feuille){
        
        if (triangle->is_root){

            Polyline polyroot;

            polyroot.append(Point(triangle->sommet[0].x, triangle->sommet[0].y)); 
            polyroot.append(Point(triangle->sommet[1].x, triangle->sommet[1].y));
            polyroot.append(Point(triangle->sommet[2].x, triangle->sommet[2].y));

            poly->push_back(polyroot);
        }
        

        if (triangle->children[2] != nullptr)
        {
            Polyline polychild2;
            for (int i = 0; i<3; i++){
                polychild2.append(triangle->children[2]->sommet[i]);
            }
            polychild2.append(triangle->children[2]->sommet[0]);
            poly->push_back(polychild2);
        }
        else{
            for (int i = 0; i < 4; ++i) {
                if (triangle->children[i] != nullptr) {
                    switch (i) {
                        case 0:{
                            Polyline polychild0;
                            polychild0.append(triangle->children[i]->sommet[1]);
                            polychild0.append(triangle->children[i]->sommet[2]);
                            poly->push_back(polychild0);
                            break;
                        }
                        case 1:{
                            Polyline polychild1;
                            polychild1.append(triangle->children[i]->sommet[0]);
                            polychild1.append(triangle->children[i]->sommet[1]);
                            poly->push_back(polychild1);
                            break;
                        }
                        case 3:{
                            Polyline polychild3;
                            polychild3.append(triangle->children[i]->sommet[2]);
                            polychild3.append(triangle->children[i]->sommet[0]);
                            poly->push_back(polychild3);
                            break;
                        }
                    }
                }
            }
        };
        for (int i = 0; i<4; i++){
            if (triangle->children[i] != nullptr){
                exportPolylines(poly, triangle->children[i]);
            }
        }
    }  


}

void Tree::export_polylines_to_txt(const Polylines* polylines, const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        // Écrire chaque polyline dans le fichier avec une séparation entre les polylines
        for (const Polyline& polyline : *polylines) {
            for (const Point& p : polyline.points) {
                file << p.x << " " << p.y << "\n";
            }
            file << "#\n"; // Utilise "#" pour marquer la fin d'une polyline
        }

        file.close();
        std::cout << "Polylines exported to " << filename << std::endl;
    }



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    ////DEBUT DE LA METHODE DE REMPLISSAGE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
Fillpavage::_fill_surface_single(
    unsigned int                    thickness_layers,
    const direction_t               &direction,
    ExPolygon                       &expolygon,
    Polylines*                      polylines_out)
{
    std::lock_guard<std::mutex> lock(console_mutex_pavage);


    if (expolygon.contour.points.empty()) return ; // si le polygone est vide on ne fait rien 

    BoundingBox b_box = expolygon.bounding_box();

    /*
    cout << "taille de l'array bounding box" << sizeof b_box.polygon() <<endl;
    cout<<"bounding box " << b_box.polygon()[0].x <<   " , " <<b_box.polygon()[0].y << "/ "
    <<b_box.polygon()[1].x <<   " , " <<b_box.polygon()[1].y << "/ " 
    << b_box.polygon()[2].x <<   " , " <<b_box.polygon()[2].y << "/ "
    << b_box.polygon()[3].x <<   " , " <<b_box.polygon()[3].y << endl; 
    */
    long hauteur = b_box.polygon()[3].y - b_box.polygon()[0].y;
    long largeur = b_box.polygon()[1].x - b_box.polygon()[0].x;
    cout << "largeur = " << largeur << endl;
    cout << "hauteur  = " << hauteur << endl;
    long x_min = b_box.polygon()[0].x;
    long x_max = b_box.polygon()[1].x;

    //Tree arbre_triangulaire = Tree([](Point) -> int { return 7; }, 10);//Tree(densityFunction_pavage, 5);

    auto variable_density = [x_min, x_max](const Point& P) -> int {
        return int((4.0 / (x_max - x_min)) * (P.x - x_min) + 3);
    };

    Tree arbre_triangulaire = Tree(variable_density, 10);

    arbre_triangulaire.createRoot(Point(static_cast<coord_t>(b_box.polygon()[0].x-hauteur/tan(PI/3)), static_cast<coord_t>(b_box.polygon()[0].y)), 
                                  Point(static_cast<coord_t>(b_box.polygon()[1].x+hauteur/tan(PI/3)), static_cast<coord_t>(b_box.polygon()[1].y)),
                                  Point(static_cast<coord_t>((b_box.polygon()[0].x+b_box.polygon()[1].x)/2), b_box.polygon()[0].y+ static_cast<coord_t>(tan(PI/3)*(largeur/2+hauteur/tan(PI/3)))));

    cout << static_cast<coord_t>(b_box.polygon()[0].x-hauteur/tan(PI/3)) << endl; 
    cout << static_cast<coord_t>(b_box.polygon()[0].y) << endl;
    //arbre_triangulaire.createRoot(Point(0, 0), Point(0, 1), Point(1, 1));
    //Triangle(int d, Point p0, Point p1, Point p2, function<int(Point)> densFunc = [](Point) { return 0; }, bool is_f = false, int max_dete = 10);

    //Triangle t = Triangle(0, Point(1, 0), Point(0, 1), Point(1, 1), arbre_triangulaire.density, false, arbre_triangulaire.max_depth);
    
    arbre_triangulaire.root->update();
    arbre_triangulaire.saveToFile("triangle.txt");
    Polylines path_out;
    cout << "avant export polylines" << endl;
    
    arbre_triangulaire.exportPolylines(&path_out,arbre_triangulaire.root);
    arbre_triangulaire.export_polylines_to_txt(&path_out, "path_out.txt");
    cout << "apres export polylines" << endl;
    path_out = intersection_pl(path_out, expolygon);
    arbre_triangulaire.export_polylines_to_txt(&path_out, "path_out_intersection.txt");
    *polylines_out=path_out;
    
}

}
