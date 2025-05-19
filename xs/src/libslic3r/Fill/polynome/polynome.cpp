#include <SFML/System/Vector2.hpp>

#include <boost/date_time/time_defs.hpp>
#include <fstream>
#include <vector>
#include <iostream>
#include <glpk.h>
#include <optional> 
#include <utility>
#include <iomanip> 

#include <cmath>


#include "polynome.hpp"
#include "json.hpp"


using json = nlohmann::json;
using Tensor2D = std::vector<std::vector<double>>;
using Tensor3D = std::vector<Tensor2D>;

////////////////// UTILITAIRE //////////////////////////////////////


void print_vector(const std::vector<double>& vec, int precision) {
    std::cout << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << std::fixed << std::setprecision(precision) << vec[i];
        if (i != vec.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

void print_matrix(const std::vector<std::vector<double>>& mat, int precision, int width) {
    for (const auto& row : mat) {
        std::cout << "[";
        for (size_t i = 0; i < row.size(); ++i) {
            std::cout << std::fixed << std::setw(width) << std::setprecision(precision) << row[i];
            if (i != row.size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
}

void print_l_point(const std::vector<std::vector<std::pair<double, double>>>& points, 
                   int precision = 3, 
                   const std::string& point_separator = " -> ")
{
    if(points.empty()) {
        std::cout << "Aucun point à afficher\n";
        return;
    }

    std::cout << std::fixed << std::setprecision(precision);

    for(size_t block_idx = 0; block_idx < points.size(); ++block_idx) {
        const auto& block = points[block_idx];
        
        std::cout << "Bloc " << block_idx + 1 << ":\n";
        
        if(block.empty()) {
            std::cout << "  [vide]\n";
            continue;
        }

        for(size_t i = 0; i < block.size(); ++i) {
            const auto& pt = block[i];
            std::cout << "(" << pt.first << ", " << pt.second << ")";
            
            if(i != block.size() - 1) {
                std::cout << point_separator;
                
                // Retour à la ligne tous les 3 points pour la lisibilité
                if((i + 1) % 3 == 0) std::cout << "\n ";
            }
        }
        std::cout << "\n\n";
    }
}

Tensor2D multiply_by_eta(const Tensor2D& input, double eta_x) {
    Tensor2D result = input;
    for (auto& row : result) {
        for (auto& value : row) {
            value *= eta_x;
        }
    }
    return result;
}

double binomial(int k, int n) {
    if (k < 0 || k > n) return 0;
    return std::tgamma(n + 1) / (std::tgamma(k + 1) * std::tgamma(n - k + 1));
}

double evaluate_polynome(std::vector<double> polynome, double x){
    float value = 0;
    for (int i=0; i< polynome.size(); i++){
        value+= polynome[i]*pow(x, i);
    }
    return value;
}

/////// ACTUAL ALGO ///////////////////////////////////

Tensor3D open_polynome_file(std::string filename) {
    /*
    Ouvre un JSON contenant un polynôme à 3 inconnues sous forme d'un tableau
    et retourne un tenseur 3D (polynome_[i][j][k] correspond au coeff de X^i * Y^j * Z^k)
    */
    Tensor3D polynome_tensor;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier : " + filename);
    }

    json json_container;
    file >> json_container;

    // Récupérer uniquement la partie "tensor"
    polynome_tensor = json_container["tensor"].get<Tensor3D>();

    return polynome_tensor;
}


Tensor2D open_bbox_file(std::string filename) {
    /*
    Ouvre un fichier JSON contenant un polynôme et retourne la bounding box (bbox),
    définie comme une liste de 4 points [x, y].
    */
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier : " + filename);
    }

    json json_container;
    file >> json_container;

    // Extraction de la bbox
    Tensor2D bbox = json_container["bbox"].get<Tensor2D>();
    return bbox;
}

Tensor2D slice_polynome_on_z_axis(Tensor3D polynome, double z){
    /*
    polynome : polynome a 3 inconnues
    -> et value ce polynome en Z = z et retourne ainsi un polynome a 2 inconnues
    */
    Tensor2D slice_polynome;
    for (size_t i = 0; i < polynome.size(); ++i){
        std::vector<double> ligne; 

        for (size_t j = 0; j < polynome[i].size(); j++){
            double coeff = 0;
            
            for (size_t k = 0; k< polynome[i][j].size(); k++){
                coeff += polynome[i][j][k]*std::pow(z, k);
            }
            ligne.push_back(coeff);
        }
        slice_polynome.push_back(ligne);

    }
    return slice_polynome;
}

std::vector<double> evaluate_on_x_axis(Tensor2D polynome, double x) {
    /*
        polynome[i][j] correspond à X^i * Y^j
        retourne la liste des coefficients du polynôme en Y obtenu en évaluant X = x
    */

    if (polynome.empty()) return {};

    // Taille maximale pour Y (nombre de colonnes)
    size_t maxY = 0;
    for (const auto& row : polynome) {
        if (row.size() > maxY) {
            maxY = row.size();
        }
    }

    std::vector<double> result(maxY, 0.0);

    for (size_t i = 0; i < polynome.size(); ++i) {
        double xi = pow(x, i); // x^i
        for (size_t j = 0; j < polynome[i].size(); ++j) {
            result[j] += polynome[i][j] * xi;
        }
    }

    return result;
}

std::vector<double> evaluate_on_y_axis(Tensor2D polynome, double y) {
    /*
        polynome[i][j] correspond à X^i * Y^j
        retourne la liste des coefficients du polynôme en X obtenu en évaluant Y = y
    */

    if (polynome.empty()) return {};

    // Taille maximale pour X (nombre de lignes)
    size_t maxX = polynome.size();

    std::vector<double> result(maxX, 0.0);

    for (size_t i = 0; i < polynome.size(); ++i) {
        for (size_t j = 0; j < polynome[i].size(); ++j) {
            double yj = pow(y, j); // y^j
            result[i] += polynome[i][j] * yj;
        }
    }

    return result;
}

Eigen::MatrixXd Matrice_angle(std::vector<double> angle_vect, int k){
    /*
    Voire la demonstration du theoreme de décomposition : retourne la matrice M_{k}(theta0....theta1)
    */
    int n = angle_vect.size()-1;
    Eigen::MatrixXd mat_transfert(k+1, n+1);

    for (int i=0; i<=k; i++){
        for (int j=0; j<=n; j++){
            mat_transfert(i, j) = pow(std::sin(angle_vect[j]), k-i)*pow(std::cos(angle_vect[j]), i); 
        }
    }
    return mat_transfert;

}

Eigen::VectorXd find_max_min_solution(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    /*
    cherche/ trouve x tel que Ax = b et tel que min (b) soit le plus grand possible
    (on cherche ainsi a ne pas avoir de valeur negative du polynome)
    */
    
    int m = A.rows(); // Nb d'équations
    int n = A.cols(); // Nb de variables

    int total_vars = n + 1;
    int total_rows = m + n; 

    glp_prob* lp = glp_create_prob();
    glp_term_out(GLP_OFF); // on lui demande gentillement de se taire
    glp_set_obj_dir(lp, GLP_MAX); // On maximise t

    glp_add_rows(lp, total_rows);
    // Contraintes A * x = b
    for (int i = 0; i < m; ++i) {
        glp_set_row_bnds(lp, i + 1, GLP_FX, b(i), b(i));
    }


    for (int i = 0; i < n; ++i) {
        glp_set_row_bnds(lp, m + 1 + i, GLP_LO, 0.0, 0.0);
    }

    glp_add_cols(lp, total_vars);
    // Variables x
    for (int j = 0; j < n; ++j) {
        glp_set_col_bnds(lp, j + 1, GLP_FR, 0.0, 0.0);
    }
    // Variable t
    int t_idx = total_vars;
    glp_set_col_bnds(lp, t_idx, GLP_FR, 0.0, 0.0); 
    glp_set_obj_coef(lp, t_idx, 1.0); 

    // Préparer matrice creuse (1-based indexing)
    std::vector<int> ia(1), ja(1); 
    std::vector<double> ar(1);

    // Remplir A * x = b
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            ia.push_back(i + 1);
            ja.push_back(j + 1);
            ar.push_back(A(i, j));
        }
    }

  
    for (int i = 0; i < n; ++i) {
        ia.push_back(m + 1 + i); ja.push_back(i + 1); ar.push_back(1.0);      
        ia.push_back(m + 1 + i); ja.push_back(t_idx);  ar.push_back(-1.0);     
    }

    glp_load_matrix(lp, ia.size() - 1, ia.data(), ja.data(), ar.data());
    glp_simplex(lp, nullptr);

    Eigen::VectorXd x(n);
    for (int i = 0; i < n; ++i) {
        x(i) = glp_get_col_prim(lp, i + 1);
    }

    glp_delete_prob(lp);
    return x;
}



Tensor2D decomposition_polynomial(Tensor2D polynome, std::vector<double> angle_vector){
    /*
    polynome_phi[i][k] correpond au coeff X^k du i éme polynome
    */
    int n = angle_vector.size()-1;
    
    Tensor2D polynome_phi;
    
    for (int i=0; i<=n ;i++){
        polynome_phi.push_back({});
    }
    
    for (int k=0; k<=n; k++){
        Eigen::VectorXd b(k+1);
        for (int i=0; i<=k; i++){
            b(i) = polynome[i][k-i]/binomial(i, k);
        } 
        Eigen::MatrixXd M = Matrice_angle(angle_vector, k);
        Eigen::VectorXd x = find_max_min_solution(M, b);

        
        for (int i=0; i<=n ;i++){
            polynome_phi[i].push_back(x(i));
        }
    }
    
    return polynome_phi;

}

Tensor2D primitive(Tensor2D polynome_phi){
    /*
    Les polynomes sont representer par leur coeffs en ligne 
    dans la matrice polynome_phi
    -> la primitive de ces mêmes polynome dans le meme format
    */
    Tensor2D primitive_phi;
    
    for (int i = 0; i<polynome_phi.size(); i++){
        primitive_phi.push_back({0});
        for (int j = 0; j<polynome_phi[i].size(); j++){
            primitive_phi[i].push_back(polynome_phi[i][j]/(j+1));
        }
    }
    return primitive_phi;
}


Tensor2D rehausse_polynome(Tensor2D poly_phi, std::vector<double> a, std::vector<double> b, int nb_point, double epsilon){
    /*
    poly_phi : vecteur contenant touts s polynome phi : il est calculer le minum 
    du minium des evaluation de chaque polynome sur le segment [a, b] en les evaluant sur 
    nb_point equirepartie.
    Si ce minimum est negatifs chaque coeff constant du polynome est rehausser pour que le polynome soit supérieur a 0 sur [a, b]
    */
    double minimum = 0;
    int nb_poly = poly_phi.size();
    //double pas = (b-a)/nb_point;

    for (int i=0; i<poly_phi.size(); i++){
        std::vector<double> poly = poly_phi[i];
        double step = (b[i]-a[i])/nb_point;
        
        for (double x=a[i]; x<b[i]; x+=step){
        
            double value = evaluate_polynome(poly, x);
            if (value<minimum){
                minimum = value;
            }
        }
    }

    if (minimum < 0){
        for (int i=0; i<nb_poly; i++){
            std::cout<< "on rehausse de " << -minimum << std::endl; 
            poly_phi[i][0] += epsilon - minimum; 
        }
    }
    return poly_phi; 
}
/*
bool solve_polynomial_in_interval_eigen(
    const std::vector<double>& coeffs_std,
    double y,
    double a,
    double b,
    double & result,
    double tol =1e-12
) {
    int deg = coeffs_std.size() - 1;
    while (deg >= 0 && std::abs(coeffs_std[deg]) < tol) {
        --deg;
    }
    if (deg < 0) {
        return false; // Le polynôme est nul, donc aucune solution utile
    }

    Eigen::VectorXd coeffs(deg + 1);
    for (int i = 0; i <= deg; ++i) {
        coeffs[i] = coeffs_std[i];
    }

    // Soustraire y au terme constant
    coeffs[0] -= y;

    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    solver.compute(coeffs);

    for (int i = 0; i < solver.roots().size(); ++i) {
        std::complex<double> root = solver.roots()[i];
        if (std::abs(root.imag()) < tol) {
            double x = root.real();
            if (x >= a && x <= b) {
                result = x;
                return true;
            }
        }
    }

    return false;
}
*/

bool solve_polynomial_in_interval_eigen(
    const std::vector<double>& coeffs_std,
    double y,
    double a,
    double b,
    double & result,
    double tol = 1e-12
) {
    int deg = coeffs_std.size() - 1;

    // Ignore trailing near-zero coefficients (high degrees)
    while (deg >= 0 && std::abs(coeffs_std[deg]) < tol) {
        --deg;
    }

    if (deg < 0) {
        std::cerr << "Polynôme nul ou proche de zéro." << std::endl;
        return false;
    }

    // Crée un vecteur Eigen avec les coefficients valides
    Eigen::VectorXd coeffs(deg + 1);
    for (int i = 0; i <= deg; ++i) {
        coeffs[i] = coeffs_std[i];
    }

    // Soustrait y du terme constant pour résoudre P(x) - y = 0
    coeffs[0] -= y;

    // Vérifie que le polynôme change de signe sur [a, b]
    auto eval_poly = [&](double x) -> double {
        double val = 0;
        double x_pow = 1.0;
        for (int i = 0; i <= deg; ++i) {
            val += coeffs[i] * x_pow;
            x_pow *= x;
        }
        return val;
    };

    double fa = eval_poly(a);
    double fb = eval_poly(b);

    if (std::abs(fa) < tol) {
        result = a;
        return true;
    }
    if (std::abs(fb) < tol) {
        result = b;
        return true;
    }

    if (fa * fb > 0) {
        // Aucun changement de signe → peu probable qu'une racine existe
        std::cerr << "Pas de changement de signe entre a et b." << std::endl;
        return false;
    }

    // Recherche de toutes les racines avec Eigen
    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    solver.compute(coeffs);

    bool found = false;
    double closest_root = std::numeric_limits<double>::quiet_NaN();
    double min_distance = std::numeric_limits<double>::max();

    for (int i = 0; i < solver.roots().size(); ++i) {
        std::complex<double> root = solver.roots()[i];
        if (std::abs(root.imag()) < tol) {
            double x = root.real();
            if (x >= a && x <= b) {
                double dist = std::abs(x - a);  // Choisir le plus proche de a
                if (dist < min_distance) {
                    closest_root = x;
                    min_distance = dist;
                    found = true;
                }
            }
        }
    }

    if (found) {
        result = closest_root;
        return true;
    }

    return false;
}

bool solve_polynomial_in_interval_newtonne(
    const std::vector<double>& coeffs_std,
    double y,
    double a,
    double b,
    double & result,
    double tol=1e-8
) {
    std::cout << "a = " << a << " b = " << b << " y = " << y << std::endl;
    print_vector(coeffs_std);

    const int max_iter = 100;
    const double h = 1e-6; // pour la dérivée numérique

    // Point de départ : milieu de l'intervalle
    double x = 0.5 * (a + b);

    for (int iter = 0; iter < max_iter; ++iter) {
        double fx = evaluate_polynome(coeffs_std, x) - y;

        if (std::abs(fx) < tol) {
            result = x;
            std::cout << "solution x = " << x << std::endl;
            return true;
        }

        double fxh_plus = evaluate_polynome(coeffs_std, x + h);
        double fxh_minus = evaluate_polynome(coeffs_std, x - h);
        double dfx = (fxh_plus - fxh_minus) / (2 * h);

  
        if (std::abs(dfx) < 1e-12) {
            break; // dérivée trop petite, on arrête
        }


        x = x - fx / dfx;

        if (x < a || x > b) {
            break;
        }
    }

    return false; // Pas de solution trouvée dans l'intervalle
}

Tensor2D compute_point_integration(Tensor2D poly_phi, float ratio_infill, std::vector<double> a, std::vector<double> b, float extrusion_witdh){
    /*
    poly_phi : liste des polynomes de découplage
    ratio_infill : pourcentage de remplissage entre 0 et 1  
    a, b : borne des points d'intégration 
    */
    //std::cout << "l308" << std::endl;
    int n = 1;//poly_phi.size()-1; /!\ for curve 
    poly_phi = rehausse_polynome(poly_phi, a, b);
    Tensor2D primitive_phi = primitive(poly_phi);
    Tensor2D tensor_point;

    float all_integrale = 0;
    float longueur_totale = 0;


    //for (int i=0; i<primitive_phi.size(); i++){
    for (int i=0; i<1; i++){  // /!!!!!!!!\ for curve
        all_integrale += evaluate_polynome(primitive_phi[i], b[i]) - evaluate_polynome(primitive_phi[i], a[i]);
        longueur_totale += b[i] - a[i]; 
    }
    float step = all_integrale/(longueur_totale*(n+1)*ratio_infill)*extrusion_witdh;// 10e6 permet de convertir les mm de extrusion_width en coord_t de slic3r 
    //float step = all_integrale/()
    /*
    std::cout << "step = " << step <<  "longueur total = " << longueur_totale <<  std::endl;
    print_matrix(poly_phi);
    std::cout<< "density = " <<ratio_infill << std::endl;
    print_vector(a);
    print_vector(b);
    //step = 20000000000; 
    */
    for (int i=0; i<primitive_phi.size(); i++){
        tensor_point.push_back({a[i]});
        double y = evaluate_polynome(primitive_phi[i], a[i]);

        while (true){
            y += step;
            //print_vector(primitive_phi[i]);
            //std::cout << "P(a[i]) = " << evaluate_polynome(primitive_phi[i], a[i]) << "P(b[i]) = " << evaluate_polynome(primitive_phi[i], b[i]) << std::endl;
            //std::cout << "y = " << y  << " a[i] = " << a[i] << " b[i] = " << b[i]  << "step = " << step << std::endl; 
            double root;
            if (!solve_polynomial_in_interval_eigen(primitive_phi[i], y, a[i], b[i], root)) break;
            //std::cout << "root = " << root << std::endl;
            tensor_point[i].push_back(root);

        }
    }
    //std::cout << "sortie " << std::endl; 
    //print_matrix(tensor_point);
    return tensor_point;

}



std::vector<std::vector<std::pair<double, double>>> extract_blocks_with_xy(
    const Tensor2D& matrix,
    const std::vector<double>& l_axis,
    Axis_maki axis_mode)
{
    std::vector<std::vector<std::pair<double, double>>> blocks;

    if (matrix.empty() || matrix.size() != l_axis.size())
        return blocks;

    // Déterminer la longueur maximale des colonnes
    size_t max_cols = 0;
    for (const auto& row : matrix) {
        max_cols = std::max(max_cols, row.size());
    }

    // Parcours colonne par colonne pour les deux axes
    for (size_t col = 0; col < max_cols; ++col) {
        std::vector<std::pair<double, double>> current_block;
        bool in_block = false;

        for (size_t i = 0; i < matrix.size(); ++i) {
            // Vérifier si la colonne existe dans cette ligne
            const bool has_value = (col < matrix[i].size());

            if (has_value) {
                // Créer la paire selon l'orientation
                const auto point = (axis_mode == Axis_maki::X) 
                    ? std::make_pair(l_axis[i], matrix[i][col])  // X fixe -> (l_axis, valeur)
                    : std::make_pair(matrix[i][col], l_axis[i]); // Y fixe -> (valeur, l_axis)

                current_block.push_back(point);
                in_block = true;
            } else {
                // Fin de bloc détectée
                if (in_block) {
                    blocks.push_back(current_block);
                    current_block.clear();
                    in_block = false;
                }
            }
        }

        // Ajouter le dernier bloc de la colonne
        if (in_block) {
            blocks.push_back(current_block);
        }
    }

    return blocks;
}