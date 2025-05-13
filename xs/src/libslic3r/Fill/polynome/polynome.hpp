#ifndef POLYNOME_HPP_
#define POLYNOME_HPP_

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>
#include <complex>
#include <iomanip>



using Tensor2D = std::vector<std::vector<double>>;
using Tensor3D = std::vector<Tensor2D>;

void print_vector(const std::vector<double>& vec, int precision = 6);
void print_matrix(const std::vector<std::vector<double>>& mat, int precision = 6, int width = 8); 

Tensor3D open_polynome_file(std::string filename);
Tensor2D slice_polynome_on_z_axis(Tensor3D polynome, double z);

Eigen::MatrixXd Matrice_angle(std::vector<double> angle_vect, int k);
Eigen::VectorXd find_max_min_solution(const Eigen::MatrixXd& A, const Eigen::VectorXd& b); 

Tensor2D decomposition_polynomial(Tensor2D polynome, std::vector<double> angle_vector); 
Tensor2D primitive(Tensor2D polynome_phi);
Tensor2D rehausse_polynome(Tensor2D poly_phi, std::vector<double> a, std::vector<double> b, int nb_point=100, double epsilon=0.1); 

bool solve_polynomial_in_interval_eigen(
    const std::vector<double>& coeffs_std,
    double y,
    double a,
    double b,
    double & result,  
    double tol = 1e-8
);

Tensor2D compute_point_integration(Tensor2D poly_phi, float ratio_infill,  std::vector<double> a, std::vector<double> b, float extrusion_witdh=0.4); 

#endif //POLYNOME_HPP