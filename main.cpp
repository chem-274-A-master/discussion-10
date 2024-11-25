#include <iostream>
#include <fstream>
#include <cmath>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>

const std::map<std::string, double> atomic_weights = {
    {"H", 1.007825032},
    {"C", 12.011},
    {"N", 14.007},
    {"O", 15.99491462}
};

int main(void)
{
    // O H H
    Eigen::MatrixXd h2o{
       {0.000000000000,  0.000000000000, -0.134503695264}, 
       {0.000000000000, -1.684916670000,  1.067335684736},
       {0.000000000000,  1.684916670000,  1.067335684736},
    };

    // C O H H
    Eigen::MatrixXd form{
        { 0.000025165297,  0.000000000000,  0.144571523302}, 
        {-0.000038305955,  0.000000000000,  1.343510833886},
        { 0.938708677255,  0.000000000000, -0.443151260635},
        {-0.938658598164,  0.000000000000, -0.443084756552}
    };

    std::cout << "Water coordinates" << std::endl;
    std::cout << h2o << std::endl << std::endl;

    std::cout << "Formaldehyde coordinates" << std::endl;
    std::cout << form << std::endl << std::endl;


    return 0;
}
