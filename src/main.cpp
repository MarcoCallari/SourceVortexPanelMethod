/* ----------
 *               TODO
 * - Costruire matrice velocità [u,v]
 * - Salvare il risultato in un file txt e plottare con python
 * ----------
 */
#ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
#endif
#include <cmath>
//#include <matplot/matplot.h>
#include "panel.hpp"
#include "airfoil.hpp"
#include <functional>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "profiledeserializer.hpp"

//s \in [0,Panel.length()]
double integral(const Panel& panelI, const Panel& panelJ, const double s) {
    const double xi = panelI.controlPoint().x;
    const double yi = panelI.controlPoint().y;
    const double deltai = panelI.delta();
    const double xj = panelJ.startPoint().x + s*std::cos(panelJ.angleFromHorizontal());
    const double yj = panelJ.startPoint().y + s*std::sin(panelJ.angleFromHorizontal());
    return ((xi-xj)*std::cos(deltai)+(yi-yj)*std::sin(deltai)) /
            (std::pow(xi-xj,2)+std::pow(yi-yj,2));
}

//s \in [0,Panel.length()]
double velocityIntegralX(const Point& point, const Panel& panelJ, const double s) {
    const double xi = point.x;
    const double yi = point.y;
    const double xj = panelJ.startPoint().x + s*std::cos(panelJ.angleFromHorizontal());
    const double yj = panelJ.startPoint().y + s*std::sin(panelJ.angleFromHorizontal());
    return (xi-xj) /
            (std::pow(xi-xj,2)+std::pow(yi-yj,2));
}

//s \in [0,Panel.length()]
double velocityIntegralY(const Point& point, const Panel& panelJ, const double s) {
    const double xi = point.x;
    const double yi = point.y;
    const double xj = panelJ.startPoint().x + s*std::cos(panelJ.angleFromHorizontal());
    const double yj = panelJ.startPoint().y + s*std::sin(panelJ.angleFromHorizontal());
    return (yi-yj) /
            (std::pow(xi-xj,2)+std::pow(yi-yj,2));
}

//@TODO move to a Gaussian method, parallelize, use fast delegate
double trapezoidal(const std::function<double(double x)>& func, const double min, const double max, const double step) {
    double result = 0;
    for(double val = min; val < max - step; val += step) {
            result += ((func(val) + func(val+step)) / 2.0)*(step);
    }
    return result;
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& objects) {
    for (const auto& obj : objects) {
        stream << obj << "\n";
    }
    return stream;
}

int main(void) {
    Airfoil test(*ProfileDeserializer::open("/home/noon/C++/integral2/naca0012.dat"), 40);
    /* Compute sigma values */
    Eigen::MatrixXd A;
    A.resize(test.panels().size(), test.panels().size());
    Eigen::MatrixXd b;
    b.resize(test.panels().size(), 1);
    for(int i = 0; i < test.panels().size(); ++i){ 
        for(int j = 0; j < test.panels().size(); ++j) {
            if(i!=j) {
                const auto& panelI = test.panels().at(i);
                const auto& panelJ = test.panels().at(j);
                A(i, j) = 1/(2.0*M_PI) * trapezoidal([&panelI, &panelJ](const double s) { 
                        return integral(panelI, panelJ, s);
                        }, 0.0, panelJ.length(), 1.0/1'000.0); 
            }
            else {
                A(i, j)= 0.5;
            }
        }
        b(i) = -1.0*std::cos(test.panels().at(i).delta());
    }
    const auto factorizedMatrix =  A.fullPivLu();
    const auto sigmaMatrix = factorizedMatrix.solve(b);
    /* ---- */
    int i = 0;
    double sum = 0;
    for(auto& panel : test.panels()) {
      panel.setSigma(sigmaMatrix(i, 0));
      sum += sigmaMatrix(i,0);
      ++i;
    }
    //std::cout << test.panels() << std::endl;
    return 0;
    /* Mesh grid */
    //x \in [-3.0,3.0]
    //y \in [-3.0,3.0]
    constexpr int numPoints = 100; 
    constexpr double xStep = (1.25 - (-0.25)) / (numPoints - 1); 
    constexpr double yStep = (0.2 - (-0.2)) / (numPoints - 1); 
    double x[numPoints];
    double y[numPoints];
    for (int index = 0; index < numPoints; ++index) {
      x[index] = -0.25 + index * xStep; 
    }
    for (int index = 0; index < numPoints; ++index) {
      y[index] = -0.2 + index * yStep; 
    }
    const double stepX = (1.25 - (-0.25)) / 100.0;
    const double stepY = (0.2 - (-0.2)) / 100.0;
    Eigen::MatrixXd u;
    Eigen::MatrixXd v;
    u.resize(numPoints,numPoints);
    v.resize(numPoints,numPoints);
    for (const auto& panel : test.panels()) {
        for (size_t indexX = 0; indexX < numPoints; ++indexX) {
            for (size_t indexY = 0; indexY < numPoints; ++indexY) {
            u(indexX, indexY) =  1.0 + panel.sigma() / (2.0 * M_PI) * trapezoidal([&panel, x, y, indexX, indexY](const double s) {
                    return velocityIntegralX(Point(x[indexX],y[indexY]),panel,s);
                    }, 0.0, panel.length(), 1.0/ 1'000.0);
            v(indexX, indexY) = panel.sigma() / (2.0 * M_PI) * trapezoidal([&panel, x, y, indexX, indexY](const double s) {
                    return velocityIntegralY(Point(x[indexX],y[indexY]),panel,s);
                    }, 0.0, panel.length(), 1.0/ 1'000.0);
            }
        }
    }
    std::cout << "U" << std::endl;
    std::cout << u << std::endl;
    std::cout << "END U" << std::endl;
    std::cout << "V" << std::endl;
    std::cout << v << std::endl;
    std::cout << "END V" << std::endl;
    /* ---- */
    return 0;
}
