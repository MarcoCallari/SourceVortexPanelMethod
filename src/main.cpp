#ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
#endif
#include <cmath>
#include "panel.hpp"
#include "airfoil.hpp"
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "profiledeserializer.hpp"
#include <fstream>


//s \in [0,Panel.length()]
inline double integral(const Panel& panelI, const Panel& panelJ, const double s) {
    const double xi = panelI.controlPoint().x;
    const double yi = panelI.controlPoint().y;
    const double deltai = panelI.delta();
    const auto t = (s / panelJ.length());
    const double xj = (1 - t) * panelJ.startPoint().x + t * panelJ.endPoint().x;
    const double yj = (1 - t) * panelJ.startPoint().y + t * panelJ.endPoint().y;
    return ((xi-xj)*std::cos(deltai)+(yi-yj)*std::sin(deltai)) /
            ((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
}

//s \in [0,Panel.length()]
inline double velocityIntegrandX(const Point& point, const Panel& panelJ, const double s) {
    const double xi = point.x;
    const double yi = point.y;
    const auto t = (s / panelJ.length());
    const double xj = (1 - t) * panelJ.startPoint().x + t * panelJ.endPoint().x;
    const double yj = (1 - t) * panelJ.startPoint().y + t * panelJ.endPoint().y;
    return (xi-xj) /
            ((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
}

//s \in [0,Panel.length()]
inline double velocityIntegrandY(const Point& point, const Panel& panelJ, const double s) {
    const double xi = point.x;
    const double yi = point.y;
    const auto t = (s / panelJ.length());
    const double xj = (1 - t) * panelJ.startPoint().x + t * panelJ.endPoint().x;
    const double yj = (1 - t) * panelJ.startPoint().y + t * panelJ.endPoint().y;
    return (yi-yj) /
            ((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
}

template <typename F>
inline double trapezoidal(const F& func, const double min, const double max, const double step) {
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
    const auto points = *ProfileDeserializer::open("/home/noon/C++/integral2/naca0012.dat");
    constexpr int panelsNumber = 40;
    auto test = *Airfoil::points_into_panels(points, panelsNumber);
    /* Compute sigma values */
    Eigen::Matrix<double, panelsNumber, panelsNumber> A;
    Eigen::Matrix<double, panelsNumber, 1> b;
    std::ofstream file_a("./out_a.txt");
    assert(file_a.good());
    for(size_t i = 0; i < test.panels().size(); ++i){ 
        for(size_t j = 0; j < test.panels().size(); ++j) {
            if(i!=j) {
                const auto& panelI = test.panels().at(i);
                const auto& panelJ = test.panels().at(j);
                A(i, j) = 1/(2.0*M_PI) * trapezoidal([&panelI, &panelJ](const double s) { 
                        return integral(panelI, panelJ, s);
                        }, 0.0, panelJ.length(), 1.0/10'000.0); 
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
    file_a << sigmaMatrix << "\n";
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
    Eigen::Matrix<double, numPoints, numPoints> u;
    u.setOnes(numPoints, numPoints);
    Eigen::Matrix<double, numPoints, numPoints> v;
    v.setZero(numPoints, numPoints);
    #pragma omp parallel for 
    for (size_t indexY = 0; indexY < numPoints; ++indexY) {
        for (size_t indexX = 0; indexX < numPoints; ++indexX) {
            for (const auto& panel : test.panels()) {
                const auto velIntX = [&panel, x, y, indexX, indexY](const double s) {
                        return velocityIntegrandX(Point(x[indexX],y[indexY]),panel,s); };
                u(indexY, indexX) += (panel.sigma() / (2.0 * M_PI)) * trapezoidal(velIntX, 0, panel.length(), 1.0/10'000);
                const auto velIntY = [&panel, x, y, indexX, indexY](const double s) {
                        return velocityIntegrandY(Point(x[indexX],y[indexY]),panel,s); };
                v(indexY, indexX) += panel.sigma() / (2.0 * M_PI) *  trapezoidal(velIntY, 0, panel.length(), 1.0/10'000);
            }
        }
    }
    std::ofstream file_u("./out_u.txt");
    assert(file_u.good());
    file_u << u;
    std::ofstream file_v("./out_v.txt");
    assert(file_v.good());
    file_v << v;
    /* ---- */
    return 0;
}
