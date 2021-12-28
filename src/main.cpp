#ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include "panel.hpp"
#include "airfoil.hpp"
#include "profiledeserializer.hpp"
#include <gsl/gsl_integration.h>


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

struct VelocityParameters {
    const Panel* panel;
    const double* x;
    const double* y;
    size_t indexX;
    size_t indexY;
};

double velocityFunctionX(double s, void * params) {
    const auto p = reinterpret_cast<VelocityParameters*>(params);
    return velocityIntegrandX(Point(p->x[p->indexX], p->y[p->indexY]), *p->panel, s);
}

double velocityFunctionY(double s, void * params) {
    const auto p = reinterpret_cast<VelocityParameters*>(params);
    return velocityIntegrandY(Point(p->x[p->indexX], p->y[p->indexY]), *p->panel, s);
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
    double maxerrX = 0.0;
    double maxerrY = 0.0;
    //#pragma omp parallel for 
    for (size_t indexY = 0; indexY < numPoints; ++indexY) {
        for (size_t indexX = 0; indexX < numPoints; ++indexX) {
            for (const auto& panel : test.panels()) {
                gsl_function integrand;
                integrand.function = velocityFunctionX;
                VelocityParameters p;
                p.panel = &panel;
                p.x = x;
                p.y = y;
                p.indexX = indexX;
                p.indexY = indexY;
                integrand.params = reinterpret_cast<void*>(&p);
                double result = 0;
                double abserr = 0;
                size_t calls = 0;
                gsl_integration_qng(&integrand, 0, panel.length(), 100, 100, &result, &abserr, &calls);
                u(indexY, indexX) += (panel.sigma() / (2.0 * M_PI)) * result;
                if (abserr > maxerrX) maxerrX = abserr;
                assert(abserr < 1);
                integrand.function = velocityFunctionY;
                gsl_integration_qng(&integrand, 0, panel.length(), 100, 100, &result, &abserr, &calls);
                if (abserr > maxerrY) maxerrY = abserr;
                assert(abserr < 1);
                v(indexY, indexX) += panel.sigma() / (2.0 * M_PI) * result;
            }
        }
    }
    printf("Max err v_x : %f, max err v_y : %f\n", maxerrX, maxerrY);
    std::ofstream file_u("./out_u.txt");
    assert(file_u.good());
    file_u << u;
    std::ofstream file_v("./out_v.txt");
    assert(file_v.good());
    file_v << v;
    /* ---- */
    return 0;
}
