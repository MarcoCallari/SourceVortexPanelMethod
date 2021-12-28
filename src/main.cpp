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
#include "stackworkspace.hpp"
#include <gsl/gsl_integration.h>


struct SigmaIntegralParameters {
  const Panel* panelI;
  const Panel* panelJ;
};

//s \in [0,Panel.length()]
inline double sigmaIntegral(const double s, void* params) {
    const auto p = reinterpret_cast<SigmaIntegralParameters*>(params);
    const auto panelI = p->panelI;
    const auto panelJ = p->panelJ;
    const double xi = panelI->controlPoint().x;
    const double yi = panelI->controlPoint().y;
    const double deltai = panelI->delta();
    const auto t = (s / panelJ->length());
    const double xj = (1 - t) * panelJ->startPoint().x + t * panelJ->endPoint().x;
    const double yj = (1 - t) * panelJ->startPoint().y + t * panelJ->endPoint().y;
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

struct VelocityIntegralParameters {
    const Panel* panel;
    Point point;
};

double velocityFunctionX(double s, void * params) {
    const auto p = reinterpret_cast<VelocityIntegralParameters*>(params);
    return velocityIntegrandX(p->point, *p->panel, s);
}

double velocityFunctionY(double s, void * params) {
    const auto p = reinterpret_cast<VelocityIntegralParameters*>(params);
    return velocityIntegrandY(p->point, *p->panel, s);
}

//@TODO: create a workspace entirely on the stack. avoid dynamic allocations!
int main(int argv, char** argc) {
    if (argv < 2) {
        printf("Please provide a file containing the airfoil data! Format must be the same used by the airfoiltools.com website.\n");
        return -1;
    }
    const auto points = ProfileDeserializer::open(argc[1]);
    if (!points) {
      printf("Could not parse file %s.\n", argc[1]);
    }
    constexpr int panelsNumber = 40;
    auto test = Airfoil::points_into_panels(*points, panelsNumber);
    if (!test) {
      printf("Could not create a closed profile using file %s.\n", argc[1]);
    }
    /* Compute sigma values */
    Eigen::Matrix<double, panelsNumber, panelsNumber> A;
    Eigen::Matrix<double, panelsNumber, 1> b;
    std::ofstream file_a("./out_a.txt");
    assert(file_a.good());
    for(int i = 0; i < test->panels().size(); ++i){ 
        for(int j = 0; j < test->panels().size(); ++j) {
            if(i!=j) {
                gsl::StackWorkspace<1000> w;
                SigmaIntegralParameters p { &test->panels()[i], &test->panels()[j] };
                gsl_function integrand;
                integrand.function = sigmaIntegral;
                integrand.params = reinterpret_cast<void*>(&p);
                double result = 0;
                double abserr = 0;
                gsl_integration_qags(&integrand, 0, test->panels()[j].length(), 1.49e-8, 1.49e-8, 1'000, w.ptr(), &result, &abserr);
                A(i, j) = 1 / (2.0 * M_PI) * result;
            }
            else {
                A(i, j)= 0.5;
            }
        }
        b(i) = -1.0*std::cos(test->panels().at(i).delta());
    }
    const auto factorizedMatrix =  A.fullPivLu();
    const auto sigmaMatrix = factorizedMatrix.solve(b);
    /* ---- */
    int i = 0;
    double sum = 0;
    for(auto& panel : test->panels()) {
      panel.setSigma(sigmaMatrix(i, 0));
      sum += sigmaMatrix(i,0);
      ++i;
    }
    file_a << sigmaMatrix << "\n";
    /* Mesh grid */
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
    for (int indexY = 0; indexY < numPoints; ++indexY) {
        for (int indexX = 0; indexX < numPoints; ++indexX) {
            for (const auto& panel : test->panels()) {
                gsl::StackWorkspace<1000> w;
                gsl_function integrand;
                integrand.function = velocityFunctionX;
                VelocityIntegralParameters p = { &panel, Point(x[indexX], y[indexY]) };
                integrand.params = reinterpret_cast<void*>(&p);
                double result = 0;
                double abserr = 0;
                gsl_integration_qags(&integrand, 0, panel.length(), 1.49e-8, 1.49e-8, 1'000, w.ptr(), &result, &abserr);
                u(indexY, indexX) += (panel.sigma() / (2.0 * M_PI)) * result;
                integrand.function = velocityFunctionY;
                gsl_integration_qags(&integrand, 0, panel.length(), 1.49e-8, 1.49e-8, 1'000, w.ptr(), &result, &abserr);
                v(indexY, indexX) += panel.sigma() / (2.0 * M_PI) * result;
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
