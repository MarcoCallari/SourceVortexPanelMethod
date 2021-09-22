/* ----------
 *               TODO
 *
 *
 *
 * ----------
 */
#ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
#endif
#include <cmath>
// #include <matplot/matplot.h>
#include "panel.hpp"
#include "airfoil.hpp"
#include <functional>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "profiledeserializer.hpp"
#include <assert.h>

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

// std::vector<Point> discretize() {
    // constexpr int R = 1;
    // std::vector<Point> points;
    // points.reserve(11);
    // const double step = (2*M_PI) / 10; 
    // for(int i=0; i<11; ++i) {
        // points.emplace_back(R*std::cos(step*i), R*std::sin(step*i));
    // }
    // return points;
// }

int main(void) {
    const auto points = ProfileDeserializer::open("/home/noon/C++/integral2/naca0012.dat");
    std::cout << *points << std::endl;
    return 0;
}
// int main(void) {
    // Airfoil test(discretize());
    // /* Compute sigma values */
    // Eigen::Matrix<double,10,10> A;
    // Eigen::Matrix<double,10,1> b;
    // for(int i = 0; i < test.panels().size(); ++i){ 
        // for(int j = 0; j < test.panels().size(); ++j) {
            // if(i!=j) {
                // const auto& panelI = test.panels().at(i);
                // const auto& panelJ = test.panels().at(j);
                // A(i, j) = 1/(2.0*M_PI) * trapezoidal([&panelI, &panelJ](const double s) { 
                        // return integral(panelI, panelJ, s);
                        // }, 0.0, panelJ.length(), 1.0/1'000.0); 
            // }
            // else {
                // A(i, j)= 0.5;
            // }
        // }
        // b[i] = -1.0*std::cos(test.panels().at(i).delta());
    // }
    // const auto sigmaMatrix = A.fullPivLu().solve(b);
    // /* ---- */
    // int i = 0;
    // for(auto& panel : test.panels()) {
        // panel.setSigma(sigmaMatrix[i]);
        // ++i;
    // }
    // /* Mesh grid */
    // //x \in [-3.0,3.0]
    // //y \in [-3.0,3.0]
    // constexpr int numPoints = 100; 
    // constexpr double xStep = (3.0 - (-3.0)) / (numPoints - 1); 
    // constexpr double yStep = (3.0 - (-3.0)) / (numPoints - 1); 
    // double x[numPoints];
    // double y[numPoints];
    // for (int index = 0; index < numPoints; ++index) {
        // x[index] = -3.0 + index * xStep; 
    // }
    // for (int index = 0; index < numPoints; ++index) {
        // y[index] = -3.0 + index * yStep; 
    // }
    // /* ---- */
    // auto xM = matplot::linspace(-3.0, 3.0, 100);
    // auto yM = matplot::linspace(-3.0, 3.0, 100);
    // auto [X, Y] = matplot::meshgrid(xM, yM);
    // /* Compute x velocity */
    // /*{
        // double vX[numPoints][numPoints];
        // std::fill(*vX, *vX + numPoints*numPoints, 1); 
        // for (int indexY = 0; indexY < numPoints; ++indexY) {
            // std::cout << "[";
            // for (int indexX = 0; indexX < numPoints; ++indexX) {
                // for (const auto& panel : test.panels()) {
                    // vX[indexY][indexX] += panel.sigma() / (2.0*M_PI) * trapezoidal([&panel,&x,&y,indexX,indexY](const double s) { 
                            // return velocityIntegralX(Point(x[indexX], y[indexY]), panel, s); }, 0.0, panel.length(), 1.0/1'000); 
                // }
                // std::cout << vX[indexY][indexX] << ",";
            // }
            // std::cout << "]" << std::endl;
        // }
    // }*/
    // /* ---- */
    // /* Compute y velocity */
    // /*{
        // double vY[numPoints][numPoints];
        // std::fill(*vY, *vY + numPoints*numPoints, 1); 
        // for (int indexY = 0; indexY < numPoints; ++indexY) {
            // std::cout << "[";
            // for (int indexX = 0; indexX < numPoints; ++indexX) {
                // for (const auto& panel : test.panels()) {
                    // vY[indexY][indexX] += panel.sigma() / (2.0*M_PI) * trapezoidal([&panel,&x,&y,indexX,indexY](const double s) { 
                            // return velocityIntegralY(Point(x[indexX], y[indexY]), panel, s); }, 0.0, panel.length(), 1.0/1'000); 
                // }
                // std::cout << vY[indexY][indexX] << ",";
            // }
            // std::cout << "]" << std::endl;
        // }
    // }*/
    // auto u = matplot::transform(X,Y, [&test](const double x, const double y) { 
            // double u = 1; // V_inf
            // for (const auto& panel : test.panels()) {
                // u += panel.sigma() / (2.0 * M_PI) * trapezoidal([&panel, x,y](const double s) {
                        // return velocityIntegralX(Point(x,y),panel,s);
                        // }, 0.0, panel.length(), 1.0/ 1'000.0);
                // }
            // return  u; });
    // auto v = matplot::transform(X,Y, [&test](const double x, const double y) { 
            // double v = 0; // V_inf
            // for (const auto& panel : test.panels()) {
                // v += panel.sigma() / (2.0 * M_PI) * trapezoidal([&panel, x,y](const double s) {
                        // return velocityIntegralY(Point(x,y),panel,s);
                        // }, 0.0, panel.length(), 1.0/ 1'000.0);
                // }
            // return  v; });
    // matplot::quiver(X,Y,u, v);
    // matplot::show();
    // /* ---- */
    // return 0;
// }
