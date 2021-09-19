/* ----------
 *               TODO
 * -Invertire matrice
 * -Confrontare risultato con circonferenza python
 * -Fare il grafico
 * ----------
 */
#ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
#endif
#include <cmath>
#include "panel.hpp"
#include "airfoil.hpp"
#include <functional>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

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

std::vector<Point> discretize() {
    const int R = 1;
    std::vector<Point> points;
    points.reserve(11);
    const double step = (2*M_PI) / 10; 
    for(int i=0; i<11; ++i) {
        points.emplace_back(R*std::cos(step*i), R*std::sin(step*i));
    }
    return points;
}

int main(void) {
    Airfoil test(discretize());
    Eigen::Matrix<double,10,10> A;
    Eigen::Matrix<double,10,1> b;
    for(int i = 0; i < test.panels().size(); ++i){ 
        for(int j = 0; j < test.panels().size(); ++j) {
            if(i!=j) {
                const auto& panelI = test.panels().at(i);
                const auto& panelJ = test.panels().at(j);
                A(i, j) = 1/(2.0*M_PI) * trapezoidal([&panelI, &panelJ](const double s) { 
                        return integral(panelI, panelJ, s);
                        }, 0.0, panelJ.length(), 1.0/100'000.0); 
            }
            else {
                //@TODO move to init
                A(i, j)= 0.5;
            }
        }
        b[i] = -1.0*std::cos(test.panels().at(i).delta());
    }
    std::cout << A.fullPivLu().solve(b) << std::endl;
    return 0;
}
