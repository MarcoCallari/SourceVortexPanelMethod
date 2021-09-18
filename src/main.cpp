#include <cmath>
#include "panel.hpp"
#include <functional>
#include <vector>
#include "airfoil.hpp"
#include <iostream>
#include <ostream>

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

int main(void) {
    std::vector<Point> points {
        Point(0,0),
        Point(1,1),
        Point(2,1),
        Point(3,0),
        Point(2,-1),
        Point(1,-1),
        Point(0,0)};
    Airfoil test(points);
    const auto& firstPanel = test.panels().front();
    std::cout << test.panels() << std::endl;
    for(const auto & panel : test.panels()) { 
        std::cout << panel.delta() * 180.0 / M_PI << std::endl;
    }
    for(const auto & panel : test.panels()) { 
        std::cout << "[" << std::cos(panel.delta()) << "," << std::sin(panel.delta()) << "]" << std::endl;
    }
    std::cout << "INTEGRALS" << std::endl;
    for(const auto & panel : test.panels()) {
        std::cout << trapezoidal([&panel, &firstPanel](const double s) { return integral(firstPanel, panel, s);}, 0.0, panel.length(), 1.0/100000.0) << std::endl; 
    }
    return 0;
}
