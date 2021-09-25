#ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
#endif
#include <cmath>
#include "airfoil.hpp"
#include <algorithm>
#include <assert.h>


Airfoil::Airfoil(const std::vector<Point>& originalPoints, const int nPanels) {
    if(!originalPoints.empty()) {
        auto points = originalPoints;
        points.emplace_back(originalPoints.at(0));
        const double maxX = std::max_element(points.cbegin(), points.cend(), 
                [](const auto& point1, const auto& point2) { return point1.x < point2.x; })->x;
        const double maxY = std::max_element(points.cbegin(), points.cend(), 
                [](const auto& point1, const auto& point2) { return point1.x < point2.x; })->y;
        const double minX = std::min_element(points.cbegin(), points.cend(), 
                [](const auto& point1, const auto& point2) { return point1.x < point2.x; })->x;
        const double minY = std::min_element(points.cbegin(), points.cend(), 
                [](const auto& point1, const auto& point2) { return point1.x < point2.x; })->y;
        const double R = (maxX - minX) / 2.0; 
        const double centerX = (maxX + minX) / 2.0; 
        std::vector<double> circleX;
        circleX.reserve(nPanels + 1);
        const double thetaStep = M_PI * 2.0  / static_cast<double>(nPanels);
        for (size_t index = 0; index < nPanels + 1; ++index) {
            circleX.push_back(centerX + R*std::cos(index * thetaStep)); 
        }
        std::vector<double> circleY;
        circleY.reserve(circleX.size());
        for (const double x : circleX) {
            size_t index = 0;
            for(; index < points.size() - 1; ++index) {
                if (((x >= points.at(index).x) && (x <= points.at(index+1).x)) ||
                    ((x >= points.at(index+1).x) && (x <= points.at(index).x)) ){
                    break;
                }
            }
            /*
             * Y - Y0   X - X0
             * ------ = ------
             * Y1- Y0   X1- X0
             */
            circleY.push_back((x-points.at(index).x) / (points.at(index+1).x - points.at(index).x) * (points.at(index+1).y - points.at(index).y) + points.at(index).y);
        }
        m_panels.reserve(nPanels);
        for(size_t index = 0; index < nPanels; ++index) {
            m_panels.emplace_back(Point(circleX.at(index), circleY.at(index)), Point(circleX.at(index+1), circleY.at(index+1))); 
        }
    }
}

const std::vector<Panel>& Airfoil::panels() const {
    return m_panels;
}

std::vector<Panel>& Airfoil::panels() {
    return m_panels;
}
