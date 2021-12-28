#ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
#endif
#include <cmath>
#include "airfoil.hpp"
#include <algorithm>
#include <assert.h>
#include <iostream>


std::optional<Airfoil> Airfoil::points_into_panels(const std::vector<Point>& originalPoints, const size_t nPanels) {
    if (originalPoints.empty()) {
      return {};
    }
    auto points = originalPoints;
    points.emplace_back(originalPoints.at(0));
    const double maxX = std::max_element(points.cbegin(), points.cend(), 
            [](const auto& point1, const auto& point2) { return point1.x < point2.x; })->x;
    const double minX = std::min_element(points.cbegin(), points.cend(), 
            [](const auto& point1, const auto& point2) { return point1.x < point2.x; })->x;
    const double R = (maxX - minX) / 2.0; 
    const double centerX = (maxX + minX) / 2.0; 
    std::vector<double> circleX;
    circleX.reserve(static_cast<size_t>(nPanels) + 1);
    const double thetaStep = M_PI * 2.0  / static_cast<double>(nPanels);
    //@TODO parallelize
    for (size_t index = 0; index < static_cast<size_t>(nPanels) + 1; ++index) {
        circleX.push_back(centerX + R*std::cos(index * thetaStep)); 
    }
    std::vector<double> circleY;
    //@TODO clear up this mess
    size_t index = 0;
    circleY.reserve(circleX.size());
    for (const double x : circleX) {
        for(; index < points.size() - 1; ++index) {
            if (((points.at(index).x <= x) && (x <= points.at(index+1).x)) ||
                ((x >= points.at(index+1).x) && (x <= points.at(index).x)) ){
                break;
            }
        }
        circleY.push_back((x-points.at(index).x) / (points.at(index+1).x - points.at(index).x) * (points.at(index+1).y - points.at(index).y) + points.at(index).y);
    }
    Airfoil foil;
    foil.m_panels.reserve(nPanels);
    for(size_t pointIndex = 0; pointIndex < nPanels; ++pointIndex) {
        foil.m_panels.emplace_back(Point(circleX.at(pointIndex), circleY.at(pointIndex)), Point(circleX.at(pointIndex+1), circleY.at(pointIndex+1))); 
    }
    return foil;
}

const std::vector<Panel>& Airfoil::panels() const {
    return m_panels;
}

std::vector<Panel>& Airfoil::panels() {
    return m_panels;
}
