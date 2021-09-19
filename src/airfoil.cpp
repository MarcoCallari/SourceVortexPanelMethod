#include "airfoil.hpp"

Airfoil::Airfoil(const std::vector<Point>& points) {
    m_panels.reserve(points.size() - 1);
    for(size_t index = 0; index < points.size() - 1; ++index) {
        m_panels.emplace_back(points.at(index), points.at(index+1)); 
    }
}

const std::vector<Panel>& Airfoil::panels() const {
    return m_panels;
}

std::vector<Panel>& Airfoil::panels() {
    return m_panels;
}
