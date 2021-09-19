#pragma once
#include <vector>
#include <cstddef>
#include "panel.hpp"

class Airfoil {
public:
    Airfoil(const std::vector<Point>& points);
    const std::vector<Panel>& panels() const;
    std::vector<Panel>& panels();
private:
    std::vector<Panel> m_panels;
};
