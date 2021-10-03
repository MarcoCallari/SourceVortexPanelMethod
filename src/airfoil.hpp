#pragma once
#include <vector>
#include <cstddef>
#include "panel.hpp"

class Airfoil {
public:
    static std::optional<Airfoil> points_into_panels(const std::vector<Point>& points, const int nPanels);
    const std::vector<Panel>& panels() const;
    std::vector<Panel>& panels();
private:
    Airfoil() = default;
    std::vector<Panel> m_panels;
};
