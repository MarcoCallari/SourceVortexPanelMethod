#include "point.hpp"
#include <iostream>

std::ostream& operator<<(std::ostream& stream, const Point& point) {
    stream << "[" << point.x << "," << point.y << "]";
    return stream;
}
