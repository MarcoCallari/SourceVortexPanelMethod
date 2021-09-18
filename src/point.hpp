#pragma once
#include <iosfwd>

struct Point {
    Point(const double x, const double y) : x(x), y(y) {}
    double x = 0.0;
    double y = 0.0;
    friend std::ostream& operator<<(std::ostream& stream, const Point& point);
};
