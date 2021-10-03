#pragma once 
#include "point.hpp"
#include <optional>
#include <iosfwd>

class Panel {
public:
    Panel(const Point& start, const Point& end);
    double delta() const;
    double length() const;
    const Point& startPoint() const;
    const Point& controlPoint() const;
    const Point& endPoint() const;
    double angleFromHorizontal() const;
    void setSigma(const double sigma);
    double sigma() const;
private:
    Point m_start;
    Point m_end;
    Point m_control;
    double m_delta;
    double m_angleFromHorizontal;
    double m_length;
    std::optional<double> m_sigma;
    friend std::ostream& operator<<(std::ostream& stream, const Panel& panel);
};
