#ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
#endif
#include <cmath>
#include "panel.hpp"
#include <ostream>

Panel::Panel(const Point& start, const Point& end) :
    m_start (start),
    m_end (end),
    m_control(0,0) {
    m_length = std::sqrt((end.x-start.x)*(end.x-start.x) + (end.y - start.y)*(end.y-start.y));
    m_angleFromHorizontal = std::atan2(end.y-start.y, end.x-start.x);
    m_angleFromHorizontal = std::fmod(std::fmod(m_angleFromHorizontal, M_PI * 2.0) + M_PI*2.0, M_PI*2.0);
    m_delta = m_angleFromHorizontal - M_PI / 2.0;
    m_control = Point(m_start.x + m_length/2.0*std::cos(m_angleFromHorizontal), 
                      m_start.y + m_length/2.0*std::sin(m_angleFromHorizontal));
}

double Panel::delta() const {
    return m_delta;
}

double Panel::length() const {
    return m_length;
}

const Point& Panel::startPoint() const {
    return m_start;
}

const Point& Panel::controlPoint() const {
    return m_control;
}

const Point& Panel::endPoint() const {
    return m_end;
}

double Panel::angleFromHorizontal() const {
    return m_angleFromHorizontal;
}

std::ostream& operator<<(std::ostream& stream, const Panel& panel) {
    stream << "[" << panel.m_start << "," << panel.m_end << "]";
    return stream;
}
