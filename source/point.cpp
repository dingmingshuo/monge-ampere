#include "ma.hpp"
#include <cmath>

Point::Point() {}

Point::Point(real x, real y) {
    this->x = x;
    this->y = y;
}

Point Point::operator+(Point p) {
    return Point(this->x + p.x, this->y + p.y);
}

Point Point::operator-(Point p) {
    return Point(this->x - p.x, this->y - p.y);
}

Point Point::operator*(real s) {
    return Point(this->x * s, this->y * s);
}

Point operator*(real s, Point p) {
    return Point(p.x * s, p.y * s);
}

real Point::norm() {
    return std::sqrt(this->x * this->x + this->y * this->y);
}

real Point::norm2() {
    return this->x * this->x + this->y * this->y;
}

void Point::normalize() {
    real norm = std::sqrt(this->x * this->x + this->y * this->y);
    this->x /= norm;
    this->y /= norm;
}