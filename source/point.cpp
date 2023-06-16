#include "ma.hpp"
#include <cmath>

Point::Point() {}

Point::Point(ma_real x, ma_real y) {
    this->x = x;
    this->y = y;
}

Point Point::operator+(Point p) {
    return Point(this->x + p.x, this->y + p.y);
}

Point Point::operator-(Point p) {
    return Point(this->x - p.x, this->y - p.y);
}

Point Point::operator*(ma_real s) {
    return Point(this->x * s, this->y * s);
}

Point operator*(ma_real s, Point p) {
    return Point(p.x * s, p.y * s);
}

ma_real Point::norm() {
    return std::sqrt(this->x * this->x + this->y * this->y);
}

ma_real Point::norm2() {
    return this->x * this->x + this->y * this->y;
}

void Point::normalize() {
    ma_real norm = std::sqrt(this->x * this->x + this->y * this->y);
    this->x /= norm;
    this->y /= norm;
}