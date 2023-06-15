#include "ma.hpp"

#include <stdexcept>
#include <algorithm>
#include <cmath>

MeshFunction::MeshFunction(Mesh& mesh): mesh(mesh) {
    this->values.resize(mesh.num_points);
}

MeshFunction::~MeshFunction() {}

MeshFunction MeshFunction::operator-(MeshFunction& u) {
    if (&this->mesh != &u.mesh) {
        throw std::runtime_error("Mesh not equal!");
    }
    MeshFunction ret(this->mesh);
    for (int i = 0; i < (int)values.size(); i++) {
        ret[i] = this->values[i] - u[i];
    }
    return ret;
}

real& MeshFunction::operator[](int i){
    return this->values[i];
}

void MeshFunction::init(real (*f)(Point)) {
    for (int i = 0; i < this->mesh.num_points; i++) {
        this->values[i] = f(this->mesh.points[i]);
    }
}

real MeshFunction::interpolate(Point p) {
    int id = this->mesh.point_location(p);
    if (id == -1) {
        throw std::runtime_error("Point not in mesh!");
    }
    // Locate element
    Element e = this->mesh.elements[id];
    Point v0 = this->mesh.points[e.v[0]];
    Point v1 = this->mesh.points[e.v[1]];
    Point v2 = this->mesh.points[e.v[2]];
    // Compute barycentric coordinates
    real alpha = ((v1.y - v2.y)*(p.x - v2.x) + (v2.x - v1.x)*(p.y - v2.y)) /
        ((v1.y - v2.y)*(v0.x - v2.x) + (v2.x - v1.x)*(v0.y - v2.y));
    real beta = ((v2.y - v0.y)*(p.x - v2.x) + (v0.x - v2.x)*(p.y - v2.y)) /
        ((v1.y - v2.y)*(v0.x - v2.x) + (v2.x - v1.x)*(v0.y - v2.y));
    real gamma = 1.0 - alpha - beta;
    // Interpolate
    return alpha*this->values[e.v[0]] + beta*this->values[e.v[1]] +
        gamma*this->values[e.v[2]];
}

real MeshFunction::interpolate(Point p, int id) {
    // Locate element
    Element e = this->mesh.elements[id];
    Point v0 = this->mesh.points[e.v[0]];
    Point v1 = this->mesh.points[e.v[1]];
    Point v2 = this->mesh.points[e.v[2]];
    // Compute barycentric coordinates
    real alpha = ((v1.y - v2.y)*(p.x - v2.x) + (v2.x - v1.x)*(p.y - v2.y)) /
        ((v1.y - v2.y)*(v0.x - v2.x) + (v2.x - v1.x)*(v0.y - v2.y));
    real beta = ((v2.y - v0.y)*(p.x - v2.x) + (v0.x - v2.x)*(p.y - v2.y)) /
        ((v1.y - v2.y)*(v0.x - v2.x) + (v2.x - v1.x)*(v0.y - v2.y));
    real gamma = 1.0 - alpha - beta;
    // Interpolate
    return alpha*this->values[e.v[0]] + beta*this->values[e.v[1]] +
        gamma*this->values[e.v[2]];
}

real MeshFunction::second_difference(int id, Point v, real delta) {
    real rho = 1;
    Point p = this->mesh.points[id];
    v.normalize();
    // Determine rho; see equation (2.1) in the paper.
    if (v.x != 0) {
        rho = std::min(rho, sgn(v.x) * (1 - p.x) / (delta * v.x));
        rho = std::min(rho, sgn(v.x) * p.x / (delta * v.x));

    }
    if (v.y != 0) {
        rho = std::min(rho, sgn(v.y) * (1 - p.y) / (delta * v.y));
        rho = std::min(rho, sgn(v.y) * p.y / (delta * v.y));
    }

    if (std::fabs(rho) < EPS) {
        throw std::runtime_error("Second order difference on border is undefined!");
    }

    // Compute second order difference
    Point p_plus = p + rho * delta * v;
    Point p_minus = p - rho * delta * v;
    real diff = (this->interpolate(p_plus) - 2 * this->values[id] +
        this->interpolate(p_minus)) / (rho * rho * delta * delta);
    return diff;
}

real MeshFunction::second_difference(int id, Point p_plus, int e_plus_id,
                                     Point p_minus, int e_minus_id, real rho, real delta) {
    real diff = (this->interpolate(p_plus, e_plus_id) - 2 * this->values[id] +
        this->interpolate(p_minus, e_minus_id)) / (rho * rho * delta * delta);
    return diff;
}

real MeshFunction::norm_inf() {
    real ret = 0;
    for (int i = 0; i < (int)values.size(); i++) {
        ret = std::max(ret, std::fabs(this->values[i]));
    }
    return ret;
}

real MeshFunction::norm2() {
    real ret = 0;
    for (int i = 0; i < (int)values.size(); i++) {
        ret += this->values[i] * this->values[i];
    }
    return std::sqrt(ret);
}