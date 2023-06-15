#include "ma.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

Mesh::Mesh(std::string filename) {
    // Read mesh from .obj file
    std::ifstream in_file(filename);
    if (!in_file.is_open()) {
        throw std::runtime_error("Unable to open file " + filename);
    }

    std::string line;
    while (std::getline(in_file, line)) {
        if (line[0] == 'v') {
            std::istringstream iss(line.substr(2));
            Point p;
            real z;
            // Ignore z coordinate
            iss >> p.x >> p.y >> z;
            points.push_back(p);
        } else if (line[0] == 'f') {
            std::istringstream iss(line.substr(2));
            Element e;
            iss >> e.v[0] >> e.v[1] >> e.v[2];
            // Indices start at 1 in .obj files
            e.v[0]--;
            e.v[1]--;
            e.v[2]--;
            elements.push_back(e);
        }
    }
    this->num_points = (int)points.size();
    this->num_elements = (int)elements.size();
#if KD_TREE
    // Build kd-tree
    std::vector<int> all_elements;
    for (int i = 0; i < this->num_elements; i++) {
        all_elements.push_back(i);
    }
    // Initilize root node
    this->tree.push_back(TreeNode());
    this->build_tree(0, all_elements, 0, 0, 1, 0, 1);
#endif
}

Mesh::~Mesh() {
    // Nothing to do here
}

void Mesh::build_tree(int node_id, std::vector<int> elements,
                      int dim, real lx, real rx, real ly, real ry) {
    this->tree[node_id].dim = dim;
    if ((int)elements.size() < MAXIMUM_ELEMENTS_IN_LEAF) {
        this->tree[node_id].is_leaf = true;
        this->tree[node_id].elements.assign(elements.begin(), elements.end());
        return;
    }
    this->tree[node_id].is_leaf = false;
    this->tree[node_id].par = (dim == 0) ? (lx + rx) / 2 : (ly + ry) / 2;
    this->tree.push_back(TreeNode());
    this->tree.push_back(TreeNode());
    this->tree[node_id].lson = (int)this->tree.size() - 2;
    this->tree[node_id].rson = (int)this->tree.size() - 1;
    if (this->tree[node_id].dim == 0) {
        std::vector<int> l_elements, r_elements;
        for (int i = 0; i < (int)elements.size(); i++) {
            Element e = this->elements[elements[i]];
            Point p1 = this->points[e.v[0]];
            Point p2 = this->points[e.v[1]];
            Point p3 = this->points[e.v[2]];
            if (p1.x < this->tree[node_id].par &&
                p2.x < this->tree[node_id].par &&
                p3.x < this->tree[node_id].par) {
                l_elements.push_back(elements[i]);
            } else if (p1.x > this->tree[node_id].par &&
                       p2.x > this->tree[node_id].par &&
                       p3.x > this->tree[node_id].par) {
                r_elements.push_back(elements[i]);
            } else {
                l_elements.push_back(elements[i]);
                r_elements.push_back(elements[i]);
            }
        }
        this->build_tree(this->tree[node_id].lson, l_elements,
                            1 - dim, lx, this->tree[node_id].par, ly, ry);
        this->build_tree(this->tree[node_id].rson, r_elements,
                            1 - dim, this->tree[node_id].par, rx, ly, ry);
    } else {
        std::vector<int> l_elements, r_elements;
        for (int i = 0; i < (int)elements.size(); i++) {
            Element e = this->elements[elements[i]];
            Point p1 = this->points[e.v[0]];
            Point p2 = this->points[e.v[1]];
            Point p3 = this->points[e.v[2]];
            if (p1.y < this->tree[node_id].par &&
                p2.y < this->tree[node_id].par &&
                p3.y < this->tree[node_id].par) {
                l_elements.push_back(elements[i]);
            } else if (p1.y > this->tree[node_id].par &&
                       p2.y > this->tree[node_id].par &&
                       p3.y > this->tree[node_id].par) {
                r_elements.push_back(elements[i]);
            } else {
                l_elements.push_back(elements[i]);
                r_elements.push_back(elements[i]);
            }
        }
        this->build_tree(this->tree[node_id].lson, l_elements,
                            1 - dim, lx, rx, ly, this->tree[node_id].par);
        this->build_tree(this->tree[node_id].rson, r_elements,
                            1 - dim, lx, rx, this->tree[node_id].par, ry);
    }
}

std::vector<int> Mesh::N0() {
    // Return point id on the boarderline.
    std::vector<int> ret;
    for (int i = 0; i < this->num_points; i++) {
        if (!this->is_on_boarder(this->points[i])) {
            ret.push_back(i);
        }
    }
    return ret;
}

std::vector<int> Mesh::Nb() {
    // Return point id on the boarderline.
    std::vector<int> ret;
    for (int i = 0; i < this->num_points; i++) {
        if (this->is_on_boarder(this->points[i])) {
            ret.push_back(i);
        }
    }
    return ret;
}

bool Mesh::is_on_boarder(Point p) {
    return std::fabs(p.x) < EPS || std::fabs(p.y) < EPS ||
           std::fabs(1 - p.x) < EPS || std::fabs(1 - p.y) < EPS;
}

bool Mesh::is_in_element(Point p, Element e) {
    // Check if p is in e
    Point v0 = this->points[e.v[0]];
    Point v1 = this->points[e.v[1]];
    Point v2 = this->points[e.v[2]];
    // Compute barycentric coordinates
    real alpha = ((v1.y - v2.y)*(p.x - v2.x) + (v2.x - v1.x)*(p.y - v2.y)) /
        ((v1.y - v2.y)*(v0.x - v2.x) + (v2.x - v1.x)*(v0.y - v2.y));
    real beta = ((v2.y - v0.y)*(p.x - v2.x) + (v0.x - v2.x)*(p.y - v2.y)) /
        ((v1.y - v2.y)*(v0.x - v2.x) + (v2.x - v1.x)*(v0.y - v2.y));
    real gamma = 1.0 - alpha - beta;
    // Check if p is in e
    return (alpha >= 0.0 - EPS && alpha <= 1.0 + EPS &&
            beta >= 0.0 - EPS && beta <= 1.0 + EPS &&
            gamma >= 0.0 - EPS && gamma <= 1.0 + EPS);
}

int Mesh::point_location(Point p) {
#if KD_TREE
    // Search the tree
    int node_id = 0;
    while (this->tree[node_id].is_leaf == false) {
        if (this->tree[node_id].dim == 0) {
            if (p.x < this->tree[node_id].par) {
                node_id = this->tree[node_id].lson;
            } else {
                node_id = this->tree[node_id].rson;
            }
        } else {
            if (p.y < this->tree[node_id].par) {
                node_id = this->tree[node_id].lson;
            } else {
                node_id = this->tree[node_id].rson;
            }
        }
    }
    // Check if x is inside any element
    int id = -1;
    for (int i = 0; i < (int)this->tree[node_id].elements.size(); i++) {
        if (this->is_in_element(p, this->elements[this->tree[node_id].elements[i]])) {
            id = this->tree[node_id].elements[i];
        }
    }
    // If not in this area, return -1
    return id;
#else
    // Check if x is inside any element
    int id = -1;
    for (int i = 0; i < this->num_elements; i++) {
        if (this->is_in_element(p, this->elements[i])) {
            id = i;
        }
    }
    // If not in this area, return -1
    return id;
#endif
}