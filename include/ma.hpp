#ifndef __MA_HPP__
#define __MA_HPP__

#define real double

#include <string>
#include <vector>
#include <omp.h>

#define EPS (1e-10)
#define PI (3.14159265358979323846264)
#define INF (1e18)

#define TWO_SCALE_PREPARE 1
#define KD_TREE 1

#ifdef KD_TREE
    #define MAXIMUM_ELEMENTS_IN_LEAF 20
#endif

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct Point {
    real x;
    real y;
    Point();
    Point(real x, real y);
    Point operator+(Point p);
    Point operator-(Point p);
    Point operator*(real s);
    friend Point operator*(real s, Point p);

    real norm();
    real norm2();
    void normalize();
};

struct Element {
    int v[3];
};

struct Mesh {
    std::vector<Point> points;
    std::vector<Element> elements;
    int num_points;
    int num_elements;
    
#if KD_TREE
    struct TreeNode {
        bool is_leaf;
        int dim; // 0 for x, 1 for y
        real par; // Partition value
        std::vector<int> elements; // Only valid for leaf nodes
        int lson, rson;
    };
    std::vector<TreeNode> tree;

    void build_tree(int node_id, std::vector<int> elements,
                    int dim, real lx, real rx, real ly, real ry);
#endif

    Mesh(std::string filename);
    ~Mesh();


    std::vector<int> N0();
    std::vector<int> Nb();
    bool is_on_boarder(Point p);
    bool is_in_element(Point p, Element e);
    int point_location(Point p);
};

struct MeshFunction {
    std::vector<real> values;
    Mesh& mesh;

    MeshFunction(Mesh& mesh);
    ~MeshFunction();

    real& operator[](int i);
    MeshFunction operator-(MeshFunction& u);

    void init(real (*f)(Point));
    real interpolate(Point p);
    real interpolate(Point p, int id);
    real second_difference(int id, Point v, real delta);
    real second_difference(int id, Point p_plus, int e_plus_id, Point p_minus, int e_minus_id, real rho, real delta);
    real norm_inf();
    real norm2();
};

namespace two_scale {
    std::vector<Point> generate_s_theta(real theta);
    real T_epsilon(MeshFunction& u, int id, std::vector<Point> s_theta, real delta);
    MeshFunction perron(real (*f)(Point), real (*g)(Point), Mesh& mesh, real delta, real theta, int p = -1);
}

#endif