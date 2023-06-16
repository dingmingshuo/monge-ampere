#ifndef __MA_HPP__
#define __MA_HPP__

#define ma_real double

#include <string>
#include <vector>
#include <omp.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

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
    ma_real x;
    ma_real y;
    Point();
    Point(ma_real x, ma_real y);
    Point operator+(Point p);
    Point operator-(Point p);
    Point operator*(ma_real s);
    friend Point operator*(ma_real s, Point p);

    ma_real norm();
    ma_real norm2();
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
        ma_real par; // Partition value
        std::vector<int> elements; // Only valid for leaf nodes
        int lson, rson;
    };
    std::vector<TreeNode> tree;

    void build_tree(int node_id, std::vector<int> elements,
                    int dim, ma_real lx, ma_real rx, ma_real ly, ma_real ry);
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
    std::vector<ma_real> values;
    Mesh& mesh;

    MeshFunction(Mesh& mesh);
    ~MeshFunction();

    ma_real& operator[](int i);
    MeshFunction operator-(MeshFunction& u);

    void init(ma_real (*f)(Point));
    ma_real interpolate(Point p);
    ma_real interpolate(Point p, int id);
    ma_real second_difference(int id, Point v, ma_real delta);
    ma_real second_difference(int id, Point p_plus, int e_plus_id, Point p_minus, int e_minus_id, ma_real rho, ma_real delta);
    ma_real norm_inf();
    ma_real norm2();
};

namespace solver {
VectorXd GMRES(const SparseMatrix<double> &A,
               const VectorXd &b, int m, int max_iter,
               double tolerance);
}

namespace two_scale {
    std::vector<Point> generate_s_theta(ma_real theta);
    ma_real T_epsilon(MeshFunction& u, int id, std::vector<Point> s_theta, ma_real delta);
    MeshFunction perron(ma_real (*f)(Point), ma_real (*g)(Point), Mesh& mesh, ma_real delta, ma_real theta, int p = -1);
}

#endif