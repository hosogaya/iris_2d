#pragma once

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <vector>

namespace iris_2d
{

using Matrix = Eigen::Matrix2d;
using PlaneMatrix = Eigen::Matrix<double, -1, 2>;
using Vector = Eigen::Vector2d;
using PlaneVector = Eigen::VectorXd;
using Obstacle = Eigen::Matrix<double, 2, -1>;
using HyperPlane = std::pair<PlaneVector, double>;

const double ELLIPSOID_C_EPSILON = 1e-4;

class Polyhedron
{
public: 
    Polyhedron() {}
    Polyhedron(const PlaneMatrix& A, const PlaneVector& b) 
    : A_(A), b_(b) {}
    ~Polyhedron() {}

    void setA(const PlaneMatrix& A) {A_ = A;}
    void setB(const PlaneVector& b) {b_ = b;}
    PlaneMatrix& getA() {return A_;}
    PlaneVector& getB() {return b_;}

private:
    PlaneMatrix A_;
    PlaneVector b_;
};

class Ellipsoid
{
public:
    Ellipsoid() {}
    Ellipsoid(const Matrix& C, const Vector& d)
    : C_(C), d_(d) {}
    ~Ellipsoid() {}

    void setC(const Matrix& C) {C_ = C;}
    void setD(const Vector& d) {d_ = d;}
    Matrix& getC() {return C_;}
    Vector& getD() {return d_;}
    double getVolume() const {return C_.determinant();}
private:
    Matrix C_;
    Vector d_;
};

class Region
{
public:
    Region() {}
    ~Region() {}

    Polyhedron& getPolyhedron() {return polyhedron_;}
    Ellipsoid& getEllipsoid() {return ellipsoids_;}
private:
    Polyhedron polyhedron_;
    Ellipsoid ellipsoids_;
};

}