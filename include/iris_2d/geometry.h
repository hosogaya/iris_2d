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
using HyperPlane = std::pair<Vector, double>;

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
    const PlaneMatrix& getA() const {return A_;}
    const PlaneVector& getB() const {return b_;}
    PlaneMatrix& getARef() {return A_;}
    PlaneVector& getBRef() {return b_;}

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
    const Matrix& getC() const {return C_;}
    const Vector& getD() const {return d_;}

    Matrix& getCRef() {return C_;}
    Vector& getDRef() {return d_;}
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

    const Polyhedron& getPolyhedron() const {return polyhedron_;}
    const Ellipsoid& getEllipsoid() const {return ellipsoids_;}
    Polyhedron& getPolyhedronRef() {return polyhedron_;}
    Ellipsoid& getEllipsoidRef() {return ellipsoids_;}

    const PlaneMatrix& getA() const {return polyhedron_.getA();}
    const PlaneVector& getB() const {return polyhedron_.getB();}
    const Matrix& getC() const {return ellipsoids_.getC();}
    const Vector& getD() const {return ellipsoids_.getD();}

    PlaneMatrix& getARef() {return polyhedron_.getARef();}
    PlaneVector& getBRef() {return polyhedron_.getBRef();}
    Matrix& getCRef() {return ellipsoids_.getCRef();}
    Vector& getDRef() {return ellipsoids_.getDRef();}

private:
    Polyhedron polyhedron_;
    Ellipsoid ellipsoids_;
};

}