#pragma once

#include <iris_2d/geometry.h>
#include <iris_2d/ccp_solver.h>
#include <iris_2d/ie_solver.h>
#include <vector>

namespace iris_2d
{

class Problem
{
public: 
    Problem() {}
    ~Problem() {}

    bool initialize(const Vector& seed, const std::vector<Obstacle>& obs);
    void addObstacle(const Obstacle& obs) {obs_.push_back(obs);}
    void setObstacle(const std::vector<Obstacle>& obs) {obs_ = obs;}
    bool solve();

    bool separatingPlanes();
    bool inflateRegion(double& volume);
    HyperPlane tangentPlane(const Vector& point, const Matrix& Civn2);
    Region& getRegion() {return region_;}
    Polyhedron& getPolyhedron() {return region_.getPolyhedron();}
    Ellipsoid& getEllipsoid() {return region_.getEllipsoid();}
    PlaneMatrix& getA() {return region_.getPolyhedron().getA();}
    PlaneVector& getB() {return region_.getPolyhedron().getB();}
    Matrix& getC() {return region_.getEllipsoid().getC();}
    Vector& getD() {return region_.getEllipsoid().getD();}
    
private:
    std::vector<Obstacle> obs_;
    Region region_;
    ccp_solver::CcpSovler ccp_solver_;
    ie_solver::IeSovler ie_solver_;

    double best_vol_;
};

}