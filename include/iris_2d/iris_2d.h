#pragma once

#include <iris_2d/geometry.h>
#include <iris_2d/ccp_solver.h>
#include <iris_2d/ie_solver.h>
#include <vector>
#include <chrono>
#include <numeric>

namespace iris_2d
{

class Problem
{
public: 
    Problem() {}
    ~Problem() {}

    void reset();
    void initialize(const Vector& seed, const std::vector<Obstacle>& obs);
    void addObstacle(const Obstacle& obs) {obs_.push_back(obs);}
    void setObstacle(const std::vector<Obstacle>& obs) {obs_ = obs;}
    void setSeedPos(const Vector& seed) {getEllipsoidRef().setD(seed);}
    bool solve();

    bool separatingPlanes();
    bool inflateRegion(double& volume);
    HyperPlane tangentPlane(const Vector& point, const Matrix& Civn2);

    const Region& getRegion() const {return region_;}
    Polyhedron& getPolyhedronRef() {return region_.getPolyhedronRef();}
    Ellipsoid& getEllipsoidRef() {return region_.getEllipsoidRef();}
    PlaneMatrix& getARef() {return region_.getPolyhedronRef().getARef();}
    PlaneVector& getBRef() {return region_.getPolyhedronRef().getBRef();}
    Matrix& getCRef() {return region_.getEllipsoidRef().getCRef();}
    Vector& getDRef() {return region_.getEllipsoidRef().getDRef();}

    const Polyhedron& getPolyhedron() const {return region_.getPolyhedron();}
    const Ellipsoid& getEllipsoid() const {return region_.getEllipsoid();}
    const PlaneMatrix& getA() const {return region_.getPolyhedron().getA();}
    const PlaneVector& getB() const {return region_.getPolyhedron().getB();}
    const Matrix& getC() const {return region_.getEllipsoid().getC();}
    const Vector& getD() const {return region_.getEllipsoid().getD();}

    double getIeTime() const {return mean_ie_time_;}
    double getCcpTime() const {return mean_ccp_time_;}
    int getIteration() const {return iter_;}

    void setMaxIteration(const int iter) {max_iteration_ = iter;}
private:
    std::vector<Obstacle> obs_;
    Region region_;
    ccp_solver::CcpSovler ccp_solver_;
    ie_solver::IeSovler ie_solver_;

    double best_vol_;
    double mean_ie_time_;
    double mean_ccp_time_;
    int max_iteration_ = 20;
    int iter_;

    template <typename T>
    std::vector<size_t> arg_sort(const std::vector<T> &vec) {
        std::vector<size_t> idx(vec.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&vec](size_t i0, size_t i1) {return vec[i0] < vec[i1];});
        return idx;
    }
};

}