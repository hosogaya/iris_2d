#include <iris_2d/iris_2d.h>

namespace iris_2d
{

bool Problem::initialize(const Vector& seed, const std::vector<Obstacle>& obs)
{
    getEllipsoidRef().setD(seed);
    getCRef().coeffRef(0,0) = ELLIPSOID_C_EPSILON;
    getCRef().coeffRef(1,0) = 0.0;
    getCRef().coeffRef(0,1) = 0.0;
    getCRef().coeffRef(1,1) = ELLIPSOID_C_EPSILON;

    obs_ = obs;

    best_vol_ = std::pow(ELLIPSOID_C_EPSILON, 2.0);

    return true;
}

bool Problem::solve()
{
    double volume = 0.0;
    while (true)
    {
        if (!separatingPlanes()) return false;
        if (!inflateRegion(volume)) return false;
        if (std::abs(volume - best_vol_)/best_vol_ < 2e-2) break;
        best_vol_ = volume;
    }

    return true;
}

HyperPlane Problem::tangentPlane(const Vector& point, const Matrix& Cinv2)
{
    Vector normal = ((Cinv2 + Cinv2.transpose())*(point-getD())).normalized();
    return std::make_pair(normal, normal.transpose()*point);
}

bool Problem::separatingPlanes()
{
    Matrix Cinv = getC().inverse();
    Matrix Cinv2 = Cinv*Cinv.transpose();
    Polyhedron poly;

    std::vector<Obstacle> img_obs(obs_.size());
    for (size_t i=0; i<obs_.size(); ++i)
    {
        img_obs[i] = Cinv*(obs_[i].colwise() - getD()); 
    }

    std::vector<HyperPlane> planes;

    for (size_t i=0; i<obs_.size(); ++i)
    {
        if (!ccp_solver_.solve(img_obs[i])) return false;

        Vector y_star = ccp_solver_.getX();
        if (y_star.squaredNorm() < 1e-6)
        {
            return false;
        }
        else 
        {
            Vector x_star = getC()*y_star + getD();
            // std::cout << "x_star: " << x_star << std::endl;
            planes.emplace_back(tangentPlane(x_star, Cinv2));
        }
    }

    PlaneMatrix A_mat(planes.size(), 2);
    PlaneVector b_vec(planes.size());
    for (size_t i=0; i<planes.size(); ++i)
    {
        A_mat.row(i) = planes[i].first.transpose();
        b_vec[i] = planes[i].second;
    }

    // std::cout << "A: " << A_mat << std::endl;
    // std::cout << "b: " << b_vec << std::endl;

    getPolyhedronRef().setA(A_mat);
    getPolyhedronRef().setB(b_vec);

    return true;
}

bool Problem::inflateRegion(double& volume)
{
    if (!ie_solver_.solve(getA(), getB(), getC(), getD())) return false;
    getEllipsoidRef().setC(ie_solver_.getC());
    getEllipsoidRef().setD(ie_solver_.getD());

    volume = getEllipsoidRef().getVolume();

    // std::cout << "C: " << getC() << std::endl;
    // std::cout << "d: " << getD() << std::endl;

    return true;
}
}