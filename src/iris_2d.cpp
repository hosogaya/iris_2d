#include <iris_2d/iris_2d.h>

namespace iris_2d
{
void Problem::reset()
{
    getCRef().coeffRef(0,0) = ELLIPSOID_C_EPSILON;
    getCRef().coeffRef(1,0) = 0.0;
    getCRef().coeffRef(0,1) = 0.0;
    getCRef().coeffRef(1,1) = ELLIPSOID_C_EPSILON;

    best_vol_ = std::pow(ELLIPSOID_C_EPSILON, 2.0);
    mean_ie_time_ = 0;
    mean_ccp_time_ = 0;
    iter_ = 0;
}

void Problem::initialize(const Vector& seed, const std::vector<Obstacle>& obs)
{
    getEllipsoidRef().setD(seed);
    obs_ = obs;

    reset();
}

bool Problem::solve()
{
    double volume = 0.0;
    int i=0;
    for (; i<max_iteration_; ++i)
    {
        if (!separatingPlanes()) return false;
        if (!inflateRegion(volume)) return false;
        if (std::abs(volume - best_vol_)/best_vol_ < 2e-2) break;
        best_vol_ = volume;
        ++iter_;
    }
    if (i == max_iteration_) return false;

    return true;
}

HyperPlane Problem::tangentPlane(const Vector& point, const Matrix& Cinv2)
{
    Vector normal = ((Cinv2 + Cinv2.transpose())*(point-getD())).normalized();
    return std::make_pair(normal, normal.transpose()*point);
}

bool Problem::separatingPlanes()
{
    auto start = std::chrono::system_clock::now();
    Matrix Cinv = getC().inverse();
    Matrix Cinv2 = Cinv*Cinv.transpose();
    Polyhedron poly;

    std::vector<Obstacle> img_obs(obs_.size());
    for (size_t i=0; i<obs_.size(); ++i)
    {
        img_obs[i] = Cinv*(obs_[i].colwise() - getD()); 
    }

    // std::vector<Vector> image_squared_dists(obs_.size());
    // for (size_t i=0; i<obs_.size(); ++i) image_squared_dists[i] = img_obs[i].colwise().squaredNorm();

    // std::vector<double> min_squared_dists(obs_.size());
    // for (size_t i=0; i<obs_.size(); ++i) min_squared_dists[i] = image_squared_dists[i].minCoeff();

    // std::vector<size_t> obs_sort_idx = arg_sort(min_squared_dists);

    std::vector<HyperPlane> planes(0);

    for (size_t i=0; i<obs_.size(); ++i)
    {
        for (size_t j=0; j<planes.size(); ++j)
        {
            if (((planes[j].first.transpose()*obs_[i]).array() - planes[j].second > 0.0).all()) continue;
        }
        Vector y_star;
        if (obs_[i].cols() == 2) 
        {
            // const Vector d = getD();
            const Vector p0 = img_obs[i].col(0);
            const Vector p1 = img_obs[i].col(1);
            double alpha = (p1-p0).dot(-p0) / (p1-p0).norm();
            if (alpha > 1.0) alpha = 1.0;
            else if (alpha < 0.0) alpha = 0.0;
            y_star = alpha* (p1-p0) + p0;
        }
        else 
        {
            if (!ccp_solver_.solve(img_obs[i])) return false;
            y_star = ccp_solver_.getX();
        }
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
    auto end = std::chrono::system_clock::now();
    double elasped = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    mean_ccp_time_ = (iter_*mean_ccp_time_ + elasped)/(iter_+1);
    return true;
}

bool Problem::inflateRegion(double& volume)
{
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    if (!ie_solver_.solve(getA(), getB(), getC(), getD())) return false;
    getEllipsoidRef().setC(ie_solver_.getC());
    getEllipsoidRef().setD(ie_solver_.getD());

    volume = getEllipsoidRef().getVolume();

    // std::cout << "C: " << getC() << std::endl;
    // std::cout << "d: " << getD() << std::endl;

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    double elasped = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    mean_ie_time_ = (iter_*mean_ie_time_ + elasped)/(iter_+1);

    return true;
}
}