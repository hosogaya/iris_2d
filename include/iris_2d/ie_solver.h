// Optimization problem for inscribed ellipsoid (10)

#pragma once

#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <iris_2d/geometry.h>
#include <memory>

namespace iris_2d {

namespace ie_solver {
// vars_num
// eigen value = dim
// eigen vector = dim*dim
// d = dim
class Variables : public ifopt::VariableSet
{
public:
    Variables(const int dim) : VariableSet(1+dim + dim, "ie_vars"), 
		dim_(dim)
    {
		x_.resize(GetRows());
		x_.setZero();
		bounds_.resize(GetRows());
		for (int i=0; i<dim_; ++i)
		{
			bounds_[i] = ifopt::Bounds(0.0, ifopt::inf);
		}
		bounds_[dim_] = ifopt::Bounds(0.0, M_PI);
		for (size_t i=dim_+1; i<bounds_.size(); ++i)
		{
			bounds_[i] = ifopt::Bounds(-ifopt::inf, ifopt::inf);
		}
    }

    void SetVariables(const VectorXd& x) override
    {
		assert(x_.size() == x.size());
		x_ = x;
    }

    VectorXd GetValues() const override
    {
		return x_;
    }

    VecBound GetBounds() const override
    {
		return bounds_;
    }

    void SetBounds(const VecBound& b)
    {
		assert(b.size() == bounds_.size());
		bounds_ = b;
    }

private:
    Eigen::VectorXd x_;
    VecBound bounds_;
	int dim_;
}; // class Variables


inline autodiff::VectorXreal f_constraints(const autodiff::Vector2real& eigen_value, const autodiff::real& rad, const autodiff::Vector2real& d, const PlaneMatrix& A, const PlaneVector& b)
{
	// std::cout << "f_constraints called" << std::endl;
	int plane_num = b.size();
	autodiff::VectorXreal cons(plane_num+1);
	autodiff::Matrix2real c_mat;
	autodiff::Matrix2real eigen_value_diagonal;
	eigen_value_diagonal << eigen_value[0], 0.0,
							0.0, eigen_value[1];
	// create C matrix
	autodiff::Matrix2real eigen_vector_matrix;
	eigen_vector_matrix << cos(rad), -sin(rad), 
							sin(rad), cos(rad);
	c_mat = eigen_vector_matrix*eigen_value_diagonal*(eigen_vector_matrix.transpose());

	// |aC| + ad - b < 0 Eq.(10)
	for (int i=0; i<plane_num; ++i)
	{
		cons[i] =  (A.row(i)*c_mat).norm() + A.row(i).dot(d) - autodiff::real1st(b[i]);
	}
	cons[plane_num] = eigen_value[1] - eigen_value[0];

	// std::cout << "f_constraints called" << std::endl;

	return cons;
}

class Constraints : public ifopt::ConstraintSet
{
public:
  	Constraints(int plane_num, const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
	 : ConstraintSet(plane_num+1, "ie_cons"), plane_num_(plane_num), A_(A), b_(b)
	{
		// constrains < 0
		bounds_.resize(GetRows());
		for (int i=0; i<GetRows(); ++i)
			bounds_[i] = ifopt::Bounds(-ifopt::inf, 0.0);
	}

	VectorXd GetValues() const override
	{
		// std::cout << "GetConstraints called" << std::endl;
		Eigen::VectorXd x = GetVariables()->GetComponent("ie_vars")->GetValues();
		autodiff::Vector2real eigen_value = x.topRows(2);
		autodiff::real rad = x(2);
		autodiff::Vector2real d = x.bottomRows(2);

		Eigen::VectorXd cons(GetRows());
		autodiff::VectorXreal cons_autodiff = f_constraints(eigen_value, rad, d, A_, b_);
		for (int i=0; i<GetRows(); ++i) cons[i] = cons_autodiff(i).val();
		// std::cout << "GetConstraints finished" << std::endl;
		return cons;
	};

	VecBound GetBounds() const override
	{
		return bounds_;
	}

	void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override
	{
		// std::cout << "FillJacobian for constraints called" << std::endl;
		if (var_set == "ie_vars")
		{
			Eigen::VectorXd x = GetVariables()->GetComponent("ie_vars")->GetValues();
			autodiff::Vector2real eigen_value = x.topRows(2);
			autodiff::real rad = x(2);
			autodiff::Vector2real d = x.bottomRows(2);

			// size of jac: cons*vars_num
			jac_block.resize(GetRows(), x.size());

			Eigen::MatrixXd jac = autodiff::jacobian(f_constraints, autodiff::wrt(eigen_value, rad, d), autodiff::at(eigen_value, rad, d, A_, b_));

			for (int i=0; i<GetRows(); ++i)
			{
				for (int j=0; j<x.size(); ++j)
				{
					jac_block.coeffRef(i,j) = jac.coeffRef(i, j);
				}
			}
		}
		// std::cout << "FillJacobian for constraints finished" << std::endl;
	}
protected:
	const int plane_num_;
	VecBound bounds_;
	const Eigen::MatrixXd A_;
	const Eigen::VectorXd b_;
};

// cost function
// inline autodiff::real f_cost(const autodiff::Vector2real& eigen_value) 
// {
// 	// std::cout << "f const called" << std::endl;
// 	// std::cout << "f const finished" << std::endl;

// 	return -log(eigen_value[0]*eigen_value[1]);
// }

class Cost : public ifopt::CostTerm
{
public:
	Cost() : CostTerm("ie_cost") {}


	double GetCost() const override
	{
		// std::cout << "GetCost called" << std::endl;
		Eigen::VectorXd x = GetVariables()->GetComponent("ie_vars")->GetValues();
		// autodiff::Vector2real eigen_value = x.topRows(2);
		// autodiff::real rad = x(2);
		// autodiff::Vector2real d = x.bottomRows(2);
		
		// std::cout << "GetCost finished" << std::endl;
		return -std::log(x(0)*x(1));
	};

	void FillJacobianBlock(std::string var_set, Jacobian& jac) const override
	{
		// std::cout << "FillJacobian for cost called" << std::endl;
		if (var_set == "ie_vars")
		{
			Eigen::VectorXd x = GetVariables()->GetComponent("ie_vars")->GetValues();
			// autodiff::Vector2real eigen_value = x.topRows(2);
			// autodiff::real rad = x(2);
			// autodiff::Vector2real d = x.bottomRows(2);

			// Eigen::VectorXd grad = autodiff::gradient(f_cost, autodiff::wrt(eigen_value), autodiff::at(eigen_value));
		
			jac.resize(1, x.size());
			jac.toDense().fill(0.0);

			jac.coeffRef(0,0) =-1/x(0);
			jac.coeffRef(0,1) =-1/x(1);
		}
		// std::cout << "FillJacobian for cost finished" << std::endl;
	}
private:
};

class IeSovler
{
public:
	IeSovler() 
	{
		ipopt_.SetOption("tol", 1.0e-4); // tolerance
		ipopt_.SetOption("max_iter", 50); // max iteration
		ipopt_.SetOption("max_wall_time", 1); // maximum computation time
		ipopt_.SetOption("constr_viol_tol", 0.001); // constraint violation tolerance
		ipopt_.SetOption("print_level", 4); // supress log

		ipopt_.SetOption("dependency_detection_with_rhs", "no"); // check constrains (yes) or only check gradients of constraints (no).
		ipopt_.SetOption("jacobian_approximation", "exact"); // how to compute jacobian
		ipopt_.SetOption("gradient_approximation", "exact"); // how to compute gradient
		ipopt_.SetOption("diverging_iterates_tol", 10.0e20); // Threshold of the number of primal iterations
		ipopt_.SetOption("nlp_scaling_method", "gradient-based"); // how to solve nlp
		ipopt_.SetOption("nlp_scaling_max_gradient", 1000.0); // how to solve nlp
		ipopt_.SetOption("sb", "yes");

	}
	~IeSovler() {}

	// set A,b, C,d
	bool solve(const PlaneMatrix& A, const PlaneVector& b, const Matrix& C, const Vector& d)
	{	
		dim_ = d.size();
		plane_num_ = b.size();
		// std::cout << "dim: " << dim_ << std::endl;
		// std::cout << "plane_num: " << plane_num_ << std::endl;
		Eigen::SelfAdjointEigenSolver<Matrix> solver(C);
		if (solver.info() != Eigen::Success) return false;

		Vector eigen_value = solver.eigenvalues();
		Matrix eigen_vector = solver.eigenvectors();
		double rad = std::acos(eigen_vector.col(1).dot(Vector::UnitX()));
		int frist, second;
		if (rad > M_PI_2) {
			frist = 0;
			second = 1;
			rad = 3*M_PI_2 - rad;
		}
		else {
			frist = 1;
			second = 0;
		}

		// eigen value, rad, d
		Eigen::VectorXd x(dim_+1+dim_);
		x(0) = eigen_value(frist);
		x(1) = eigen_value(second);
		x(2) = eigen_value(rad);
		for (int i=0; i<dim_; ++i) x(dim_+1 + i) = d(i);
		vars_ = std::make_shared<Variables>(dim_);
		vars_->SetVariables(x);
		// std::cout << "init vars: " << vars_->GetValues() << std::endl;
		// std::cout << "bounds num; " << vars_->GetBounds().size() << std::endl;

		cons_ = std::make_shared<Constraints>(plane_num_, A, b);
		// std::cout << "A: " << A << std::endl;
		// std::cout << "b: " << b << std::endl;
		// std::cout << "cons num: " << cons_->GetRows() << std::endl;

		cost_ = std::make_shared<Cost>();

		ifopt::Problem prob_;
		prob_.AddVariableSet(vars_);
		prob_.AddConstraintSet(cons_);
		prob_.AddCostSet(cost_);

		ipopt_.Solve(prob_);
		if (ipopt_.GetReturnStatus() != 0) return false;

		setOptValues(prob_.GetOptVariables()->GetValues());
		return true;
	}
	const Matrix& getC() const {return C_;}
	const Vector& getD() const {return d_;}

private:
	ifopt::IpoptSolver ipopt_;
	std::shared_ptr<Variables> vars_;
	std::shared_ptr<Constraints> cons_;
	std::shared_ptr<Cost> cost_;
	int dim_;
	int plane_num_;
	Matrix C_;
	Vector d_;

	void setOptValues(const Eigen::VectorXd& x)
	{
		// std::cout << "setOptValues called. size of x: " << x.size() << std::endl;
		Matrix eigen_value_diagonal;
		Matrix eigen_vector_matrix;
		for (int i=0; i<dim_; ++i)
			eigen_value_diagonal.coeffRef(i, i) = x(i);

		// std::cout << "eigen value diagonal: " << std::endl << eigen_value_diagonal << std::endl;
		
		double rad = x(dim_);
		eigen_vector_matrix << std::cos(rad), -std::sin(rad),
								std::sin(rad), std::cos(rad);

		// std::cout << "eigen vector matrix: " << std::endl << eigen_vector_matrix << std::endl;

		C_ = eigen_vector_matrix*eigen_value_diagonal*eigen_vector_matrix.transpose();

		d_ = x.bottomRows(dim_);

		// std::cout << "setOptValues finished" << std::endl;
	}
};
}

}