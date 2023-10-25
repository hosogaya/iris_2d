// closet point on an obstacle 
// optimization problem Eq.(4)
#pragma once

#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <iris_2d/geometry.h>

namespace iris_2d
{
namespace ccp_solver {
// vars_num: dim + vertex_num
// size of x: dim
// size of w: vertex_num
// 0.0 <= w <= 1.0
class Variables : public ifopt::VariableSet
{
public:
    Variables(const int vertex_num) : VariableSet(vertex_num, "ccp_vars"), 
		vertex_num_(vertex_num)
    {
		x_.resize(GetRows());
		x_.setZero();
		bounds_.resize(GetRows());
		// bounds for w
		for (int i=0; i<vertex_num_; ++i)
		{
			bounds_[i] = ifopt::Bounds(0.0, 1.0);
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
	int vertex_num_;
}; // class Variables



// vars num; dim + vertex_num
// cons num: dim + 1
// [vi1, vi2, ... , vim]w = x (dim)
// sum(w) = 1

class Constraints : public ifopt::ConstraintSet
{
public:
  	Constraints(const iris_2d::Obstacle& obs)
	 : ConstraintSet(1, "ccp_cons"), obs_(obs)
	{
		// constrains = 0;
		bounds_.resize(GetRows());
		for (size_t i=0; i<bounds_.size(); ++i)
		{
			bounds_[i] = ifopt::Bounds(0.0, 0.0);
		}
	}

	VectorXd GetValues() const override
	{
		Eigen::VectorXd x = GetVariables()->GetComponent("ccp_vars")->GetValues();

		Eigen::VectorXd cons(GetRows());
		cons[0] = 1.0 - x.sum();
		return cons;
	};

	VecBound GetBounds() const override
	{
		return bounds_;
	}

	void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override
	{
		if (var_set == "ccp_vars")
		{
			Eigen::VectorXd x = GetVariables()->GetComponent("ccp_vars")->GetValues();

			jac_block.resize(GetRows(), x.size());
			jac_block.setZero();
			for (int j=0; j<x.size(); ++j) jac_block.coeffRef(0, j) = -1.0;
		}
	}
private:
	iris_2d::Obstacle obs_;
	VecBound bounds_;
	const autodiff::real p_ = 1.0;
};


class Cost : public ifopt::CostTerm
{
public:
	Cost(const iris_2d::Obstacle& obs) : CostTerm("ccp_cost"), obs2_(obs.transpose()*obs)
	{}


	double GetCost() const override
	{
		Eigen::VectorXd x = GetVariables()->GetComponent("ccp_vars")->GetValues();

		return x.transpose()*obs2_*x;
	};

	void FillJacobianBlock(std::string var_set, Jacobian& jac) const override
	{
		if (var_set == "ccp_vars")
		{
			Eigen::VectorXd x = GetVariables()->GetComponent("ccp_vars")->GetValues();

			jac.resize(1, x.size());
			jac.setZero();
			Eigen::VectorXd grad = (obs2_ + obs2_.transpose())*x;
			for (int i=0; i<x.size(); ++i) jac.coeffRef(0, i) = grad[i];
		}
	}
private:
	Eigen::MatrixXd obs2_;
	int dim_;
	int vertex_num_;
};

class CcpSovler
{
public:
	CcpSovler() 
	{
		ipopt_.SetOption("tol", 1.0e-8); // tolerance
		ipopt_.SetOption("max_iter", int(1e5)); // max iteration
		// ipopt_.SetOption("max_wall_time", 1e20); // maximum computation time
		ipopt_.SetOption("constr_viol_tol", 0.001); // constraint violation tolerance
		ipopt_.SetOption("print_level", 1); // supress log

		ipopt_.SetOption("dependency_detection_with_rhs", "no"); // check constrains (yes) or only check gradients of constraints (no).
		ipopt_.SetOption("jacobian_approximation", "exact"); // how to compute jacobian
		// ipopt_.SetOption("gradient_approximation", "exact"); // how to compute gradient
		ipopt_.SetOption("diverging_iterates_tol", 100.0);
		ipopt_.SetOption("nlp_scaling_method", "gradient-based"); // how to solve nlp
		ipopt_.SetOption("sb", "yes");
	}

	bool solve(const iris_2d::Obstacle& obs)
	{	
		dim_ = obs.rows();
		vertex_num_ = obs.cols();
		
		// std::cout << "dim: " << dim_ << std::endl;
		// std::cout << "vertex_num: " << vertex_num_ << std::endl;

		vars_ = std::make_shared<Variables>(vertex_num_);
		Eigen::VectorXd x(vertex_num_);
		x.fill(1.0/vertex_num_);
		vars_->SetVariables(x);
		// std::cout << "Set variables: " << x << std::endl;

		cons_ = std::make_shared<Constraints>(obs);
		cost_ = std::make_shared<Cost>(obs);
		// std::cout << "Set costs" << std::endl;

		ifopt::Problem prob;
		prob.AddVariableSet(vars_);
		prob.AddConstraintSet(cons_);
		prob.AddCostSet(cost_);

		ipopt_.Solve(prob);
		if (ipopt_.GetReturnStatus() != 0) return false;

		w_star_ = prob.GetOptVariables()->GetValues();
		x_star_ = obs*w_star_;
		return true;
	}

	const Vector& getX() {return x_star_;}
	const Eigen::VectorXd& getW() {return w_star_;}

private:
	ifopt::IpoptSolver ipopt_;
	std::shared_ptr<Variables> vars_;
	std::shared_ptr<Constraints> cons_;
	std::shared_ptr<Cost> cost_;
	int dim_;
	int vertex_num_;
	Vector x_star_;
	Eigen::VectorXd w_star_;
};
}
}