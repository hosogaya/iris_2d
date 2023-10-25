#include <iris_2d/ccp_solver.h>
#include <iostream>

int main() 
{
    iris_2d::ccp_solver::CcpSovler solver;
    iris_2d::Obstacle v(2, 4);
    v << 0, 1, 1, 0,
        0, 0, 1, 1;
    solver.solve(v);

    std::cout << "x: " << solver.getX() << std::endl;

    return 1;
}