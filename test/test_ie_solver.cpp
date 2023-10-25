#include <iris_2d/ie_solver.h>

int main()
{
    iris_2d::PlaneMatrix A(3,2);
    A << -1, 0,
        0, -1,
        1, 1;
    iris_2d::PlaneVector b(3);
    b << 0, 0, 1;
    iris_2d::Matrix C;
    C << 0.1, 0.0,
        0.0, 0.1;
    iris_2d::Vector d;
    d.fill(b.mean());
    
    iris_2d::ie_solver::IeSovler solver;
    solver.solve(A, b, C, d);
    std::cout << "finish solver.solve()" << std::endl;

    std::cout << "C: " << solver.getC() <<  std::endl;
    std::cout << "d: " << solver.getD() << std::endl;

    return 0;
}