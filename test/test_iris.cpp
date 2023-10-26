#include <iris_2d/iris_2d.h>

int main()
{
    iris_2d::Problem problem;
    std::vector<iris_2d::Obstacle> obstacles;
    iris_2d::Obstacle obs(2,2);
    obs << 0, 1,
         0, 0;
    obstacles.push_back(obs);
    // obs << 1, 1,
    //         0, 1;
    // obstacles.push_back(obs);
    obs << 1, 0,
            0, 1;
    obstacles.push_back(obs);
    obs << 0, 0,
            1, 0; 
    obstacles.push_back(obs);

    iris_2d::Vector seed;
    seed << 0.1, 0.1;

    problem.initialize(seed, obstacles);
    problem.solve();

    std::cout << "C: " << problem.getC() << std::endl;
    std::cout << "d: " << problem.getD() << std::endl;
    std::cout << "A: " << problem.getA() << std::endl;
    std::cout << "b: " << problem.getB() << std::endl;
    std::cout << "ccpTime: " << problem.getCcpTime() << std::endl;
    std::cout << "ieTime: " << problem.getIeTime() << std::endl;
    std::cout << "Iteratoin: " << problem.getIteration() << std::endl;
    
    return 0;
}