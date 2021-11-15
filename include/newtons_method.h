#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g,
Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H)
{
    for (unsigned int i = 0; i < maxSteps; ++i)
    {
        // Calculate gradient and hessian
        // T and V are not blowing up, only H and G are!
        H(tmp_H, x0);
        g(tmp_g, x0);

        std::cout << "g: " << tmp_g.squaredNorm() << std::endl;
        std::cout <<"h: " << tmp_H.squaredNorm() << std::endl;

        if (tmp_g.squaredNorm() < std::numeric_limits<double>::epsilon()){
            std::cout <<"return tmp_g"<< std::endl;
            return 0.0;

        }
        // LDLT soler for Hd = -g
        Eigen::VectorXd d;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        ldlt.compute(tmp_H);
        d = -1 * ldlt.solve(tmp_g);


      double alpha = 1;
      double p = 0.5;
      double c = 1e-8;
      double energy_before = f(x0);

      while(true)
      {
        double energy_now = f(x0 + alpha * d);
         if ( energy_now <= energy_before + c * d.transpose() * tmp_g){
//            std::cout <<"return alpha"<< std::endl;
            break;
         }
         alpha *= p;
         if (alpha < std::numeric_limits<double>::epsilon()){
//            std::cout <<"return alpha too small"<< std::endl;
            return 0.0;
         }
      }

      x0 += alpha * d;
//      std::cout <<"energy_before: " << energy_before << std::endl;
//      std::cout <<"energy_after: " << f(x0) << std::endl;
   }
   return 0.0;
}
