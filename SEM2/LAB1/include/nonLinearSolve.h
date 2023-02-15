#ifndef LAB1_NONLINEARSOLVE_H
#define LAB1_NONLINEARSOLVE_H

//std::vector<TT> nonLinearSolve(const std::function<vector<TT>(const std::vector<TT>&)>&,
//                 const std::vector<TT> &initPoints, TT eps, size_t &iterCount);

std::vector<TT> Newton(const std::string &method, size_t dim, const std::vector<TT> &y0);

std::vector<TT> f(const std::vector<TT>& x);

#endif // LAB1_NONLINEARSOLVE_H
