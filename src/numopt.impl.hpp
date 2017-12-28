#ifndef __NUMOPT_IMPL_HPP__
#define __NUMOPT_IMPL_HPP__

#include <algorithm>

arma::vec numopt::qpez(
  arma::mat const & G, arma::vec const & c, arma::mat const & A)
{
  auto const tmp1 = arma::solve(G, A.t());
  auto const tmp2 = arma::solve(G, c);

  arma::mat L, U, P;
  arma::lu(L, U, P, A*tmp1);

  auto const mu = arma::solve(U, arma::solve(L, P*A*tmp2));

  return tmp1*mu - tmp2;
}

arma::vec numopt::qpi(
  arma::mat const & G, arma::vec const & c, arma::mat const & A,
  arma::vec const & b, double tol, int niters)
{
  int m = A.n_rows;
  int n = G.n_rows;

  arma::mat A_active;
  arma::vec x(n), y(n), p(n), lam(m), numer(m), denom(m), mask(m), masked(m);
  arma::uvec W = arma::find(arma::abs(A*x - b) <= tol), ind, V;

  auto const complement = [m] (arma::uvec const & W) {
    arma::uvec V {m - W.n_rows};
    int i = 0, j = 0;
    for (int k: W) {
      while (i < k) {
        V[j++] = i++;
      }
      ++i;
    }
    while (i < m) {
      V[j++] = i++;
    }
    return V;
  };

  double alpha;
  int k = 0;
  while (true) {
    A_active = A.rows(W);
    y = G*x + c;
    p = qpez(G, y, A_active);

    if (arma::norm(p, "inf") <= tol) {
      lam.fill(arma::datum::nan);
      lam(W) = arma::solve(A_active.t(), y);
      if (all(lam(W) >= 0)) {
        break;
      }
      ind = arma::find(W == lam.index_min(), 1);
      W.shed_row(ind(0));
    } else {
      V = complement(W);
      alpha = 1;
      denom.fill(arma::datum::nan);
      denom(V) = A.rows(V)*p;
      ind = denom < 0;
      if (any(ind)) {
        numer.fill(arma::datum::nan);
        numer(V) = b(V) - A.rows(V)*x;
        mask.fill(arma::fill::zeros);
        mask(1 - ind).fill(arma::datum::nan);
        masked = mask + numer/denom;
        alpha = std::max(0., std::min(masked.min(), alpha));
        if (alpha < 1) {
          arma::uvec tmp = {masked.index_min()};
          W = arma::sort(arma::join_vert(W, tmp));
        }
      }
      x += alpha*p;
    }

    ++k;
  }
  if (k == niters) {
    // TODO: handle error
  }
  return x;
}

#endif // __NUMOPT_IMPL_HPP__
