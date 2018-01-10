#include "numopt.hpp"

#include <algorithm>
#include <cassert>

arma::uvec numopt::setdiff(unsigned low, unsigned high, arma::uvec const & u) {
  int n = high - low - u.n_rows;
  assert(n >= 0);
  if (u.n_rows == 0) {
    arma::uvec v(n);
    for (unsigned i = 0, j = low; j < high; ++i, ++j) {
      v(i) = j;
    }
    return v;
  } else {
    assert(u.is_sorted());
    assert(low <= u.min());
    assert(u.max() < high);
    arma::uvec v(n);
    unsigned j = low, k = 0;
    for (unsigned i = 0; i < u.n_rows; ++i) {
      while (j < u(i)) {
        v(k++) = j++;
      }
      ++j;
    }
    while (j < high) {
      v(k++) = j++;
    }
    return v;
  }
}

arma::vec numopt::qpez_schur(
  arma::mat const & G, arma::vec const & c, arma::mat const & A)
{
  using arma::solve;
  auto const tmp1 = solve(G, A.t());
  auto const tmp2 = solve(G, c);
  return tmp1*solve(A*tmp1, A*tmp2) - tmp2;
}

template <>
void
numopt::qpe_baryplex<2, 0>(double const * G, double const * c, double * x)
{
  x[0] = 0.0;
  x[1] = -c[1]/G[3];
}

template <>
void
numopt::qpe_baryplex<2, 1>(double const * G, double const * c, double * x)
{
  x[0] = -c[0]/G[0];
  x[1] = 0.0;
}

template <>
void
numopt::qpe_baryplex<2, 2>(double const * G, double const * c, double * x)
{
  double pZ = ((G[3] - G[0])/2 + c[1] - c[0])/(G[0] - G[1] - G[2] + G[3]);
  x[0] = 0.5 + pZ;
  x[1] = 0.5 - pZ;
}

arma::vec numopt::qpi(
  arma::mat const & G, arma::vec const & c, arma::mat const & A,
  arma::vec const & b, arma::vec const * x0, bool * error, double tol,
  int niters)
{
  if (error) {
    *error = false;
  }

  int m = A.n_rows;
  int n = G.n_rows;

  arma::mat A_active;
  arma::vec x(n), y(n), p(n), lam(m), numer(m), denom(m), mask(m), masked(m);
  arma::uvec W, ind, V;

  if (x0) {
    x = *x0;
  } else {
    x.fill(arma::fill::zeros);
  }
  W = arma::find(arma::abs(A*x - b) <= tol);

  double alpha;
  int k = 0;
  while (true) {
    A_active = A.rows(W);
    y = G*x + c;
    p = qpez_schur(G, y, A_active);

    if (arma::norm(p, "inf") <= tol) {
      lam.fill(arma::datum::nan);
      lam(W) = arma::solve(A_active.t(), y);
      if (all(lam(W) >= 0)) {
        break;
      }
      ind = arma::find(W == lam.index_min(), 1);
      W.shed_row(ind(0));
    } else {
      V = setdiff(0, m, W);
      alpha = 1;
      denom.fill(arma::datum::nan);
      denom(V) = A.rows(V)*p;
      ind = denom < 0;
      if (any(ind)) {
        numer.fill(arma::datum::nan);
        numer(V) = b(V) - A.rows(V)*x;
        mask.fill(arma::fill::zeros);
        mask.elem(arma::find(ind == 0)).fill(arma::datum::nan);
        masked = mask + numer/denom;
        alpha = std::max(0., std::min(masked.min(), alpha));
        if (alpha < 1) {
          arma::uvec tmp = {masked.index_min()};
          W = arma::sort(arma::join_vert(W, tmp));
        }
      }
      x += alpha*p;
    }

    if (++k == niters && error) {
      *error = true;
      break;
    }
  }
  return x;
}

arma::vec numopt::sqp(
  numopt::field_t const & f, numopt::grad_t const & df,
  numopt::hess_t const & d2f, arma::mat const & A, arma::vec const & b,
  arma::vec const & xinit, bool * error, double tol, int niters)
{
  if (error) {
    *error = false;
  }

  arma::mat H;
  arma::vec x0, x1 = xinit, c, g;
  double f0, f1 = f(x1), lambda_min, qpi_tol, c1 = 1e-4, alpha;

  bool qpi_error;
  int k = 0, found_opt, qpi_niters = 10;
  while (true) {
    // Compute the Hessian for the current iterate and ensure that
    // it's symmetric
    H = d2f(x1);
    H = (H + H.t())/2;

    // Perturb the Hessian if it isn't positive definite
    lambda_min = arma::min(arma::eig_sym(H));
    if (lambda_min < 0) {
      H -= 1.1*lambda_min*arma::eye<arma::mat>(arma::size(H));
      assert(arma::min(arma::eig_sym(H)) > 0);
    }
    
    // Compute load vector for quadratic program
    c = df(x1) - H*x1;

    // Compute descent direction by solving inequality-constrained
    // quadratic program
    found_opt = false;
    qpi_tol = EPS(double);
    while (!found_opt) {
      g = qpi(H, c, A, b, &x1, &qpi_error, qpi_tol, qpi_niters) - x1;
      if (qpi_error) {
        qpi_tol *= 10;
      } else {
        found_opt = true;
      }
    }

    // Compute step size
    alpha = 1;
    if (arma::norm(g, "inf") > tol) {
      while (f(x1 + alpha*g) > f1 + c1*alpha*arma::dot(df(x1), g)) {
        alpha /= 2;
      }
    }

    x0 = x1;
    x1 += alpha*g;
    f0 = f1;
    f1 = f(x1);

    if (arma::norm(x1 - x0, "inf") <= tol || fabs(f1 - f0) <= tol) {
      break;
    }
    
    if (++k == niters && error) {
      *error = true;
      break;
    }
  }
  return x1;
}
