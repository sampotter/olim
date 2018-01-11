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

#define __compute_p() do {                          \
    double det = G[0]*G[3] - G[1]*G[2];             \
    p[0] = -(x[0] + (G[3]*c[0] - G[1]*c[1])/det);   \
    p[1] = -(x[1] + (G[0]*c[1] - G[2]*c[0])/det);   \
  } while (0)


#define __compute_y() do {                      \
    y[0] = G[0]*x[0] + G[2]*x[1] + c[0];        \
    y[1] = G[1]*x[0] + G[3]*x[1] + c[1];        \
  } while (0)

#define __num_active()                                                  \
  ((active[0] ? 1 : 0) + (active[1] ? 1 : 0) + (active[2] ? 1 : 0))

template <>
void
numopt::qpi_baryplex<2>(double const * G, double const * c, double const * x0,
                        double * x, bool * error, double tol, int niters)
{
  using std::max;
  using std::min;

  assert(x != nullptr);

  if (error) {
    *error = false;
  }

  double xprev[2], p[2], y[2], alpha, alpha_new;
  if (x0) {
    x[0] = x0[0];
    x[1] = x0[1];
  } else {
    x[0] = x[1] = 0.0;
  }

  bool active[3], ind[3];
  active[0] = fabs(x[0]) <= tol;
  active[1] = fabs(x[1]) <= tol;
  active[2] = fabs(1 - x[0] - x[1]) <= tol;

  int k = 0, num_active, argmin;
  while (true) {
    xprev[0] = x[0];
    xprev[1] = x[1];
    num_active = __num_active();

    __compute_y();

    if (num_active == 0) {
      __compute_p();
    } else if (num_active == 1) {
      // TODO: fix this so that we just immediately compute p
      double x_[2];
      if (active[0]) {
        qpe_baryplex<2, 0>(G, c, x_);
      } else if (active[1]) {
        qpe_baryplex<2, 1>(G, c, x_);
      } else if (active[2]) {
        qpe_baryplex<2, 2>(G, c, x_);
      } else {
        assert(false);
      }
      p[0] = x_[0] - xprev[0];
      p[1] = x_[1] - xprev[1];
    } else if (num_active == 2) {
      p[0] = p[1] = 0.0;
    }

    if (max(fabs(p[0]), fabs(p[1])) <= tol) {
      if (num_active == 0) {
        break;
      } else if (num_active == 1) {
        if (active[0]) {
          if (y[0] >= 0) break;
          else active[0] = false;
        } else if (active[1]) {
          if (y[1] >= 0) break;
          else active[1] = false;
        } else if (active[2]) {
          if (y[0] + y[1] <= 0) break;
          else active[2] = false;
        } else {
          assert(false);
        }
      } else if (num_active == 2) {
        if (active[0] && active[1]) {
          if (y[0] >= 0 && y[1] >= 0) break;
          else active[y[0] < y[1] ? 0 : 1] = false;
        } else if (active[0] && active[2]) {
          if (y[1] <= 0 && y[1] <= y[0]) break;
          else active[y[0] < 0 ? 0 : 2] = false;
        } else if (active[1] && active[2]) {
          if (y[0] <= 0 && y[0] <= y[1]) break;
          else active[y[1] < 0 ? 1 : 2] = false;
        } else {
          assert(false);
        }
      } else {
        assert(false);
      }
    } else {
      alpha = 1;
      argmin = -1;
      ind[0] = !active[0] && p[0] < 0;
      ind[1] = !active[1] && p[1] < 0;
      ind[2] = !active[2] && p[0] + p[1] > 0;
      if (ind[0] || ind[1] || ind[2]) {
        if (ind[0]) {
          alpha_new = -x[0]/p[0];
          if (alpha_new < alpha) {
            alpha = alpha_new;
            argmin = 0;
          }
        }
        if (ind[1]) {
          alpha_new = -x[1]/p[1];
          if (alpha_new < alpha) {
            alpha = alpha_new;
            argmin = 1;
          }
        }
        if (ind[2]) {
          alpha_new = (1 - x[0] - x[1])/(p[0] + p[1]);
          if (alpha_new < alpha) {
            alpha = alpha_new;
            argmin = 2;
          }
        }
        alpha = max(0.0, min(alpha, 1.0));
        if (alpha < 1) {
          active[argmin] = true;
        }
      }
      x[0] += alpha*p[0];
      x[1] += alpha*p[1];
    }

    if (++k == niters) {
      if (error) *error = true;
      break;
    }
  }
}

#undef __compute_p
#undef __compute_y
#undef __num_active

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

#define __compute_lambda_min() do {                                 \
    double half_tr = (G[0] + G[3])/2, det = G[0]*G[3] - G[1]*G[2];  \
    lambda_min = half_tr - sqrt(half_tr*half_tr - det);             \
  } while (0)                                                       \

template <>
void
numopt::sqp_baryplex<2>(cost_func<2> * func, double * x, bool * error,
                        double tol, int niters)
{
  using std::max;

  if (error) *error = false;

  double G[4], x0[2], x1[2] = {1./3, 1./3}, c[2], g[2], f0, f1,
    lambda_min, qpi_tol, c1 = 1e-4, alpha;
  bool qpi_error, found_opt;
  int k = 0, qpi_niters = 10;

  func->set_lambda(x1);
  func->eval(f1);

  while (true) {
    // Compute Hessian and perturb it if it isn't positive definite
    func->hess((double (*)[2]) &G);
    __compute_lambda_min();
    if (lambda_min < 0) {
      G[0] -= 1.1*lambda_min;
      G[3] -= 1.1*lambda_min;
    }

    // Compute load vector for quadratic program
    func->grad(c);
    c[0] -= G[0]*x1[0] + G[1]*x1[1];
    c[1] -= G[2]*x1[0] + G[3]*x1[1];

    // Compute descent direction by solving inequality-constrained
    // quadratic program
    found_opt = false;
    qpi_tol = tol;
    while (!found_opt) {
      qpi_baryplex<2>(G, c, x1, g, &qpi_error, qpi_tol, qpi_niters);
      if (qpi_error) qpi_tol *= 10;
      else found_opt = true;
    }
    g[0] -= x1[0];
    g[1] -= x1[1];

    // Compute step size... This is a bit of a convoluted dance, but
    // is efficient
    alpha = 1;
    if (max(fabs(g[0]), fabs(g[1])) > tol) {
      // TODO: use x0 instead of tmp to save space
      double tmp[2], lhs, rhs;
recompute:
      func->grad(tmp);
      rhs = f1 + c1*(tmp[0]*g[0] + tmp[1]*g[1]);
      tmp[0] = x1[0] + alpha*g[0];
      tmp[1] = x1[1] + alpha*g[1];
      func->set_lambda(tmp);
      func->eval(lhs);
      if (lhs > rhs) {
        alpha /= 2;
        func->set_lambda(x1);
        goto recompute;
      }
    }

    // Save current values for next iteration
    x0[0] = x1[0];
    x0[1] = x1[1];
    x1[0] += alpha*g[0];
    x1[1] += alpha*g[1];
    f0 = f1;
    func->set_lambda(x1);
    func->eval(f1);

    // Check for convergence
    if (max(fabs(x1[0] - x0[0]), fabs(x1[1] - x0[1])) <= tol ||
        fabs(f1 - f0) <= tol) {
      break;
    }

    // Check if we've reached our max number of iterations
    if (++k == niters) {
      if (error) *error = true;
      break;
    }
  }

  x[0] = x1[0];
  x[1] = x1[1];
}

#undef __compute_lambda_min
