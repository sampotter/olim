import numpy as np

sqrt2 = np.sqrt(2)

r = lambda x, y: np.sqrt(np.power(x, 2) + np.power(y, 2))

s0 = lambda x, y: np.ones(x.shape)
f0 = r

s1 = lambda x, y: 1 - np.sin(r(x, y))
f1 = lambda x, y: np.cos(r(x, y)) + r(x, y) - 1;

s2 = r;
f2 = lambda x, y: (np.power(x, 2) + np.power(y, 2))/2

a, b, c = 2, 1, 0.5
A = np.array([[a, b], [b, c]])
alpha = np.pi/2
C = lambda x, y: np.array([np.cos(alpha*x), np.cos(alpha*y)])
S = lambda x, y: np.array([np.sin(alpha*x), np.sin(alpha*y)])
s3 = lambda x, y: np.linalg.norm(alpha*np.diag(C(x, y))@(A + A.T)@S(x, y))
f3 = lambda x, y: S(x, y)@A@S(x, y)

d, e, f = 2/3, -0.5, 1
Bsqrt = np.array([[d, e], [e, f]])
s4 = lambda x, y: np.linalg.norm(Bsqrt@[x, y])
f4 = lambda x, y: [x, y]@Bsqrt@[x, y]/2

s_qv = lambda x, y: 1/(0.5 + 0.5*x)
f_qv = lambda x, y: 2*np.arccosh(1 + 0.25*np.sqrt(x*x + y*y)*s_qv(x, y))

def get_field(g, X, Y):
    assert(X.shape == Y.shape)
    m, n = X.shape
    return np.array([[g(X[i, j], Y[i, j]) for j in range(n)] for i in range(m)])

def get_fields(u, s, X, Y):
    return get_field(u, X, Y), get_field(s, X, Y)

_slowness_funcs = [s0, s1, s2, s3, s4]

def slowness_funcs():
    return _slowness_funcs

_soln_funcs = [f0, f1, f2, f3, f4]

for func in _slowness_funcs + _soln_funcs:
    func.dim = 2

def slowness_funcs():
    return _slowness_funcs

_slowness_func_names = {
    s0: 's0',
    s1: 's1',
    s2: 's2',
    s3: 's3',
    s4: 's4'
}

def slowness_func_names():
    return list(_slowness_func_names.values())

def get_slowness_func_name(s):
    return _slowness_func_names[s]

_soln_func_map = {
    _slowness_funcs[i]: _soln_funcs[i] for i in range(len(_slowness_funcs))}
def get_soln_func(s):
    return _soln_func_map[s]

_slowness_funcs_by_name = {v: k for k, v in _slowness_func_names.items()}
def get_slowness_func_by_name(name):
    if name not in _slowness_funcs_by_name:
        raise Exception("'%s' is not a valid slowness function" % name)
    return _slowness_funcs_by_name[name]

