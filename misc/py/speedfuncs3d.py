import numpy as np

r = lambda x, y, z: np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))

def s0(x, y, z): return 1
s0 = np.vectorize(s0)
f0 = r

s1 = lambda x, y, z: 1 - np.sin(r(x, y, z))
f1 = lambda x, y, z: np.cos(r(x, y, z)) + r(x, y, z) - 1;

s2 = r;
f2 = lambda x, y, z: (x**2 + y**2 + z**2)/2

A = np.array([
    [1,   1/4, 1/8],
    [1/4, 1,   1/4],
    [1/8, 1/4, 1]])
alpha = np.pi/5
C = lambda x, y, z: np.array([np.cos(alpha*x), np.cos(alpha*y), np.cos(alpha*z)])
S = lambda x, y, z: np.array([np.sin(alpha*x), np.sin(alpha*y), np.sin(alpha*z)])
s3 = lambda x, y, z: np.linalg.norm(alpha*np.diag(C(x, y, z))@(A + A.T)@S(x, y, z))
s3 = np.vectorize(s3)
f3 = lambda x, y, z: S(x, y, z)@A@S(x, y, z)
f3 = np.vectorize(f3)

Bsqrt = A
s4 = np.vectorize(lambda x, y, z: np.linalg.norm(Bsqrt@[x, y, z]))
f4 = np.vectorize(lambda x, y, z: [x, y, z]@Bsqrt@[x, y, z]/2)

_slowness_funcs = [s0, s1, s2, s3, s4]

_soln_funcs = [f0, f1, f2, s3, s4]

def get_field(g, X, Y, Z):
    assert(X.shape == Y.shape)
    assert(X.shape == Z.shape)
    return g(X, Y, Z)

def get_fields(u, s, X, Y, Z):
    return get_field(u, X, Y, Z), get_field(s, X, Y, Z)

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
    return list(slowness_func_names.values())

def get_slowness_func_name(s):
    return _slowness_func_names[s]

_slowness_funcs_by_name = {
    _slowness_func_names[k]: k for k in _slowness_func_names.keys()}
def get_slowness_func_by_name(name):
    return _slowness_funcs_by_name[name]

for func in _slowness_funcs + _soln_funcs:
    func.dim = 3

_soln_func_map = {
    _slowness_funcs[i]: _soln_funcs[i] for i in range(len(_slowness_funcs))}

def get_soln_func(s):
    return _soln_func_map[s]
