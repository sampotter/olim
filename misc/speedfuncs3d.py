import numpy as np

r = lambda x, y, z: np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))

s0 = lambda x, y, z: 1
f0 = r

s1 = lambda x, y, z: 1 - np.sin(r(x, y, z))
f1 = lambda x, y, z: np.cos(r(x, y, z)) + r(x, y, z) - 1;

s2 = lambda x, y, z: np.abs(x + y + z)
f2 = lambda x, y, z: np.power(x + y + z, 2)/(2*np.sqrt(3))

s3 = r;
f3 = lambda x, y, z: np.power(r(x, y, z), 2)/2

s5 = lambda x, y, z: np.power(x, 18) + np.power(y, 18) + np.power(z, 18)
f5 = lambda x, y, z: (np.power(x, 10) + np.power(y, 10) + np.power(z, 10))/10

s6 = lambda x, y, z: 4*np.sqrt(
    np.power(-np.power(x, 2) + np.power(y, 2) + np.power(z, 2), 2)*
    (np.power(x, 2) + np.power(y, 2) + np.power(z, 2)))
f6 = lambda x, y, z: np.power(
    np.power(x, 2) - np.power(y, 2) - np.power(z, 2), 2)

s8 = lambda x, y, z: np.power(
    1 + np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2)), -2)
f8 = lambda x, y, z: 1 - np.divide(
    1,
    1 + np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2)))

a = 10
s9 = lambda x, y, z: 2*a*np.sqrt(
    np.multiply(
        np.exp(-2*a*(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))),
        np.power(x, 2) + np.power(y, 2) + np.power(z, 2)))
f9 = 1 - np.divide(1, np.exp(np.power(x, 2) + np.power(y, 2) + np.power(z, 2)))

_speed_funcs = [s0, s1, s2, s3, s5, s6, s8, s9]
_soln_funcs = [f0, f1, f2, f3, f5, f6, f8, f9]

for func in _speed_funcs + _soln_funcs:
    func.dim = 3

_soln_func_map = {
    _speed_funcs[i]: _soln_funcs[i] for i in range(len(_speed_funcs))}

def get_soln_func(s):
    return _soln_func_map[s]
