import numpy as np

r = lambda x, y: np.sqrt(np.power(x, 2) + np.power(y, 2))

s1 = lambda x, y: 1 - np.sin(r(x, y))
f1 = lambda x, y: np.cos(r(x, y)) + r(x, y) - 1;

s2 = lambda x, y: np.abs(x + y)
f2 = lambda x, y: np.power(x + y, 2)/(2*np.sqrt(2))

s3 = r;
f3 = lambda x, y: (np.power(x, 2) + np.power(y, 2))/2

_speed_funcs = [s1, s2, s3]
_soln_funcs = [f1, f2, f3]

for s in _speed_funcs:
    s.dim = 2

for f in _soln_funcs:
    f.dim = 2

_soln_func_map = {
    _speed_funcs[i]: _soln_funcs[i] for i in range(len(_speed_funcs))}

def get_soln_func(s):
    return _soln_func_map[s]
