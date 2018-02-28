import numpy as np

r = lambda x, y, z: np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))

def s0(x, y, z): return 1
s0 = np.vectorize(s0)
f0 = r

s1 = lambda x, y, z: 1 - np.sin(r(x, y, z))
f1 = lambda x, y, z: np.cos(r(x, y, z)) + r(x, y, z) - 1;

s2 = lambda x, y, z: np.abs(x + y + z)
f2 = lambda x, y, z: np.power(x + y + z, 2)/2

s3 = r;
f3 = lambda x, y, z: np.power(r(x, y, z), 2)/2

s5 = lambda x, y, z: np.sqrt(x**18 + y**18 + z**18)
f5 = lambda x, y, z: (x**10 + y**10 + z**10)/10

_speed_funcs = [s0, s1, s2, s3, s5]

_soln_funcs = [f0, f1, f2, f3, f5]

def speed_funcs():
    return _speed_funcs

_speed_func_names = {
    s0: 's0',
    s1: 's1',
    s2: 's2',
    s3: 's3',
    s5: 's5'
}

def get_speed_func_name(s):
    return _speed_func_names[s]

_speed_funcs_by_name = {
    _speed_func_names[k]: k for k in _speed_func_names.keys()}
def get_speed_func_by_name(name):
    return _speed_funcs_by_name[name]

for func in _speed_funcs + _soln_funcs:
    func.dim = 3

_soln_func_map = {
    _speed_funcs[i]: _soln_funcs[i] for i in range(len(_speed_funcs))}

def get_soln_func(s):
    return _soln_func_map[s]
