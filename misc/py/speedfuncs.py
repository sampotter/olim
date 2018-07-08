import numpy as np

sqrt2 = np.sqrt(2)

r = lambda x, y: np.sqrt(np.power(x, 2) + np.power(y, 2))

s0 = lambda x, y: np.ones(x.shape)
f0 = r

s1 = lambda x, y: 1 - np.sin(r(x, y))
f1 = lambda x, y: np.cos(r(x, y)) + r(x, y) - 1;

s2 = lambda x, y: np.abs(x + y)
f2 = lambda x, y: np.power(x + y, 2)/(2*sqrt2)

s3 = r;
f3 = lambda x, y: (np.power(x, 2) + np.power(y, 2))/2

s4xy = lambda x, y: np.sqrt(10*x*x - 6*sqrt2*x*y + 37*y*y)/4
f4xy = lambda x, y: (6*x*x - 2*sqrt2*x*y + 11*y*y)/16

s4xz = lambda x, z: np.sqrt(10*x*x + 6*sqrt2*x*z + 37*z*z)/4
f4xz = lambda x, z: (6*x*x + 2*sqrt2*x*z + 11*z*z)/16

s4yz = lambda y, z: np.sqrt(37*y*y + 54*y*z + 37*z*z)/4
f4yz = lambda y, z: (11*y*y + 10*y*z + 11*z*z)/16

s5 = lambda x, y: np.sqrt(np.power(x, 18) + np.power(y, 18))
f5 = lambda x, y: (np.power(x, 10) + np.power(y, 10))/10

s6 = lambda x, y: 4*np.sqrt(np.power(x - y, 2)*np.power(x + y, 2)*
                            (np.power(x, 2) + np.power(y, 2)))
f6 = lambda x, y: np.power(np.power(x, 2) - np.power(y, 2), 2)

# TODO: masha_s

s8 = lambda x, y: np.power(1 + np.sqrt(np.power(x, 2) + np.power(y, 2)), -2)
f8 = lambda x, y: 1 - np.divide(1, 1 + np.sqrt(np.power(x, 2) + np.power(y, 2)))

a = 10;
s9 = lambda x, y: 2*a*np.sqrt(
    np.exp(-2*a*(np.power(x, 2) + np.power(y, 2)))*
    (np.power(x, 2) + np.power(y, 2)))
f9 = lambda x, y: 1 - np.divide(1, np.exp(a*(np.power(x, 2) + np.power(y, 2))))

s_qv = lambda x, y: 1/(0.5 + 0.5*x)
f_qv = lambda x, y: 2*np.arccosh(1 + 0.25*np.sqrt(x*x + y*y)*s_qv(x, y))

# TODO: star_s

_speed_funcs = [s0, s1, s2, s3, s4xy, s4xz, s4yz, s5, s6, s8, s9]

_soln_funcs = [f0, f1, f2, f3, f4xy, s4xz, s4yz, f5, f6, f8, f9]

for func in _speed_funcs + _soln_funcs:
    func.dim = 2

def speed_funcs():
    return _speed_funcs

_speed_func_names = {
    s0: 's0',
    s1: 's1',
    s2: 's2',
    s3: 's3',
    s4xy: 's4xy',
    s4xz: 's4xz',
    s4yz: 's4yz',
    s5: 's5',
    s6: 's6',
    s8: 's8',
    s9: 's9'
}
def get_speed_func_name(s):
    return _speed_func_names[s]

_soln_func_map = {
    _speed_funcs[i]: _soln_funcs[i] for i in range(len(_speed_funcs))}
def get_soln_func(s):
    return _soln_func_map[s]

_speed_funcs_by_name = {v: k for k, v in _speed_func_names.items()}
def get_speed_func_by_name(name):
    if name not in _speed_funcs_by_name:
        raise Exception("'%s' is not a valid speed function" % name)
    return _speed_funcs_by_name[name]

