import numpy as np

r = lambda x, y: np.sqrt(np.power(x, 2) + np.power(y, 2))

s1 = lambda x, y: 1 - np.sin(r(x, y))
f1 = lambda x, y: np.cos(r(x, y)) + r(x, y) - 1;

s2 = lambda x, y: np.abs(x + y)
f2 = lambda x, y: np.power(x + y, 2)/(2*np.sqrt(2))

s3 = r;
f3 = lambda x, y: (np.power(x, 2) + np.power(y, 2))/2

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

# TODO: star_s

_speed_funcs = [s1, s2, s3, s5, s6, s8, s9]

_soln_funcs = [f1, f2, f3, f5, f6, f8, f9]

for func in _speed_funcs + _soln_funcs:
    func.dim = 2

def speed_funcs():
    return _speed_funcs

_speed_func_names = {
    s1: 's1',
    s2: 's2',
    s3: 's3',
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
