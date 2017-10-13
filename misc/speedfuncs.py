r = lambda x, y: np.sqrt(np.power(x, 2) + np.power(y, 2))

s1 = lambda x, y: 1 - np.sin(r(x, y))
f1 = lambda x, y: np.cos(r(x, y)) + r(x, y) - 1;

s2 = lambda x, y: np.abs(x + y)
f2 = lambda x, y: np.power(x + y, 2)/(2*np.sqrt(2))

s3 = r;
f3 = lambda x, y: (np.power(x, 2) + np.power(y, 2))/2
