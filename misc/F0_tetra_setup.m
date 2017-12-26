u0 = rand;
u1 = rand;
u2 = rand;
h = rand;
s = rand;
s0 = rand;
s1 = rand;
s2 = rand;
theta = rand;

p0 = randn(3, 1);
p1 = randn(3, 1);
p2 = randn(3, 1);

du1 = u1 - u0;
du2 = u2 - u0;
du = [du1; du2];

dp1 = p1 - p0;
dp2 = p2 - p0;
dP = [dp1 dp2];

savg = (s0 + s1 + s2)/3;
stheta = (1 - theta)*s + theta*savg;
sh = h*stheta;
ds = [s1 - s0; s2 - s0];

u = @(x) (1 - sum(x))*u0 + x(1)*u1 + x(2)*u2;
p = @(x) (1 - sum(x))*p0 + x(1)*p1 + x(2)*p2;
q = @(x) p(x)'*p(x);
l = @(x) sqrt(q(x));
cprojp = @(x) eye(3) - p(x)*p(x)'/q(x);

N = length(p0);
A = [-eye(N - 1); ones(1, N - 1)]; % constraint matrix
b = [zeros(N - 1, 1); 1];