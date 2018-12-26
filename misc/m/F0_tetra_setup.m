u0 = 1.9915638315627207;
u1 = 1.4142135623730949;
u2 = 2;
s = 1;
s0 = 1;
s1 = 1;
s2 = 1;
h = 1;
sfac = h*rand;

P = rand(3);
% P = [1 1 1; 0 1 0; 0 0 1];
% P = eye(3);
% P = eye(3) + h*h*rand(3);
p0 = P(:, 1);
p1 = P(:, 2);
p2 = P(:, 3);

theta = 0.5;

% rfac = 0.01;
% pfac = (rfac/h)*randn(3, 1);
pfac = P(:, 1);

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
n = @(x) p(x)/l(x);

pfac = @(x) p(x) - pfac;
lfac = @(x) norm(pfac(x));
T = @(x) sfac*h*lfac(x);
nfac = @(x) (lfac(x) > 1e1*eps)*(pfac(x)/max(eps, lfac(x)));

T0 = T([0; 0]);
T1 = T([1; 0]);
T2 = T([0; 1]);

tau0 = u0 - T0;
tau1 = u1 - T1;
tau2 = u2 - T2;
dtau = [tau1 - tau0; tau2 - tau0];
tau = @(x) tau0 + dot(dtau, x);

N = length(p0);
I = eye(N);
A = [-eye(N - 1); ones(1, N - 1)]; % constraint matrix
b = [zeros(N - 1, 1); 1];
