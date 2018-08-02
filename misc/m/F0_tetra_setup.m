% u0 = 0.57428310090605839;
% u1 = 0.57434880146896394;
% u2 = 0.57442211198750626;
% s = 0.0028667363931870193;
% s0 = 0.0020480475244047947;
% s1 = 0.0012369806208707423;
% s2 = 0.0028667363931870193;
% h = 0.040000000000000001;

u0 = 0.5514651482575854;
u1 = 0.5419072788623589;
u2 = 0.5415495962762169;
s = 0.2382400185837108;
s0 = 0.2420057484596495;
s1 = 0.2740574756206167;
s2 = 0.2914785686731952;
h = 0.1;
theta = 0.5;
sfac = h*rand;
rfac = 0.01;

p0 = [1; 0; 0];
p1 = [0; 1; 0];
p2 = [0; 0; 1];

pfac = (rfac/h)*randn(3, 1);

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
nfac = @(x) pfac(x)/lfac(x);

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
