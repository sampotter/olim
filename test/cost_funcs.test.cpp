#include <gtest/gtest.h>
#include <random>

#include "cost_funcs.hpp"

TEST (cost_funcs, F0_works) {
  double u[3], h, s_hat, s[3], theta, p[3][3], f, df[2], d2f[2][2];
  double lam[2] = {1./3, 1./3};
  {
    u[0] = 0.1233189348351655;
    u[1] = 0.1839077882824167;
    u[2] = 0.2399525256649028;
    h = 0.4172670690843695;
    s_hat = 0.04965443032574213;
    s[0] = 0.9027161099152811;
    s[1] = 0.944787189721646;
    s[2] = 0.4908640924680799;
    theta = 0.4892526384000189;
    p[0][0] = -0.7981635845641424;
    p[0][1] = 1.018685282128575;
    p[0][2] = -0.1332174795077347;
    p[1][0] = -0.7145301637871584;
    p[1][1] = 1.351385768426657;
    p[1][2] = -0.2247710560525841;
    p[2][0] = -0.5890290307208013;
    p[2][1] = -0.2937535977354161;
    p[2][2] = -0.8479262436379339;
    F0<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.3629059042422889, 1e-14);
    ASSERT_NEAR(df[0], 0.093850695620551122, 1e-14);
    ASSERT_NEAR(df[1], -0.005830975252270762, 1e-14);
    ASSERT_NEAR(d2f[0][0], 0.013984781373164057, 1e-14);
    ASSERT_NEAR(d2f[0][1], -0.03387057681332333, 1e-14);
    ASSERT_NEAR(d2f[1][0], -0.03387057681332333, 1e-14);
    ASSERT_NEAR(d2f[1][1], 0.2802160611291684, 1e-14);
  }
  {
    u[0] = 0.6240600881736895;
    u[1] = 0.6791355408657477;
    u[2] = 0.395515215668593;
    h = 0.3674366485444766;
    s_hat = 0.9879820031616328;
    s[0] = 0.03773886623955214;
    s[1] = 0.8851680082024753;
    s[2] = 0.913286827639239;
    theta = 0.796183873585212;
    p[0][0] = -1.66416447498706;
    p[0][1] = -0.5900345642052215;
    p[0][2] = -0.2780641637653093;
    p[1][0] = 0.4227156912204783;
    p[1][1] = -1.67020069785047;
    p[1][2] = 0.4716343264163027;
    p[2][0] = -1.212847199674459;
    p[2][1] = 0.06619004842461142;
    p[2][2] = 0.652355888661374;
    F0<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.8529368893174771, 1e-14);
    ASSERT_NEAR(df[0], -0.1025834491128786, 1e-14);
    ASSERT_NEAR(df[1], -0.3596024800247133, 1e-14);
    ASSERT_NEAR(d2f[0][0], 1.272103018002429, 1e-14);
    ASSERT_NEAR(d2f[0][1], 0.135763147743216, 1e-14);
    ASSERT_NEAR(d2f[1][0], 0.135763147743216, 1e-14);
    ASSERT_NEAR(d2f[1][1], 0.2751072086015301, 1e-14);
  }
  {
    u[0] = 0.7150370784006941;
    u[1] = 0.9037205605563163;
    u[2] = 0.8909225043307892;
    h = 0.3341630527374962;
    s_hat = 0.6987458323347945;
    s[0] = 0.1978098266859292;
    s[1] = 0.03054094630463666;
    s[2] = 0.7440742603674624;
    theta = 0.5000224355902009;
    p[0][0] = -0.0549189146094067;
    p[0][1] = 0.9111272656538598;
    p[0][2] = 0.5945836974090524;
    p[1][0] = 0.3502011738745352;
    p[1][1] = 1.250251228304996;
    p[1][2] = 0.9297894585577157;
    p[2][0] = 0.23976325705858;
    p[2][1] = -0.6903611031112258;
    p[2][2] = -0.651553641750281;
    F0<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.9386573623510739, 1e-14);
    ASSERT_NEAR(df[0], 0.2848236230267062, 1e-14);
    ASSERT_NEAR(df[1], -0.1374461352207875, 1e-14);
    ASSERT_NEAR(d2f[0][0], 0.02146607964532761, 1e-14);
    ASSERT_NEAR(d2f[0][1], 0.05433382588472332, 1e-14);
    ASSERT_NEAR(d2f[1][0], 0.0543338258847233, 1e-14);
    ASSERT_NEAR(d2f[1][1], 0.2412036300850062, 1e-14);
  }
  {
    u[0] = 0.8865119330761013;
    u[1] = 0.0286741524641061;
    u[2] = 0.4899013885122239;
    h = 0.1679271456822568;
    s_hat = 0.9786806496411588;
    s[0] = 0.7126944716789141;
    s[1] = 0.500471624154843;
    s[2] = 0.4710883745419393;
    theta = 0.05961886757963919;
    p[0][0] = 0.5811723226759228;
    p[0][1] = -2.192434919965905;
    p[0][2] = -2.319280306643302;
    p[1][0] = 0.07993371029843969;
    p[1][1] = -0.9484809835705053;
    p[1][2] = 0.4114906214233742;
    p[2][0] = 0.6769778056840295;
    p[2][1] = 0.8577325452053552;
    p[2][2] = -0.6911591253829914;
    F0<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.6663944579421156, 1e-14);
    ASSERT_NEAR(df[0], -1.315913545982157, 1e-14);
    ASSERT_NEAR(df[1], -0.8745194627136645, 1e-14);
    ASSERT_NEAR(d2f[0][0], 0.1394565419621314, 1e-14);
    ASSERT_NEAR(d2f[0][1], -0.04419259725744203, 1e-14);
    ASSERT_NEAR(d2f[1][0], -0.04419259725744201, 1e-14);
    ASSERT_NEAR(d2f[1][1], 0.3964897378720262, 1e-14);
  }
  {
    u[0] = 0.6596052529083072;
    u[1] = 0.5185949425105382;
    u[2] = 0.9729745547638625;
    h = 0.6489914927123561;
    s_hat = 0.8003305753524015;
    s[0] = 0.4537977087269195;
    s[1] = 0.4323915037834617;
    s[2] = 0.8253137954020456;
    theta = 0.08346981485891403;
    p[0][0] = -0.5045864055140099;
    p[0][1] = -1.27059444980866;
    p[0][2] = -0.3825848027076484;
    p[1][0] = 0.6486792620486207;
    p[1][1] = 0.8257271492417582;
    p[1][2] = -1.014943642680137;
    p[2][0] = -0.4710699126831666;
    p[2][1] = 0.1370248741300503;
    p[2][2] = -0.2918633757535734;
    F0<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 1.012456254930158, 1e-14);
    ASSERT_NEAR(df[0], -0.1277064300798806, 1e-14);
    ASSERT_NEAR(df[1], 0.140073786123401, 1e-14);
    ASSERT_NEAR(d2f[0][0], 5.327886039147272, 1e-14);
    ASSERT_NEAR(d2f[0][1], 2.558835075182847, 1e-14);
    ASSERT_NEAR(d2f[1][0], 2.558835075182847, 1e-14);
    ASSERT_NEAR(d2f[1][1], 1.630350556702507, 1e-14);
  }
  {
    u[0] = 0.683363243294653;
    u[1] = 0.5465931145903228;
    u[2] = 0.4257288418711879;
    h = 0.6444427814313365;
    s_hat = 0.6476176301726844;
    s[0] = 0.6790167540932019;
    s[1] = 0.6357867105140836;
    s[2] = 0.9451741131094014;
    theta = 0.2089349224260229;
    p[0][0] = 0.7595683259147834;
    p[0][1] = -0.6572012990983503;
    p[0][2] = -0.6039184813761692;
    p[1][0] = 0.1769468223294113;
    p[1][1] = -0.3075034698627506;
    p[1][2] = -0.1318203529158936;
    p[2][0] = 0.5953576738841018;
    p[2][1] = 1.046832784305232;
    p[2][2] = -0.197958632611842;
    F0<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.8102532247452363, 1e-14);
    ASSERT_NEAR(df[0], -0.450287151146212, 1e-14);
    ASSERT_NEAR(df[1], -0.3755463418933428, 1e-14);
    ASSERT_NEAR(d2f[0][0], 0.113127418959527, 1e-14);
    ASSERT_NEAR(d2f[0][1], 0.493683542352749, 1e-14);
    ASSERT_NEAR(d2f[1][0], 0.4936835423527491, 1e-14);
    ASSERT_NEAR(d2f[1][1], 2.177922861450874, 1e-14);
  }
  {
    u[0] = 0.6877960851201071;
    u[1] = 0.3592282104018606;
    u[2] = 0.7363400743012017;
    h = 0.3947074752787632;
    s_hat = 0.6834158669679784;
    s[0] = 0.704047430334266;
    s[1] = 0.4423054133833708;
    s[2] = 0.01957762355331871;
    theta = 0.330857880214071;
    p[0][0] = -0.2700688126480988;
    p[0][1] = -0.4381413556021437;
    p[0][2] = -0.4086743147963257;
    p[1][0] = 0.9835452352055563;
    p[1][1] = -0.297697144009373;
    p[1][2] = 1.143678910770958;
    p[2][0] = -0.5316201175070693;
    p[2][1] = 0.9725657280086532;
    p[2][2] = -0.5222504849935489;
    F0<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.6227120026289962, 1e-14);
    ASSERT_NEAR(df[0], 0.04457909279171696, 1e-14);
    ASSERT_NEAR(df[1], 0.2139694297952391, 1e-14);
    ASSERT_NEAR(d2f[0][0], 2.644746157894539, 1e-14);
    ASSERT_NEAR(d2f[0][1], -2.763759804797053, 1e-14);
    ASSERT_NEAR(d2f[1][0], -2.763759804797053, 1e-14);
    ASSERT_NEAR(d2f[1][1], 2.951818740669762, 1e-14);
  }
  {
    u[0] = 0.5312092935824387;
    u[1] = 0.1088179382730454;
    u[2] = 0.6317663735284889;
    h = 0.1264998653293029;
    s_hat = 0.1343033043135746;
    s[0] = 0.09859409271099773;
    s[1] = 0.1420272484319284;
    s[2] = 0.1682512984915278;
    theta = 0.1962489222569553;
    p[0][0] = -0.2856861376018336;
    p[0][1] = -0.4624215443537469;
    p[0][2] = -0.4097852147401654;
    p[1][0] = -0.5035389841318734;
    p[1][1] = 1.233297112644135;
    p[1][2] = 0.6103052175739498;
    p[2][0] = 0.05907215505265444;
    p[2][1] = -1.466946730428247;
    p[2][2] = -1.625803266396234;
    F0<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.4338486230780023, 1e-14);
    ASSERT_NEAR(df[0], -0.4465438244925446, 1e-14);
    ASSERT_NEAR(df[1], 0.1218356011000501, 1e-14);
    ASSERT_NEAR(d2f[0][0], 0.0572055551104876, 1e-14);
    ASSERT_NEAR(d2f[0][1], -0.0365541536716814, 1e-14);
    ASSERT_NEAR(d2f[1][0], -0.03655415367168141, 1e-14);
    ASSERT_NEAR(d2f[1][1], 0.03065044324699815, 1e-14);
  }
  {
    u[0] = 0.0773468081126768;
    u[1] = 0.9138004107795679;
    u[2] = 0.7067152176969306;
    h = 0.5577889667548762;
    s_hat = 0.3134289899365913;
    s[0] = 0.1662035629021507;
    s[1] = 0.6224972592798952;
    s[2] = 0.9879347349524954;
    theta = 0.1704320230568833;
    p[0][0] = -0.4493973598602503;
    p[0][1] = -0.08429207275342081;
    p[0][2] = -1.991997177603269;
    p[1][0] = 0.8412456663525811;
    p[1][1] = -0.4146589174554855;
    p[1][2] = 1.912180838363861;
    p[2][0] = -0.3908987322337722;
    p[2][1] = 0.4091820825473549;
    p[2][2] = -1.142428168121445;
    F0<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.6481998345186274, 1e-14);
    ASSERT_NEAR(df[0], 0.05761143312148909, 1e-14);
    ASSERT_NEAR(df[1], 0.451516212604574, 1e-14);
    ASSERT_NEAR(d2f[0][0], 1.011444013947157, 1e-14);
    ASSERT_NEAR(d2f[0][1], -0.09267001267113532, 1e-14);
    ASSERT_NEAR(d2f[1][0], -0.09267001267113535, 1e-14);
    ASSERT_NEAR(d2f[1][1], 0.09281724916229861, 1e-14);
  }
}

TEST (cost_funcs, F1_works) {
  double u[3], h, s_hat, s[3], theta, p[3][3], f, df[2], d2f[2][2];
  double lam[2] = {1./3, 1./3};
  {
    u[0] = 0.1789824793143351;
    u[1] = 0.3389556782477182;
    u[2] = 0.210145637043552;
    h = 0.5101525197652502;
    s_hat = 0.9063643232652148;
    s[0] = 0.6289239386523179;
    s[1] = 0.1015338888123122;
    s[2] = 0.3908547527263546;
    theta = 0.05461661522365757;
    p[0][0] = 0.00490371102195903;
    p[0][1] = -0.3020322984463663;
    p[0][2] = 1.813582130904333;
    p[1][0] = 0.9148517528410741;
    p[1][1] = -0.05708071591447714;
    p[1][2] = 1.30936207782942;
    p[2][0] = -1.044735848409628;
    p[2][1] = -0.3482668103352654;
    p[2][2] = 1.412561160729294;
    F1<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.9277417461397732, 1e-14);
    ASSERT_NEAR(df[0], -0.1133734994151532, 1e-14);
    ASSERT_NEAR(df[1], -0.140282612382954, 1e-14);
    ASSERT_NEAR(d2f[0][0], 0.25858675084532861, 1e-14);
    ASSERT_NEAR(d2f[0][1], -0.2734996587352037, 1e-14);
    ASSERT_NEAR(d2f[1][0], -0.2734996587352037, 1e-14);
    ASSERT_NEAR(d2f[1][1], 0.3365816389629368, 1e-14);
  }
  {
    u[0] = 0.9174938324161169;
    u[1] = 0.7135740115943158;
    u[2] = 0.61833738362194;
    h = 0.3432878902413454;
    s_hat = 0.9360273266897697;
    s[0] = 0.1247740406604926;
    s[1] = 0.7305853615057071;
    s[2] = 0.6464774324258138;
    theta = 0.833151985669295;
    p[0][0] = -0.4258678455717095;
    p[0][1] = 0.6361398914029675;
    p[0][2] = 0.7931780813895736;
    p[1][0] = -0.8983771141007473;
    p[1][1] = 0.1562448442775494;
    p[1][2] = 1.597253911277571;
    p[2][0] = 0.1124397105937409;
    p[2][1] = -0.3086249345470284;
    p[2][2] = 0.4566596091151506;
    F1<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.9552435580352215, 1e-14);
    ASSERT_NEAR(df[0], 0.1422049831297019, 1e-14);
    ASSERT_NEAR(df[1], -0.2732978213137259, 1e-14);
    ASSERT_NEAR(d2f[0][0], 0.3654487449203132, 1e-14);
    ASSERT_NEAR(d2f[0][1], 0.1019099328885245, 1e-14);
    ASSERT_NEAR(d2f[1][0], 0.1019099328885244, 1e-14);
    ASSERT_NEAR(d2f[1][1], -0.03492141241044397, 1e-14);
  }
  {
    u[0] = 0.3606365710022027;
    u[1] = 0.7565095435019443;
    u[2] = 0.4139007486901897;
    h = 0.4923451043849377;
    s_hat = 0.6947432331326101;
    s[0] = 0.9727338850797841;
    s[1] = 0.3277549604934068;
    s[2] = 0.8378031830785756;
    theta = 0.7390722272735281;
    p[0][0] = 1.648383657046169;
    p[0][1] = -2.028361960852765;
    p[0][2] = -0.4492568114449102;
    p[1][0] = 0.2359932017119493;
    p[1][1] = -0.835172963783281;
    p[1][2] = -1.275955163408218;
    p[2][0] = 0.6170350812805879;
    p[2][1] = 0.6127015037072735;
    p[2][2] = 0.2893811552948809;
    F1<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.9354846054549475, 1e-14);
    ASSERT_NEAR(df[0], -0.3697887584066532, 1e-14);
    ASSERT_NEAR(df[1], -0.9199344696831468, 1e-14);
    ASSERT_NEAR(d2f[0][0], 1.277447423581187, 1e-14);
    ASSERT_NEAR(d2f[0][1], 0.7950709805264409, 1e-14);
    ASSERT_NEAR(d2f[1][0], 0.795070980526441, 1e-14);
    ASSERT_NEAR(d2f[1][1], 0.7491310069183262, 1e-14);
  }
  {
    u[0] = 0.6604379663126019;
    u[1] = 0.04755467311386607;
    u[2] = 0.3487848085100589;
    h = 0.4513405803557432;
    s_hat = 0.2409049971201107;
    s[0] = 0.7150450132961765;
    s[1] = 0.8561822920062879;
    s[2] = 0.2815076951185533;
    theta = 0.7310508297237415;
    p[0][0] = -1.034424717234166;
    p[0][1] = 1.339554597904117;
    p[0][2] = -0.9691403558288297;
    p[1][0] = 0.2087156012721984;
    p[1][1] = -0.6185933553605473;
    p[1][2] = 0.5120156437364544;
    p[2][0] = 0.01135417110082006;
    p[2][1] = -0.04398862626670082;
    p[2][2] = 2.949092541315588;
    F1<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.5625665483580271, 1e-14);
    ASSERT_NEAR(df[0], -0.4544144376448905, 1e-14);
    ASSERT_NEAR(df[1], 0.2456250444201228, 1e-14);
    ASSERT_NEAR(d2f[0][0], 1.937369008430461, 1e-14);
    ASSERT_NEAR(d2f[0][1], 2.219065295204201, 1e-14);
    ASSERT_NEAR(d2f[1][0], 2.2190652952042, 1e-14);
    ASSERT_NEAR(d2f[1][1], 1.65725997787241, 1e-14);
  }
  {
    u[0] = 0.3531418129389555;
    u[1] = 0.4494435565717483;
    u[2] = 0.9635302868434269;
    h = 0.0422977979145428;
    s_hat = 0.9729583341406347;
    s[0] = 0.189206843122376;
    s[1] = 0.6671203000400749;
    s[2] = 0.5864396146789179;
    theta = 0.6751124164051153;
    p[0][0] = -0.2788611163878222;
    p[0][1] = 0.2451915937240211;
    p[0][2] = 1.472513478821577;
    p[1][0] = -2.275101542347401;
    p[1][1] = -1.633290724921242;
    p[1][2] = 0.4154689369787572;
    p[2][0] = -0.6547688518580145;
    p[2][1] = -0.2963483154065983;
    p[2][2] = -1.496918876447554;
    F1<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.6216362845903155, 1e-14);
    ASSERT_NEAR(df[0], 0.1809666396544651, 1e-14);
    ASSERT_NEAR(df[1], 0.6312877785721317, 1e-14);
    ASSERT_NEAR(d2f[0][0], 0.1203376694711638, 1e-14);
    ASSERT_NEAR(d2f[0][1], 0.1268138352144878, 1e-14);
    ASSERT_NEAR(d2f[1][0], 0.1268138352144878, 1e-14);
    ASSERT_NEAR(d2f[1][1], 0.2108079144710181, 1e-14);
  }
  {
    u[0] = 0.2684388213972831;
    u[1] = 0.2578461701126047;
    u[2] = 0.3316652387426293;
    h = 0.1522340128629465;
    s_hat = 0.3480076597161135;
    s[0] = 0.1216584543077261;
    s[1] = 0.8841530577496601;
    s[2] = 0.09427839006414607;
    theta = 0.9300406261074891;
    p[0][0] = -0.3744366065079624;
    p[0][1] = -1.451740719917382;
    p[0][2] = -0.6186818075674646;
    p[1][0] = 0.934500962295567;
    p[1][1] = 1.055929256492714;
    p[1][2] = 0.1602272767486467;
    p[2][0] = 0.2874003202431548;
    p[2][1] = 0.6329055770412402;
    p[2][2] = -1.459041869083612;
    F1<3, 2> func(h, theta);
    func.set_args(u, s_hat, s, p);
    func.set_lambda(lam);
    func.eval(f);
    func.grad(df);
    func.hess(d2f);
    ASSERT_NEAR(f, 0.3251021429753333, 1e-14);
    ASSERT_NEAR(df[0], 0.0708733369651883, 1e-14);
    ASSERT_NEAR(df[1], 0.1308033314110975, 1e-14);
    ASSERT_NEAR(d2f[0][0], 0.7016186994936204, 1e-14);
    ASSERT_NEAR(d2f[0][1], 0.5563210804125802, 1e-14);
    ASSERT_NEAR(d2f[1][0], 0.5563210804125802, 1e-14);
    ASSERT_NEAR(d2f[1][1], 0.2980850527574, 1e-14);
  }
}
