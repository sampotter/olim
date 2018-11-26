#include <gtest/gtest.h>
#include <random>

#include "cost_funcs.hpp"

#define P001 1
#define P010 2
#define P011 3
#define P100 4
#define P101 5
#define P110 6
#define P111 7

TEST (cost_funcs, rhr_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.8147236863931789;
    u1 = 0.9057919370756192;
    u2 = 0.1269868162935061;
    h = 0.9133758561390194;
    s = 0.6323592462254095;
    s0 = 0.09754040499940952;
    s1 = 0.2784982188670484;
    s2 = 0.5468815192049838;
    lam[0] = 0.5361628596899671;
    lam[1] = 0.4638371403100329;
    p0[0] = 3.57839693972576;
    p0[1] = 2.769437029884877;
    p0[2] = -1.349886940156521;
    p1[0] = 3.034923466331855;
    p1[1] = 0.7254042249461056;
    p1[2] = -0.06305487318965619;
    p2[0] = 0.7147429038260958;
    p2[1] = -0.2049660582997746;
    p2[2] = -0.1241443482163119;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<RHR, 3>(w, p0, p1, p2, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 1.689762490331113, 5e-16);
    ASSERT_NEAR(df[0], -0.4282612531697417, 2e-15);
    ASSERT_NEAR(df[1], -2.608936019769593, 2e-15);
    ASSERT_NEAR(d2f[0], 1.549985998542815, 2e-15);
    ASSERT_NEAR(d2f[1], 1.812653533777521, 2e-15);
    ASSERT_NEAR(d2f[2], 2.18064566066142, 2e-15);
  }
  {
    u0 = 0.4087198461125521;
    u1 = 0.5948960740086143;
    u2 = 0.2622117477808454;
    h = 0.602843089382083;
    s = 0.7112157804336829;
    s0 = 0.2217467340172401;
    s1 = 0.1174176508558059;
    s2 = 0.2966758732183269;
    lam[0] = 0.4578702935688697;
    lam[1] = 0.5421297064311303;
    p0[0] = -0.5078175502781737;
    p0[1] = -0.3205755066002393;
    p0[2] = 0.01246904136161795;
    p1[0] = -3.029177341404146;
    p1[1] = -0.4570146408715826;
    p1[2] = 1.242448406390738;
    p2[0] = -1.06670139898475;
    p2[1] = 0.9337281626712385;
    p2[2] = 0.3503210013561116;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<RHR, 3>(w, p0, p1, p2, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 1.326701745236962, 2.22045e-16);
    ASSERT_NEAR(df[0], 1.364704831564247, 5e-16);
    ASSERT_NEAR(df[1], 0.2015691301705069, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.06713366527082466, 2.22045e-16);
    ASSERT_NEAR(d2f[1], -0.1164792876214269, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.2701885620146574, 2.22045e-16);
  }
}

TEST (cost_funcs, mp0_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.2510838579760311;
    u1 = 0.6160446761466392;
    u2 = 0.4732888489027293;
    h = 0.3516595070629968;
    s = 0.8308286278962909;
    s0 = 0.5852640911527243;
    s1 = 0.5497236082911395;
    s2 = 0.91719366382981;
    lam[0] = 0.877898161628302;
    lam[1] = 0.1221018383716979;
    p0[0] = -0.8044659563495471;
    p0[1] = 0.6966244158496073;
    p0[2] = 0.8350881650726819;
    p1[0] = -0.2437151403779522;
    p1[1] = 0.2156700864037444;
    p1[2] = -1.165843931482049;
    p2[0] = -1.147952778898594;
    p2[1] = 0.104874716016494;
    p2[2] = 0.7222540322250016;
    F_wkspc<MP0, 2> w;
    set_args<MP0, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP0, 3>(w, p0, p1, p2, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.8703896887718838, 2.22045e-16);
    ASSERT_NEAR(df[0], 0.776304163583716, 5e-16);
    ASSERT_NEAR(df[1], 0.2502826994447658, 9e-16);
    ASSERT_NEAR(d2f[0], 0.5650921770320323, 9e-16);
    ASSERT_NEAR(d2f[1], 0.04045877564956703, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.1226370621333253, 2.22045e-16);
  }
  {
    u0 = 0.2372835797715215;
    u1 = 0.4588488281799311;
    u2 = 0.963088539286913;
    h = 0.546805718738968;
    s = 0.5211358308040015;
    s0 = 0.2315943867085238;
    s1 = 0.4888977439201669;
    s2 = 0.6240600881736895;
    lam[0] = 0.4384777227796616;
    lam[1] = 0.5615222772203384;
    p0[0] = 0.426387557408945;
    p0[1] = -0.3728087417235042;
    p0[2] = -0.2364545837571863;
    p1[0] = 2.023690886603053;
    p1[1] = -2.258353970496191;
    p1[2] = 2.229445680456899;
    p2[0] = 0.3375637006131064;
    p2[1] = 1.000060819589125;
    p2[2] = -1.66416447498706;
    F_wkspc<MP0, 2> w;
    set_args<MP0, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP0, 3>(w, p0, p1, p2, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 1.04937627761328, 5e-16);
    ASSERT_NEAR(df[0], 0.8235524285361199, 2.22045e-16);
    ASSERT_NEAR(df[1], 0.5554223656236505, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 1.605676932937447, 2e-15);
    ASSERT_NEAR(d2f[1], -1.094595303883586, 2e-15);
    ASSERT_NEAR(d2f[2], 0.8037340187570788, 4e-16);
  }
}

TEST (cost_funcs, mp1_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.2077422927330285;
    u1 = 0.3012463302794907;
    u2 = 0.4709233485175907;
    h = 0.2304881602115585;
    s = 0.8443087926953891;
    s0 = 0.1947642895670493;
    s1 = 0.2259217809723988;
    s2 = 0.1707080471478586;
    lam[0] = 0.3009860874640407;
    lam[1] = 0.6990139125359593;
    p0[0] = -1.174212331456816;
    p0[1] = -0.1922395175392748;
    p0[2] = -0.2740702299326022;
    p1[0] = 1.530072514424096;
    p1[1] = -0.2490247425137138;
    p1[2] = -1.064213412889327;
    p2[0] = 1.603457298120044;
    p2[1] = 1.234679146890778;
    p2[2] = -0.2296264509631805;
    F_wkspc<MP1, 2> w;
    set_args<MP1, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP1, 3>(w, p0, p1, p2, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.6375555765979776, 2.22045e-16);
    ASSERT_NEAR(df[0], 0.3994996388285198, 2.22045e-16);
    ASSERT_NEAR(df[1], 0.6149253208800269, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.1218408945664885, 2.22045e-16);
    ASSERT_NEAR(d2f[1], -0.006797926693548648, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.0317874368395807, 2.22045e-16);
  }
  {
    u0 = 0.679727951377338;
    u1 = 0.1365531373553697;
    u2 = 0.7212274985817402;
    h = 0.1067618616072414;
    s = 0.6537573486685596;
    s0 = 0.4941739366392701;
    s1 = 0.7790517232312751;
    s2 = 0.7150370784006941;
    lam[0] = 0.5973384897627685;
    lam[1] = 0.4026615102372315;
    p0[0] = 1.082633504236756;
    p0[1] = 1.006077110819051;
    p0[2] = -0.6509077365977526;
    p1[0] = 0.2570561574339689;
    p1[1] = -0.944377806404219;
    p1[2] = -1.321788521392564;
    p2[0] = 0.9248259334937059;
    p2[1] = 4.984907525081331e-05;
    p2[2] = -0.0549189146094067;
    F_wkspc<MP1, 2> w;
    set_args<MP1, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP1, 3>(w, p0, p1, p2, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.456074668933617, 2.22045e-16);
    ASSERT_NEAR(df[0], -0.4449408124402994, 2.22045e-16);
    ASSERT_NEAR(df[1], 0.05475042537160603, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.285578449160759, 2.22045e-16);
    ASSERT_NEAR(d2f[1], 0.1262558215206283, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.09340734618321611, 2.22045e-16);
  }
}

TEST (cost_funcs, rhr111_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.6554778901775566;
    u1 = 0.1711866878115618;
    u2 = 0.7060460880196088;
    h = 0.03183284637742068;
    s = 0.27692298496089;
    s0 = 0.04617139063115394;
    s1 = 0.09713178123584754;
    s2 = 0.8234578283272926;
    lam[0] = 0.6866383302239732;
    lam[1] = 0.3133616697760268;
    p0[0] = 1;
    p0[1] = 0;
    p0[2] = 0;
    p1[0] = 0;
    p1[1] = 1;
    p1[2] = 0;
    p2[0] = 0;
    p2[1] = 0;
    p2[2] = 1;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3, P100, P010, P001>(w, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<RHR, 3, P100, P010, P001>(w, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.3454445478969491, 2.22045e-16);
    ASSERT_NEAR(df[0], -0.4762716207665988, 2.22045e-16);
    ASSERT_NEAR(df[1], 0.0542281005452294, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.01369271615092746, 2.22045e-16);
    ASSERT_NEAR(d2f[1], 0.007268089138334956, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.02134573675679246, 2.22045e-16);
  }
  {
    u0 = 0.9502220488383549;
    u1 = 0.03444608050290876;
    u2 = 0.4387443596563982;
    h = 0.3815584570930084;
    s = 0.7655167881490024;
    s0 = 0.7951999011370632;
    s1 = 0.1868726045543786;
    s2 = 0.4897643957882311;
    lam[0] = 0.4080836365614623;
    lam[1] = 0.5919163634385377;
    p0[0] = 1;
    p0[1] = 0;
    p0[2] = 0;
    p1[0] = 0;
    p1[1] = 1;
    p1[2] = 0;
    p2[0] = 0;
    p2[1] = 0;
    p2[2] = 1;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3, P100, P010, P001>(w, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<RHR, 3, P100, P010, P001>(w, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.4837561772031396, 2.22045e-16);
    ASSERT_NEAR(df[0], -0.7499842322661709, 2.22045e-16);
    ASSERT_NEAR(df[1], -0.2710004126874189, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.6816476084460316, 2.22045e-16);
    ASSERT_NEAR(d2f[1], 0.2164153342738836, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.5371594292454075, 2.22045e-16);
  }
}

TEST (cost_funcs, rhr_123_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.6160446761466392;
    u1 = 0.4732888489027293;
    u2 = 0.3516595070629968;
    h = 0.8308286278962909;
    s = 0.5852640911527243;
    s0 = 0.5497236082911395;
    s1 = 0.91719366382981;
    s2 = 0.2858390188203735;
    lam[0] = 0.5011486754471214;
    lam[1] = 0.4988513245528786;
    p0[0] = 1;
    p0[1] = 0;
    p0[2] = 0;
    p1[0] = 1;
    p1[1] = 1;
    p1[2] = 0;
    p2[0] = 1;
    p2[1] = 1;
    p2[2] = 1;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3, P100, P110, P111>(w, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<RHR, 3, P100, P110, P111>(w, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 1.141809140746093, 2.22045e-16);
    ASSERT_NEAR(df[0], 0.1814962984200262, 2.22045e-16);
    ASSERT_NEAR(df[1], 0.2216205589568347, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.180066544395251, 2.22045e-16);
    ASSERT_NEAR(d2f[1], 0.1081393761979407, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.3245833705380006, 2.22045e-16);
  }
}

TEST (cost_funcs, rhr_222_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.3804458469753567;
    u1 = 0.5678216407252211;
    u2 = 0.07585428956306361;
    h = 0.05395011866660715;
    s = 0.5307975530089727;
    s0 = 0.7791672301020112;
    s1 = 0.934010684229183;
    s2 = 0.1299062084737301;
    lam[0] = 0.5478865585019908;
    lam[1] = 0.4521134414980094;
    p0[0] = 1;
    p0[1] = 1;
    p0[2] = 0;
    p1[0] = 1;
    p1[1] = 0;
    p1[2] = 1;
    p2[0] = 0;
    p2[1] = 1;
    p2[2] = 1;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3, P110, P101, P011>(w, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<RHR, 3, P110, P101, P011>(w, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.3805226825467133, 2.22045e-16);
    ASSERT_NEAR(df[0], 0.2001667619329952, 2.22045e-16);
    ASSERT_NEAR(df[1], -0.2940365099200196, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.04203427274029952, 2.22045e-16);
    ASSERT_NEAR(d2f[1], 0.01950245491928368, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.04352034287973601, 2.22045e-16);
  }
}

TEST (cost_funcs, mp0_111_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.1492940055590575;
    u1 = 0.2575082541237365;
    u2 = 0.8407172559836625;
    h = 0.254282178971531;
    s = 0.8142848260688164;
    s0 = 0.2435249687249893;
    s1 = 0.9292636231872278;
    s2 = 0.3499837659848087;
    lam[0] = 0.4391432317015156;
    lam[1] = 0.5608567682984844;
    p0[0] = 1;
    p0[1] = 0;
    p0[2] = 0;
    p1[0] = 0;
    p1[1] = 1;
    p1[2] = 0;
    p2[0] = 0;
    p2[1] = 0;
    p2[2] = 1;
    F_wkspc<MP0, 2> w;
    set_args<MP0, 3, P100, P010, P001>(w, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP0, 3, P100, P010, P001>(w, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.7043216835245574, 2.22045e-16);
    ASSERT_NEAR(df[0], 0.2118249098616238, 2.22045e-16);
    ASSERT_NEAR(df[1], 0.8237507838152849, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.3822049572774064, 2.22045e-16);
    ASSERT_NEAR(d2f[1], 0.1214133062583364, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.3256096267854676, 2.22045e-16);
  }
}

TEST (cost_funcs, mp0_123_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.7093648308580726;
    u1 = 0.7546866819823609;
    u2 = 0.2760250769985784;
    h = 0.6797026768536748;
    s = 0.6550980039738407;
    s0 = 0.1626117351946306;
    s1 = 0.1189976815583766;
    s2 = 0.498364051982143;
    lam[0] = 0.7381909431454695;
    lam[1] = 0.2618090568545305;
    p0[0] = 1;
    p0[1] = 0;
    p0[2] = 0;
    p1[0] = 1;
    p1[1] = 1;
    p1[2] = 0;
    p2[0] = 1;
    p2[1] = 1;
    p2[2] = 1;
    F_wkspc<MP0, 2> w;
    set_args<MP0, 3, P100, P110, P111>(w, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP0, 3, P100, P110, P111>(w, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 1.076654193038217, 2.22045e-16);
    ASSERT_NEAR(df[0], 0.2615538762816094, 2.22045e-16);
    ASSERT_NEAR(df[1], -0.1604962261339898, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.1116985818210038, 2.22045e-16);
    ASSERT_NEAR(d2f[1], 0.08433077961136608, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.2660298640743628, 2.22045e-16);    
  }
}

TEST (cost_funcs, mp0_222_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.5852677509797773;
    u1 = 0.223811939491137;
    u2 = 0.7512670593056529;
    h = 0.2550951154592691;
    s = 0.5059570516651424;
    s0 = 0.699076722656686;
    s1 = 0.8909032525357985;
    s2 = 0.9592914252054443;
    lam[0] = 0.7978764021813243;
    lam[1] = 0.2021235978186757;
    p0[0] = 1;
    p0[1] = 1;
    p0[2] = 0;
    p1[0] = 1;
    p1[1] = 0;
    p1[2] = 1;
    p2[0] = 0;
    p2[1] = 1;
    p2[2] = 1;
    F_wkspc<MP0, 2> w;
    set_args<MP0, 3, P110, P101, P011>(w, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP0, 3, P110, P101, P011>(w, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.5543810063040289, 2.22045e-16);
    ASSERT_NEAR(df[0], -0.254931255983449, 2.22045e-16);
    ASSERT_NEAR(df[1], 0.1929848493894792, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.2163522824130686, 2.22045e-16);
    ASSERT_NEAR(d2f[1], 0.1206745491311252, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.2637686030153841, 2.22045e-16);
  }
}

TEST (cost_funcs, mp1_111_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.0119020695012414;
    u1 = 0.3371226443988815;
    u2 = 0.1621823081932428;
    h = 0.794284540683907;
    s = 0.3112150420448049;
    s0 = 0.5285331355062127;
    s1 = 0.1656487294997809;
    s2 = 0.6019819414016365;
    lam[0] = 0.2867577282668127;
    lam[1] = 0.7132422717331873;
    p0[0] = 1;
    p0[1] = 0;
    p0[2] = 0;
    p1[0] = 0;
    p1[1] = 1;
    p1[2] = 0;
    p2[0] = 0;
    p2[1] = 0;
    p2[2] = 1;
    F_wkspc<MP1, 2> w;
    set_args<MP1, 3, P100, P010, P001>(w, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP1, 3, P100, P010, P001>(w, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.4529430571954809, 2.22045e-16);
    ASSERT_NEAR(df[0], 0.3311834545123941, 2.22045e-16);
    ASSERT_NEAR(df[1], 0.4630909441088277, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.650101197469628, 2.22045e-16);
    ASSERT_NEAR(d2f[1], 0.143392189981853, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.5179182309286099, 2.22045e-16);
  }
}

TEST (cost_funcs, mp1_123_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.6892145031400078;
    u1 = 0.7481515928237095;
    u2 = 0.4505415985024978;
    h = 0.08382137799693257;
    s = 0.2289769687168188;
    s0 = 0.9133373615016696;
    s1 = 0.152378018969223;
    s2 = 0.8258169774895474;
    lam[0] = 0.3508311835064203;
    lam[1] = 0.6491688164935797;
    p0[0] = 1;
    p0[1] = 0;
    p0[2] = 0;
    p1[0] = 1;
    p1[1] = 1;
    p1[2] = 0;
    p2[0] = 1;
    p2[1] = 1;
    p2[2] = 1;
    F_wkspc<MP1, 2> w;
    set_args<MP1, 3, P100, P110, P111>(w, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP1, 3, P100, P110, P111>(w, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.6083344569269014, 2.22045e-16);
    ASSERT_NEAR(df[0], 0.03135544872856954, 2.22045e-16);
    ASSERT_NEAR(df[1], -0.2080235671311109, 2.22045e-16);
    ASSERT_NEAR(d2f[0], -0.02804903602251281, 2.22045e-16);
    ASSERT_NEAR(d2f[1], -0.02912624593269977, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.01155464589786689, 2.22045e-16);
  }
}

TEST (cost_funcs, mp1_222_works) {
  double u0, u1, u2, h, s, s0, s1, s2, lam[2], p0[3], p1[3], p2[3], f, df[2], d2f[3];
  {
    u0 = 0.07817552875318368;
    u1 = 0.4426782697754463;
    u2 = 0.1066527701805844;
    h = 0.9618980808550537;
    s = 0.004634224134067444;
    s0 = 0.7749104647115024;
    s1 = 0.817303220653433;
    s2 = 0.8686947053635097;
    lam[0] = 0.1743755070300453;
    lam[1] = 0.8256244929699547;
    p0[0] = 1;
    p0[1] = 1;
    p0[2] = 0;
    p1[0] = 1;
    p1[1] = 0;
    p1[2] = 1;
    p2[0] = 0;
    p2[1] = 1;
    p2[2] = 1;
    F_wkspc<MP1, 2> w;
    set_args<MP1, 3, P110, P101, P011>(w, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP1, 3, P110, P101, P011>(w, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.7091953080104028, 2.22045e-16);
    ASSERT_NEAR(df[0], 0.446582281851884, 2.22045e-16);
    ASSERT_NEAR(df[1], 0.3498090802809615, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 0.635221544278329, 2.22045e-16);
    ASSERT_NEAR(d2f[1], 0.3098742845386949, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 0.5658542787531913, 2.22045e-16);
  }
}

TEST (cost_funcs, rhr_fac_works) {
  double u0, u1, u2, h, s, s0, s1, s2, sf, lam[2], p0[3], p1[3], p2[3], pf[3],
    f, df[2], d2f[3];
  {
    u0 = 0.5078582846611182;
    u1 = 0.08551579709004398;
    u2 = 0.2624822346983327;
    h = 0.8010146227697388;
    s = 0.02922027756214629;
    s0 = 0.9288541394780446;
    s1 = 0.7303308628554529;
    s2 = 0.4886089738035791;
    sf = 0.5785250610234389;
    lam[0] = 0.04089130472927549;
    lam[1] = 0.9591086952707245;
    p0[0] = -1.565056014150725;
    p0[1] = -0.08453947981772419;
    p0[2] = 1.60394635060288;
    p1[0] = 0.09834777464010795;
    p1[1] = 0.04137361348961467;
    p1[2] = -0.7341691126967387;
    p2[0] = -0.03081373001231996;
    p2[1] = 0.2323470126244768;
    p2[2] = 0.426387557408945;
    pf[0] = -0.3728087417235042;
    pf[1] = -0.2364545837571863;
    pf[2] = 2.023690886603053;
    F_fac_wkspc<RHR, 2> w;
    set_args<RHR, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, pf, sf);
    set_lambda<RHR, 3>(w, p0, p1, p2, pf, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.2647221281773937, 2.22045e-16);
    ASSERT_NEAR(df[0], 0.00833384190278625, 3e-16);
    ASSERT_NEAR(df[1], 0.2304578531655662, 3e-16);
    ASSERT_NEAR(d2f[0], 0.6527748823302069, 8e-16);
    ASSERT_NEAR(d2f[1], 0.5773439119170971, 6e-16);
    ASSERT_NEAR(d2f[2], 0.5794797038336468, 4e-16);
  }
}

TEST (cost_funcs, mp0_fac_works) {
  double u0, u1, u2, h, s, s0, s1, s2, sf, lam[2], p0[3], p1[3], p2[3], pf[3],
    f, df[2], d2f[3];
  {
    u0 = 0.8175470920792863;
    u1 = 0.7224395923668423;
    u2 = 0.1498654424779668;
    h = 0.6596052529083072;
    s = 0.5185949425105382;
    s0 = 0.9729745547638625;
    s1 = 0.6489914927123561;
    s2 = 0.8003305753524015;
    sf = 0.4537977087269195;
    lam[0] = 0.5112413543858396;
    lam[1] = 0.4887586456141604;
    p0[0] = -0.1472014561512673;
    p0[1] = 1.007773405305439;
    p0[2] = -2.12365546241575;
    p1[0] = -0.5045864055140099;
    p1[1] = -1.27059444980866;
    p1[2] = -0.3825848027076484;
    p2[0] = 0.6486792620486207;
    p2[1] = 0.8257271492417582;
    p2[2] = -1.014943642680137;
    pf[0] = -0.4710699126831666;
    pf[1] = 0.1370248741300503;
    pf[2] = -0.2918633757535734;
    F_fac_wkspc<MP0, 2> w;
    set_args<MP0, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, pf, sf);
    set_lambda<MP0, 3>(w, p0, p1, p2, pf, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.5587087401852298, 2.22045e-16);
    ASSERT_NEAR(df[0], -0.3024733586635268, 4e-16);
    ASSERT_NEAR(df[1], -0.8839234895192173, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 7.734841306686919, 2e-15);
    ASSERT_NEAR(d2f[1], 1.537379799160529, 2e-15);
    ASSERT_NEAR(d2f[2], 1.365603642828094, 7e-16);
  }
}

TEST (cost_funcs, mp1_fac_works) {
  double u0, u1, u2, h, s, s0, s1, s2, sf, lam[2], p0[3], p1[3], p2[3], pf[3],
    f, df[2], d2f[3];
  {
    u0 = 0.2919840799617149;
    u1 = 0.4316511702487202;
    u2 = 0.01548712563601895;
    h = 0.9840637243791538;
    s = 0.167168409914656;
    s0 = 0.1062163449286638;
    s1 = 0.372409740055537;
    s2 = 0.1981184025429746;
    sf = 0.4896876380160239;
    lam[0] = 0.5127098684891452;
    lam[1] = 0.4872901315108547;
    p0[0] = -0.5606649965281846;
    p0[1] = 2.177778709184154;
    p0[2] = 1.138465387329595;
    p1[0] = -2.496886503184502;
    p1[1] = 0.4413269317276085;
    p1[2] = -1.398137875811074;
    p2[0] = -0.2550551798808058;
    p2[1] = 0.1644040733187286;
    p2[2] = 0.7477340288081226;
    pf[0] = -0.2730469494029069;
    pf[1] = 1.576300146546075;
    pf[2] = -0.4809371517788452;
    F_fac_wkspc<MP1, 2> w;
    set_args<MP1, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, pf, sf);
    set_lambda<MP1, 3>(w, p0, p1, p2, pf, lam);
    eval(w, f);
    grad(w, df);
    hess(w, d2f);
    ASSERT_NEAR(f, 0.2852053655607227, 2.22045e-16);
    ASSERT_NEAR(df[0], 1.511262232883143, 2.22045e-16);
    ASSERT_NEAR(df[1], 0.2055717498377105, 2.22045e-16);
    ASSERT_NEAR(d2f[0], 4.005379755768571, 2e-15);
    ASSERT_NEAR(d2f[1], 1.042181708635259, 2.22045e-16);
    ASSERT_NEAR(d2f[2], 1.298770103749939, 5e-16);
  }
}

// TODO: test lag_mults
