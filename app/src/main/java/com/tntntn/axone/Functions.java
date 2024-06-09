package com.tntntn.axone;

import java.util.ArrayList;
import java.util.List;
// import org.apache.commons.math3.special;

// TO DO: Add javaDoc.
public class Functions {
    private static final double[] zeta_taylor = {-0.5, -0.91893853320467274178, -1.0031782279542924256, -1.000785194477042408, -0.99987929950057116496, -1.000001940896320456, -1.0000013011460139596, -0.99999983138417361078, -1.0000000057646759799, -1.0000000009110164892, -0.9999999998502992406, -1.0000000000094068957, -1.0000000000000409258, -0.999999999999934601, -1.000000000000006544, -0.9999999999999996988};
    private static final double[] tet_2_taylor = {1, 0.889364954620976, 0.00867654896536993, 0.0952388000751818, -0.00575234854012612, 0.0129665820200372, -0.00219604962303099, 0.00199674684791144, -0.000563354814878522, 0.000348242328188164, -0.000128532441264720, 0.0000670819244205308, -0.0000282987528227980, 0.0000138001319906329, -0.00000620190939837452, 0.00000295556146480966, -0.00000136867922453470, 0.000000649057075651896, -0.000000305166939328926, 0.000000144948206151230, -0.0000000687466431137918, 0.0000000327674451778934, -0.0000000156310874679970, 0.00000000747812838918081, -0.00000000358267681323948, 0.00000000171986774579521, -0.000000000826815962468010, 0.000000000398110468382690, -1.91942992587732e-10, 9.26631867862980e-11, -4.47869882913349e-11, 2.16712032119157e-11, -1.04969742914240e-11, 5.08944818189222e-12, -2.46987831946185e-12, 1.19965565647417e-12, -5.83165889022053e-13, 2.83702361814325e-13, -1.38118252499282e-13, 6.72883824735072e-14, -3.28030864119803e-14, 1.60015058245703e-14, -7.81025880698438e-15, 3.81431246327820e-15, -1.86381177420546e-15, 9.11196869462329e-16, -4.45694121263769e-16, 2.18105634006697e-16, -1.06780883356354e-16, 5.23008408468722e-17, -2.56274120196474e-17};
    private static final double[] tet_e_taylor = {1, 1.091767351258322138, 0.271483212901696469, 0.212453248176258214, 0.069540376139988952, 0.044291952090474256, 0.014736742096390039, 0.008668781817225539, 0.002796479398385586, 0.001610631290584341, 0.000489927231484419, 0.000288181071154065, 0.000080094612538551, 0.000050291141793809, 0.000012183790344901, 0.000008665533667382, 0.000001687782319318, 0.000001493253248573, 0.000000198760764204, 0.000000260867356004, 0.000000014709954143, 0.000000046834497327, 0.000000001549241666, 0.000000008741510781, 0.000000001125787310, 0.000000001707959267};
    private static final double[] ausin_taylor = {0, 2.29163807440958, 1.96043852439688, 1.07862851256147, 0.59622997993395, 0.28333997139829, 0.14193261194548, 0.06423734271234, 0.03026687705508, 0.01351721250427};
    private static final double[] slog_e_taylor = {0, 0.91594605649953, -0.20861842957759, -0.05450400630209, 0.07134941925273, -0.02004387374438, -0.01101258023037, 0.01207268318645, -0.0027292288076, -0.00269905319156, 0.00243941500632, -0.00036220360858, -0.00070125921262, 0.0005278215538, -0.00002987943551, -0.00018614540434, 0.00011722843042};
    private static final double[] sufact_coefficient = {2, 1, .798731835172434541585621072345730147, .577880975476483235803807592348110833, .393978809662971757177848639852917378, .257533958032332679820773329133486586, .162901958103705249541496101752195514, .100282419171352371943554511785342142, .0603184725913977494512136774562415014, .0355544582258061836048059212969418417, .0205859954874424134686332481358935023, .0117302279624549548734823541033644211, .00658835541777254650743317221091667507, .00365218351418374834372649788987162842, .00200039479760669665711545138631474960, .00108362752868222808502286098449166985, .000581036636299227699924018045799185045, .000308601963223618214714523083268563975, .000162, .000084, .00004};
    private static final double[] aufact_coefficient = {Double.NaN, 1, -0.7987318351724345, 0.6980641135593670, -0.6339640557572815, 0.5884152357911399, -0.5538887519936520, 0.5265479025985924, -0.5041914604280215, 0.4854529800293392, -0.4694346809094714};
    private static final double[] erf_taylor1 = {0.84270079294971476, 0.4151074974205948, -0.8302149948411895, 0.8302149948411895, 1.660429989682379, -8.302149948411897, 3.3208599793647573, 76.37977952538947, -192.60987880315588, -684.09715574914, 4449.9523723487755, 3413.8440587869713, -95826.73556454938, 116548.90183578526, 2066743.8498776173, -7163759.147485651, -4.354130950160196E7, 3.0199539342777455E8, 7.893311171957216E8, -1.1846505610935738E10, -4.722908997174056E9, 4.5961303120990784E11};
    private static final double[] erf_taylor2 = {0.99532226501895268, 0.020666985354092053, -0.08266794141636821, 0.2893377949572888, -0.8266794141636822, 1.570690886910996, 0.3306717656654717, -17.02959593177184, 64.15032253910175, -18.186947111601008, -953.6573721792242, 4141.994536725726, 2505.1692966817695, -101144.55699469196, 344454.16485841287, 1251941.8224282868, -1.4652483905749133E7, 2.105168095014646E7, 3.846727611834013E8, -2.2544481970385017E9, -4.830426614446579E9};
    private static final double[] erf_taylor3 = {0.99997790950300133, 1.392530519467479E-4, -8.355183116804872E-4, 0.004734603766189428, -0.02506554935041462, 0.12198567350535114, -0.5313896462287899, 1.968481142319228, -5.434211099169891, 5.046530602550132, 56.66819397141701, -430.8467146744075, 1451.7164086180885, 768.3292711286032, -39451.16943360529, 216730.4555522883, -195749.98917299512, -5327413.731531855, 3.8228482042711675E7};
    private static final double[] erf_taylor4 = {0.9999999845827421, 1.2698234671866558E-7, -1.0158587737493247E-6, 7.872905496557267E-6, -5.891980887746084E-5, 4.241210380403431E-4, -0.002921609833303058, 0.019131668286021034, -0.11799402828853155, 0.6761088703039579, -3.52096650981516, 15.997772413050026, -57.562849108096984, 108.55179977767766, 513.093980372932, -6927.098637203324, 41050.15764717917, -120588.30206135192};
    private static final double[] erf_taylor5 = {0.9999999999984626, 1.5670866531017336E-11, -1.567086653101734E-10, 1.535744920039699E-9, -1.4730614539156296E-8, 1.3809167587132477E-7, -1.2630718423999975E-6, 1.1249801665286726E-5, -9.734115454406728E-5, 8.159143221266587E-4, -0.006601684748561513, 0.05133038968733525, -0.38127020190212246, 2.683433445899848, -17.68384961334751, 107.06922654007936};
    private static final double[] erf_taylor6 = {0.99999999999999954, 4.9384851409642196E-15, -5.6792579121088527E-14, 6.432376896105896E-13, -7.170063114037426E-12};
    private static final double[] ci_values = {-0.1777840788066129, 0.3374039229009681, 0.4703563171953999, 0.422980828774865, 0.2858711963653835, 0.1196297860080003, -0.03212854851248112, -0.1409816978869304, -0.1934911221017388, -0.1900297496566439, -0.1420529475515193, -0.06805724389324713, 0.01110151951493011, 0.07669527848218452, 0.1156332032379343, 0.1224338825320096, 0.09943135857342192, 0.05534753133313361, 0.002678058835650657, -0.04545643300445537};
    private static final double[] ci_reciprocal_taylor = {-4,96,-4320,322560,-3628800,5748019200d,-1220496076800d,334764638208000d,115242726703104000d};

    private static final double MAX_SAFE_INT = 9007199254740992d;
    public static final double GOLDEN_RATIO = 1.618033988749895;
    public static final double OMEGA = 0.5671432904097838;

    public static final double DOTTIE = 0.73908513321516064;

    public static final double EULER_GAMMA = 0.5772156649015329;

    public static double mod(double a, double b) {
        return b == 0 ? 0 : a - b * floor(a / b);
    }

    public static double sgn(double x) {
        return x < 0 ? -1 : x > 0 ? 1 : 0;
    }

    public static double frac(double x) {
        return x % 1 + (x < 0 ? 1 : 0);
    }

    public static double floor(double x) {
        return x - frac(x);
    }

    public static double ceil(double x) {
        return -floor(-x);
    }

    public static double round(double x) {
        return floor(x + 0.5);
    }

    public static double fact(double n) {
        double ans = 1;
        for (int i = 1; i <= n; i++) {
            ans *= i;
        }
        return ans;
    }

    public static double max(double... v) {
        double ans = v[0];
        for (double i : v) {
            if (i > ans) {
                ans = i;
            }
        }
        return ans;
    }

    public static double min(double... v) {
        double ans = v[0];
        for (double i : v) {
            if (i < ans) {
                ans = i;
            }
        }
        return ans;
    }

    public static double W0(double x) {
        double min = x < 12 ? Math.log(x + 1) - 1 : Math.log(x) - Math.log(Math.log(x));
        double max = min + 1.1;
        double mid = (min + max) / 2.0;
        double ans;
        if (x == Double.POSITIVE_INFINITY) {
            ans = Double.POSITIVE_INFINITY;
        } else if (x < -1 / Math.E) {
            ans = Double.NaN;
        } else if (Double.isNaN(x)) {
            ans = Double.NaN;
        } else {
            while (true) {
                if (mid * Math.exp(mid) < x) {
                    min = mid;
                } else {
                    max = mid;
                }
                mid = (min + max) / 2.0;
                if (Math.abs(max - min) < 2e-13) {
                    ans = mid;
                    break;
                }
            }
        }
        return ans;
    }

    public static double delta(double x) {
        return Math.pow(2, Math.max(Math.log(Math.abs(x)), -1023) - 35);
    }

    public static double gamma(double x) {
        double ans;
        x += 20;
        if (x >= 0.5) {
            ans = Math.sqrt(Math.PI * 2 * x) * Math.pow(x / Math.E, x) * (1 + 1 / (12d * x) + 1 / (288d * x * x) + 139 / (51840d * x * x * x));
            for (int i = 0; i <= 20; i++) {
                ans /= x - i;
            }
        } else {
            x = 1 - x;
            ans = Math.sqrt(Math.PI * 2 * x) * Math.pow(x / Math.E, x) * (1 + 1 / (12d * x) + 1 / (288d * x * x) + 139 / (51840d * x * x * x));
            for (int i = 0; i <= 20; i++) {
                ans /= x - i;
            }
            ans = Math.PI / (ans * Math.sin(Math.PI * x));
        }
        return ans;
    }

    public static double diGamma(double x) {
        return x<=160?(Math.log(Math.abs(gamma(x + delta(x)))) - Math.log(Math.abs(gamma(x)))) / delta(x):(x>0?Math.log(x)-1/(2*x)-1/(12*x*x)+1/(120*Math.pow(x,4)):diGamma(1-x)-Math.PI*cot(Math.PI*x));
    }

    public static double invGamma(double x) {
        double a = Math.log((x + 0.1) / Math.sqrt(Math.PI * 2));
        double max = a / W0(a / Math.E) + 1;
        double min = max(max - 1, 1.46163221);
        double mid = (min + max) / 2.0;
        double ans;
        if (x == Double.POSITIVE_INFINITY) {
            ans = Double.POSITIVE_INFINITY;
        } else if (x < 0.8856032) {
            ans = Double.NaN;
        } else if (Double.isNaN(x)) {
            ans = Double.NaN;
        } else {
            while (true) {
                if (gamma(mid) < x) {
                    min = mid;
                } else {
                    max = mid;
                }
                mid = (min + max) / 2.0;
                if (Math.abs(max - min) < 5e-14) {
                    ans = mid;
                    break;
                }
            }
        }
        return ans;
    }

    public static double primePi(double x) {
        double ans = 0;
        double tmp = 0;
        for (int i = 2; i <= floor(x); i++) {
            for (int j = 2; j <= Math.sqrt(i); j++) {
                tmp += mod(i, j) > 0 ? 0 : 1;
            }
            ans += tmp == 0 ? 1 : 0;
            tmp = 0;
        }
        return ans;
    }

    public static double zeta(double x) {
        double ans = 0;
        if (x >= 5) {
            for (int i = 1; i < (x < 6 ? 100 : (x < 7 ? 30 : (x < 10 ? 15 : (x < 15 ? 6 : 3)))); i++) {
                ans += Math.pow(i, -x);
            }
        } else if (x <= -5) {
            for (int i = 1; i <= (1 - x < 6 ? 100 : (1 - x < 7 ? 30 : (1 - x < 10 ? 15 : (1 - x < 15 ? 6 : 3)))); i++) {
                ans += Math.pow(i, 1 + x);
            }
            ans = ans * Math.pow(2, x) * Math.pow(Math.PI, x - 1) * Math.sin(Math.PI / 2 * x) * gamma(1 - x);
        } else {
            for (int i = 0; i <= 15; i++) {
                ans += zeta_taylor[i] * Math.pow(x, i);
            }
            ans += Math.pow(x, 16) / (x - 1);
        }
        return ans;
    }

    public static double muMobius(double n) {
        n = round(n);
        double ans = Double.NaN;
        if (n >= 0) {
            ans = 0;
            for (int i = 1; i <= n; i++) {
                if (gcd(i, n) == 1) {
                    ans += Math.cos(2 * Math.PI * i / n);
                }
            }
        }
        return round(ans);
    }

    public static double primeZeta(double x) {
        double ans = Double.NaN;
        if (x > 1) {
            ans = 0;
            for (int i = 1; i <= 11; i++) {
                if (muMobius(i) != 0) {
                    ans += muMobius(i) * Math.log(zeta(i * x)) / i;
                }
            }
        } else if (x == 1 || x <= 0) {
            ans = Double.POSITIVE_INFINITY;
        }
        return ans;
    }

    public static double gcd(double a, double b) {
        a = round(Math.abs(a));
        b = round(Math.abs(b));
        double c;
        if (a >= MAX_SAFE_INT || b >= MAX_SAFE_INT) {
            if (a == 0) {
                return b;
            } else if (b == 0) {
                return a;
            } else if (a == 1 || b == 1) {
                return 1;
            } else {
                return Double.NaN;
            }
        }
        if (b == 0) {
            return a;
        }
        while (b > 0) {
            c = b;
            b = a % b;
            a = c;
        }
        return a;
    }

    public static double lcm(double a, double b) {
        a = round(Math.abs(a));
        b = round(Math.abs(b));
        return a * b / gcd(a, b);
    }

    public static double divisorSigma(double n, double z) {
        double ans = 0;
        n = round(n);
        for (int i = 1; i <= n; i++) {
            if (mod(n, i) == 0) {
                ans += Math.pow(i, z);
            }
        }
        return ans;
    }

    public static double isPrime(double n) {
        n = round(n);
        double ans = 0;
        if (n <= 1) {
            return 0;
        } else {
            for (int i = 1; i <= Math.sqrt(n); i++) {
                if (mod(n, i) == 0) {
                    ans++;
                }
            }
        }
        return ans == 1 ? 1 : 0;
    }

    public static double omegaPrime(double n) {
        n = round(n);
        double ans = 0;
        for (int i = 1; i <= n; i++) {
            if (mod(n, i) == 0) {
                ans += isPrime(i);
            }
        }
        return ans;
    }

    public static double OmegaPrime(double n) {
        n = round(n);
        double ans = 0;
        for (int i = 1; i <= n; i++) {
            if (mod(n, i) == 0 && omegaPrime(i) == 1) {
                ans++;
            }
        }
        return ans;
    }

    public static double lambdaLiouville(double n) {
        n = round(n);
        double ans = 0;
        for (int i = 1; i <= Math.sqrt(n); i++) {
            if (mod(n, i * i) == 0) {
                ans += muMobius(n / (i * i));
            }
        }
        return ans;
    }

    public static double phiEuler(double n) {
        n = round(n);
        double ans = 0;
        for (int i = 1; i <= n; i++) {
            if (gcd(i, n) == 1) {
                ans++;
            }
        }
        return ans;
    }

    public static double cot(double x) {
        return 1 / Math.tan(x);
    }

    public static double sec(double x) {
        return 1 / Math.cos(x);
    }

    public static double csc(double x) {
        return 1 / Math.sin(x);
    }

    public static double coth(double x) {
        return 1 / Math.tanh(x);
    }

    public static double sech(double x) {
        return 1 / Math.cosh(x);
    }

    public static double csch(double x) {
        return 1 / Math.sinh(x);
    }

    /**
     * Returns the arc cotangent of x value; the returned angle is in the
     * range -<i>pi</i>/2 through <i>pi</i>/2. Special case:
     * <ul><li>If the argument is NaN, then the result is NaN.
     * <li>If the argument is positive zero or negative zero, the result is a <i>pi</i> with the
     * same sign as the argument.
     * <li>If the argument is {@linkplain Double#isInfinite infinite},
     * the result is a zero with the same sign as the argument.
     * </ul>
     *
     * @param   x   the value whose arc cotangent is to be returned.
     * @return  the arc cotangent of the argument.
     */
    public static double acot(double x) {
        return Math.atan(1 / x);
    }

    /**
     * Returns the arc secant of x value; the returned angle is in the
     * range 0 through <i>pi</i>. Special case:
     * <ul><li>If the argument is NaN or the absolute value of the argument
     * less than 1, then the result is NaN.
     * <li>If the argument is <i>1</i>, the result is <i>0</i>.
     * <li>If the argument is -<i>1</i>, the result is <i>pi</i>.
     * <li>If the argument is {@linkplain Double#isInfinite infinite},
     * the result is <i>pi</i>/2.
     * </ul>
     *
     * @param   x   the value whose arc secant is to be returned.
     * @return  the arc secant of the argument.
     */
    public static double asec(double x) {
        return Math.acos(1 / x);
    }

    /**
     * Returns the arc cosecant of x value; the returned angle is in the
     * range -<i>pi</i>/2 through <i>pi</i>/2. Special case:
     * <ul><li>If the argument is NaN or the absolute value of the argument
     * less than 1, then the result is NaN.
     * <li>If the absolute value of the argument is <i>1</i>, then the result is
     * <i>pi</i>/2 with the same sign as the argument.
     * <li>If the argument is {@linkplain Double#isInfinite infinite},
     * the result is <i>0</i> with the same sign as the argument.
     * </ul>
     *
     * @param   x   the value whose arc cosecant is to be returned.
     * @return  the arc cosecant of the argument.
     */
    public static double acsc(double x) {
        return Math.asin(1 / x);
    }

    /**
     * Returns the arc hyperbolic sine of x value; Special case:
     * <ul><li>If the argument is NaN, then the result is NaN.
     * <li>If the argument is {@linkplain Double#isInfinite infinite},
     * the result is the same as the argument.
     * </ul>
     *
     * @param   x   the value whose arc hyperbolic sine is to be returned.
     * @return  the arc hyperbolic sine of the argument.
     */
    public static double asinh(double x) {
        return Math.log(x + Math.sqrt(x * x + 1));
    }

    /**
     * Returns the arc hyperbolic cosine of x value; Special case:
     * <ul><li>If the argument is NaN or the argument less than 1, then the result is NaN.
     * <li>If the argument is positive infinity, the result is the same as the argument.
     * <li>If the argument is one, the result is positive zero.
     * </ul>
     *
     * @param   x   the value whose arc hyperbolic cosine is to be returned.
     * @return  the arc hyperbolic cosine of the argument.
     */
    public static double acosh(double x) {
        return Math.log(x + Math.sqrt(x * x - 1));
    }

    /**
     * Returns the arc hyperbolic tangent of x value; Special case:
     * <ul><li>If the argument is NaN or the absolute value of the argument
     * larger than 1, then the result is NaN.
     * <li>If the argument is 1 or -1, the result is infinity
     * with the same sign as the argument.
     * </ul>
     *
     * @param   x   the value whose arc hyperbolic tangent is to be returned.
     * @return  the arc hyperbolic tangent of the argument.
     */
    public static double atanh(double x) {
        return Math.log((1 + x) / (1 - x)) / 2;
    }

    /**
     * Returns the arc hyperbolic cotangent of x value; Special case:
     * <ul><li>If the argument is NaN or the absolute value of the argument
     * less than 1, then the result is NaN.
     * <li>If the argument is 1 or -1, the result is infinity
     * with the same sign as the argument.
     * <li> If the argument is {@linkplain Double#isInfinite infinite},
     * the result is zero with the same sign as the argument.
     * </ul>
     *
     * @param   x   the value whose arc hyperbolic cotangent is to be returned.
     * @return  the arc hyperbolic cotangent of the argument.
     */
    public static double acoth(double x) {
        return Math.log((x + 1) / (x - 1)) / 2;
    }

    /**
     * Returns the arc hyperbolic secant of x value; Special case:
     * <ul><li>If the argument is NaN or the argument is negative or
     * the argument less lagrer 1, then the result is NaN.
     * <li>If the argument is positive zero the result is infinity.
     * <li>If the argument is 1, the result is positive zero.
     * </ul>
     *
     * @param   x   the value whose arc hyperbolic secant is to be returned.
     * @return  the arc hyperbolic secant of the argument.
     */
    public static double asech(double x) {
        return Math.log((1 + Math.sqrt(1 - x * x)) / x);
    }
    /**
     * Returns the arc hyperbolic cosecant of x value; Special case:
     * <ul><li>If the argument is NaN then the result is NaN.
     * <li>If the argument is zero the result is infinity
     * with the same sign as the argument.
     * <li> If the argument is {@linkplain Double#isInfinite infinite},
     * the result is zero with the same sign as the argument.
     * </ul>
     *
     * @param   x   the value whose arc hyperbolic cosecant is to be returned.
     * @return  the arc hyperbolic cosecant of the argument.
     */
    public static double acsch(double x) {
        return Math.log((1 + sgn(x) * Math.sqrt(1 + x * x)) / x);
    }

    public static double sinp(double x) {
        return 3 - 2 * Math.cosh(2 / 3d * asinh((3 * x - 4) / 2d));
    }

    public static double cosp(double x) {
        return -2 * Math.sinh(1 / 3d * asinh((3 * x - 4) / 2d));
    }

    public static double tanp(double x) {
        return sinp(x) / cosp(x);
    }

    public static double cotp(double x) {
        return cosp(x) / sinp(x);
    }

    public static double secp(double x) {
        return 1 / cosp(x);
    }

    public static double cscp(double x) {
        return 1 / sinp(x);
    }

    public static double asinp(double x) {
        return 1 / 3d * (4 + 2 * Math.sinh(3 / 2d * acosh((3 - x) / 2d)));
    }

    public static double acosp(double x) {
        return 2 / 3d * (2 - Math.sinh(3 * asinh(x / 2d)));
    }

    public static double atanp(double x) {
        return 1 / 3d * (4 + 2 * Math.sinh(3 * acoth(Math.sqrt(2 * x * x + 5 - 2 * x * Math.sqrt(x * x + 4)))));
    }

    public static double acotp(double x) {
        return x == 0 ? 4 / 3d : atanp(1 / x);
    }

    public static double asecp(double x) {
        return acosp(1 / x);
    }

    public static double acscp(double x) {
        return asinp(1 / x);
    }

    public static double[] collatzSeq(double x) {
        long n = (long) round(x);
        ArrayList S = new ArrayList();
        S.add(n);
        while (n > 1) {
            if (n % 2 == 0) {
                n = n / 2;
            } else {
                n = 3 * n + 1;
            }
            S.add(n);
        }
        Object[] a = S.toArray();
        double[] ans = new double[S.size()];
        for (int i = 0; i < S.size(); i++) {
            ans[i] = Double.parseDouble(a[i].toString());
        }
        return ans;
    }

    public static double collatzLen(double x) {
        long n = (long) round(x);
        long k = 0;
        while (n != 1 && n != 0 && n != -1 && n != -5 && n != -17) {
            if (n % 2 == 0) {
                n = n / 2;
            } else {
                n = 3 * n + 1;
            }
            k++;
        }
        return k + 1;
    }

    public static double collatzMax(double x) {
        long n = (long) round(x);
        long m = n;
        while (n != 1 && n != 0 && n != -1 && n != -5 && n != -17) {
            if (n % 2 == 0) {
                n = n / 2;
            } else {
                n = 3 * n + 1;
            }
            m = n > 0 ? (long) max(m, n) : (long) min(m, n);
        }
        return m;
    }

    public static double logb(double b, double x) {
        return Math.log(x) / Math.log(b);
    }

    public static double factorial(double x) {
        return gamma(x + 1);
    }

    public static double nCr(double a, double b) {
        return factorial(a) / (factorial(b) * factorial(a - b));
    }

    public static double nPr(double a, double b) {
        return factorial(a) / factorial(a - b);
    }

    public static double mean(List<Double> in) {
        double ans = 0;
        for (int i = 0; i < in.size(); i++) {
            ans += in.get(i);
        }
        return ans / in.size();
    }

    public static double max(List<Double> in) {
        if (in.isEmpty()) {
            return Double.NaN;
        }
        double ans = in.get(0);
        for (int i = 0; i < in.size(); i++) {
            ans = max(in.get(i), ans);
        }
        return ans;
    }

    public static double min(List<Double> in) {
        if (in.isEmpty()) {
            return Double.NaN;
        }
        double ans = in.get(0);
        for (int i = 0; i < in.size(); i++) {
            ans = min(in.get(i), ans);
        }
        return ans;
    }

    public static double lcm(List<Double> in) {
        if (in.isEmpty()) {
            return Double.NaN;
        }
        double ans = in.get(0);
        for (int i = 0; i < in.size(); i++) {
            ans = lcm(in.get(i), ans);
        }
        return ans;
    }

    public static double gcd(List<Double> in) {
        if (in.isEmpty()) {
            return Double.NaN;
        }
        double ans = in.get(0);
        for (int i = 0; i < in.size(); i++) {
            ans = gcd(in.get(i), ans);
        }
        return ans;
    }

    public static double rect(Double x) {
        return Math.abs(x) < 0.5 ? 1 : Math.abs(x) > 0.5 ? 0 : 0.5;
    }

    public static double tri(Double x) {
        return max(1 - Math.abs(x), 0);
    }

    public static double random(List<Double> in) {
        if (in.size() == 1) {
            return Math.random() * in.get(0);
        }
        if (in.size() == 2) {
            return Math.random() * (in.get(1) - in.get(0)) + in.get(0);
        }
        return Math.random();
    }

    public static double SuSun(double x) {
        x += 101.3869252578924;
        double ans = Math.sqrt(3 / x) * (1 - (3 * Math.log(x) / (10 * x)));
        for (int i = 0; i < 100; i++) {
            ans = Math.asin(ans);
        }
        return ans;
    }

    public static double tetBase2(double x) {
        if (x > 4.8) {
            return Double.POSITIVE_INFINITY;
        } else if (x == -2) {
            return Double.NEGATIVE_INFINITY;
        } else if (x < -2) {
            return Double.NaN;
        } else if (-0.5 <= x && x <= 0.5) {
            double ans = 0;
            for (int i = 0; i <= 50; i++) {
                ans += tet_2_taylor[i] * Math.pow(x, i);
            }
            return ans;
        } else if (x > 0.5) {
            return Math.pow(2, tetBase2(x - 1));
        } else {
            return Math.log(tetBase2(x + 1)) / Math.log(2);
        }
    }

    public static double compositeZeta(double x) {
        return zeta(x) - 1 - primeZeta(x);
    }

    public static double xiRiemann(double x) {
        return x * (x - 1) * Math.pow(Math.PI, -x / 2) * gamma(x / 2) * zeta(x) / 2;
    }

    public static double thetaChebyshev(double x) {
        x = round(x);
        double ans = 0;
        for (int i = 1; i <= x; i++) {
            if (isPrime(i) == 1) {
                ans += Math.log(i);
            }
        }
        return ans;
    }

    public static double lambdaMangoldt(double x) {
        x = round(x);
        double ans = 0;
        for (int i = 1; i <= x; i++) {
            if (mod(x, i) == 0) {
                ans += muMobius(i) * Math.log(i);
            }
        }
        return -ans;
    }

    public static double psiChebyshev(double x) {
        x = round(x);
        double ans = 0;
        for (int i = 1; i <= x; i++) {
            ans += lambdaMangoldt(i);
        }
        return ans;
    }

    public static double beta(double x, double y) {
        return gamma(x) * gamma(y) / gamma(x + y);
    }

    public static double fibonacci(double x) {
        return (Math.pow(GOLDEN_RATIO, x - 2) / Math.sqrt(5)) * (1 + GOLDEN_RATIO) - (Math.cos(Math.PI * x) / Math.pow(GOLDEN_RATIO, x - 2)) / (Math.sqrt(5) * (1 + GOLDEN_RATIO));
    }

    public static double ssqrt(double x) {
        return Math.log(x) / (W0(Math.log(x)));
    }

    public static double trunc(double x) {
        return (int) x;
    }

    public static double invFactorial(double x) {
        return invGamma(x) - 1;
    }

    public static double Fermi(double x) {
        return sgn(x) * Math.pow(10, round(Math.log10(Math.abs(x))));
    }

    public static double Fermi(double x, double base) {
        return sgn(x) * Math.pow(base, round(logb(base, Math.abs(x))));
    }

    public static double dms(double x) {
        return sgn(x) * (floor(Math.abs(x)) + floor(60d * (Math.abs(x) % 1d)) / 100d + 0.36d * mod(Math.abs(x), 1 / 60d));
    }

    public static double deg(double x) {
        return sgn(x) * (floor(Math.abs(x)) + floor(100d * (Math.abs(x) % 1d)) / 60d + 25d * mod(Math.abs(x), 1 / 60d) / 9d);
    }

    //public static double arg(double x) {
    //    return x < 0 ? Math.PI : 0;
    //}

    public static double bernoulli(double x) {
        return x == 0 ? 1 : -x * zeta(1 - x);
    }

    public static double digits(double x) {
        return Math.max(floor(Math.log10(Math.abs(x))), 1);
    }

    public static double digits(double x, double base) {
        return Math.max(floor(logb(base, Math.abs(x))), 1);
    }

    public static double concat(double x, double y) {
        return x * Math.pow(10, digits(y)) + y;
    }

    public static double concat(double x, double y, double base) {
        return x * Math.pow(base, digits(y, base)) + y;
    }

    public static double xexp(double x) {
        return x * Math.exp(x);
    }

    public static double tetHeightInf(double x) {
        if (x == -1) {
            return -1;
        } else if (x == 1) {
            return 1;
        } else if (x == Math.exp(1 / Math.E)) {
            return Math.E;
        } else if (x > Math.exp(1 / Math.E)) {
            return Double.POSITIVE_INFINITY;
        } else if (x >= Math.exp(-Math.E)) {
            return W0(-Math.log(x)) / -Math.log(x);
        } else {
            return Double.NaN;
        }
    }

    public static double harmonic(double x) {
        return diGamma(x + 1) + EULER_GAMMA;
    }

    public static double clamp(double x, double min_val, double max_val) {
        return Math.max(Math.min(min_val, max_val), Math.min(Math.max(min_val, max_val), x));
    }

    //public static double digitsum(double x) {
    //    double ans = mod(Math.abs(x), 10);
    //    for (int i = 1; i <= digits(x); i++) {
    //        ans += floor(mod(Math.abs(x), Math.pow(10, i + 1)) / Math.pow(10, i));
    //    }
    //    return ans;
    //}

    //public static double digitsum(double x, double base) {
    //    double ans = mod(Math.abs(x), base);
    //    for (int i = 1; i <= digits(x, base); i++) {
    //        ans += floor(mod(Math.abs(x), Math.pow(base, i + 1)) / Math.pow(base, i));
    //    }
    //    return ans;
    //}

    public static double AuSin(double x) {
        if (x == 0 || x == Math.PI) {
            return Double.POSITIVE_INFINITY;
        } else if (x < Math.PI && x > 0) {
            if (Math.PI - 0.627551832049 >= x && x >= 0.627551832049) {
                double ans = 0;
                for (int i = 1; i <= 9; i++) {
                    ans += ausin_taylor[i] * Math.pow(x - Math.PI / 2, 2 * i);
                }
                return ans;
            } else {
                for (int i = 1; i <= 100; i++) {
                    x = Math.sin(x);
                }
                return 3 / (x * x) + 5 * Math.log(x) / 6d - 102.73788774542039;
            }
        } else {
            return Double.NaN;
        }
    }

    public static double sinIterated(double x, double iterations) {
        return SuSun(iterations + AuSin(x));
    }

    public static double tetBaseE(double x) {
        if (x < -2) {
            return Double.NaN;
        } else if (x == -2) {
            return Double.NEGATIVE_INFINITY;
        } else if (x > 3.7) {
            return Double.POSITIVE_INFINITY;
        } else if (-0.5 <= x && x <= 0.5) {
            double ans = 0;
            for (int i = 0; i <= 25; i++) {
                ans += tet_e_taylor[i] * Math.pow(x, i);
            }
            return ans;
        } else if (x > 0.5) {
            return Math.exp(tetBaseE(x - 1));
        } else {
            return Math.log(tetBaseE(x + 1));
        }
    }

    public static double penBaseE(double x) {
        double ans = Math.exp(1.86573322821 * (x - 10.7518255)) - 1.8503545290271812;
        for (int i = 1; i <= 13; i++) {
            ans = tetBaseE(ans);
        }
        return ans;
    }

    public static double slogBaseE(double x) {
        if (x == Double.POSITIVE_INFINITY) {
            return Double.POSITIVE_INFINITY;
        } else if (x == Double.NEGATIVE_INFINITY) {
            return -2;
        } else if (Double.isNaN(x)) {
            return Double.NaN;
        } else if (0.4983 <= x && x <= 1.646) {
            double ans = 0;
            for (int i = 1; i <= 16; i++) {
                ans += slog_e_taylor[i] * Math.pow(x - 1, i);
            }
            return ans;
        } else if (x < 0.4983) {
            return slogBaseE(Math.exp(x)) - 1;
        } else {
            return slogBaseE(Math.log(x)) + 1;
        }
    }

    public static double expIterated(double x, double iterations) {
        return tetBaseE(iterations + slogBaseE(x));
    }

    public static double penBase2(double x) {
        double ans = Math.pow(5.171798, x - 8.9249685988278) - 1.7439091761327365;
        for (int i = 1; i <= 11; i++) {
            ans = tetBase2(ans);
        }
        return ans;
    }

    //public static double deltaLn(double x, double y) {
    //    if (y > 0 && y % 1 == 0) {
    //        return Math.pow(-1, 1 + y) * Math.pow(x, -y) * gamma(y);
    //    } else {
    //        return -(Math.pow(x, -y) * (EULER_GAMMA - Math.log(x) + diGamma(1 - y))) / gamma(1 - y);
    //    }
    //}

    public static double SuFact(double x) {
        double k = Math.exp(0.6127874523307084 * (x - 2.91938596545218));
        double ans = 0;
        for (int i = 20; i >= 0; i--) {
            ans += sufact_coefficient[i];
            ans *= k;
        }
        return factorial(factorial(ans / k));
    }

    public static double AuFact(double x) {
        double k = invFactorial(invFactorial(invFactorial(invFactorial(x)))) - 2;
        double ans = 0;
        for (int i = 1; i <= 10; i++) {
            ans += aufact_coefficient[i] * Math.pow(k, i);
        }
        return logb(3 - 2 * EULER_GAMMA, ans) + 4.9200636447;
    }

    public static double factIterated(double x, double iterations) {
        return SuFact(iterations + AuSin(x));
    }

    public static double arcPenBaseE(double x) {
        for (int i = 1; i <= 13; i++) {
            x = slogBaseE(x);
        }
        return Math.log(x + 1.8503545290271812) / 1.86573322821 + 10.7518255;
    }

    public static double tetEIterated(double x, double iterations) {
        return penBaseE(iterations + arcPenBaseE(x));
    }

    /**
     * Returns the sine of x value divided by x value; Special case:
     * <ul><li>If the argument is NaN, then the result is NaN.
     * <li>If the argument is {@linkplain Double#isInfinite infinite},
     * the result is positive zero.
     * <li>If the argument is 0, the result is 1.
     * </ul>
     *
     * @param   x   the value whose sinс is to be returned.
     * @return  the arc hyperbolic cosine of the argument.
     */
    public static double sinc(double x) {
        return x == 0 ? 1 : Double.isInfinite(x) ? 0 : Math.sin(x) / x;
    }

    /**
     * Returns the cosine of x value divided by x value; Special case:
     * <ul><li>If the argument is NaN, then the result is NaN.
     * <li>If the argument is {@linkplain Double#isInfinite infinite},
     * the result is zero with the same sign as the argument.
     * <li>If the argument is zero, the result is infinity
     * with the same sign as the argument.
     * </ul>
     *
     * @param   x   the value whose cosс is to be returned.
     * @return  the arc hyperbolic cosine of the argument.
     */
    public static double cosc(double x) {
        return Double.isInfinite(x) ? 0*sgn(x) : Math.cos(x) / x;
    }

    /**
     * Returns the tangent of x value divided by x value; Special case:
     * <ul><li>If the argument is NaN or the argument is
     * {@linkplain Double#isInfinite infinite}, then the result is NaN.
     * <li>If the argument is x+<i>pi</i>/2 is multiple of <i>pi</i>,
     * but the argument is not equal to 0, then the result is infinity
     * with the same sign as the argument.
     * <li>If the argument is zero, the result is one.
     * </ul>
     *
     * @param   x   the value whose tanс is to be returned.
     * @return  the arc hyperbolic tangent of the argument.
     */
    public static double tanc(double x) {
        return x == 0 ? 1 : Math.tan(x) / x;
    }

    public static double cotc(double x) {
        return cot(x) / x;
    }

    public static double secc(double x) {
        return sec(x) / x;
    }

    public static double cscc(double x) {
        return csc(x) / x;
    }

    public static double sinhc(double x) {
        return x == 0 ? 1 : Math.sinh(x) / x;
    }

    public static double coshc(double x) {
        return Math.cosh(x) / x;
    }

    public static double tanhc(double x) {
        return x == 0 ? 1 : Math.tanh(x) / x;
    }

    public static double cothc(double x) {
        return coth(x) / x;
    }

    public static double sechc(double x) {
        return sech(x) / x;
    }

    public static double cschc(double x) {
        return csch(x) / x;
    }

    public static double sinpc(double x) {
        return x == 0 ? 1 : sinp(x) / x;
    }

    public static double cospc(double x) {
        return cosp(x) / x;
    }

    public static double tanpc(double x) {
        return x == 0 ? 1 : tanp(x) / x;
    }

    public static double cotpc(double x) {
        return cotp(x) / x;
    }

    public static double secpc(double x) {
        return secp(x) / x;
    }

    public static double cscpc(double x) {
        return cscp(x) / x;
    }

    public static double erf(double x) {
        double a = sgn(x);
        x = Math.abs(x);
        double ans = 0;
        if (x < 0.5) {
            for (int i = 0; i <= 10; i++) {
                ans += (1 - 2 * (i % 2)) * Math.pow(x, 2 * i + 1) / (fact(i) * (2 * i + 1));
            }
            ans *= 1.1283791670955126;
        } else if (x < 1.5) {
            for (int i = 0; i <= 21; i++) {
                ans += erf_taylor1[i] * Math.pow(x - 1, i) / fact(i);
            }
        } else if (x < 2.5) {
            for (int i = 0; i <= 20; i++) {
                ans += erf_taylor2[i] * Math.pow(x - 2, i) / fact(i);
            }
        } else if (x < 3.5) {
            for (int i = 0; i <= 18; i++) {
                ans += erf_taylor3[i] * Math.pow(x - 3, i) / fact(i);
            }
        } else if (x < 4.5) {
            for (int i = 0; i <= 17; i++) {
                ans += erf_taylor4[i] * Math.pow(x - 4, i) / fact(i);
            }
        } else if (x < 5.5) {
            for (int i = 0; i <= 15; i++) {
                ans += erf_taylor5[i] * Math.pow(x - 5, i) / fact(i);
            }
        } else if (x < 6) {
            for (int i = 0; i <= 4; i++) {
                ans += erf_taylor6[i] * Math.pow(x - 5.75, i) / fact(i);
            }
        } else {
            ans = 1;
        }
        return ans * a;
    }

    public static double erfc(double x) {
        return 1 - erf(x);
    }

    public static double tetee1(double x) {
        double base = 1.444667861009766;
        x += 5.7982472438146786;
        double y = -Math.log(x);
        double ans = 1 + y / (3 * x);
        ans += (y * y + y + 0.5) / (9 * x * x);
        ans += (y * y * y + 2.5 * y * y + 2.5 * y + 0.7) / (27 * x * x * x);
        ans += (y * y * y * y + 13 / 3d * y * y * y + 45 / 6d * y * y + 5.3 * y + 67 / 60d) / (81 * x * x * x * x);
        ans += (y * y * y * y * y + 77 / 12d * y * y * y * y + 101 / 6d * y * y * y + 20.75 * y * y + 653 / 60d * y + 2701 / 1680d) / (243 * x * x * x * x * x);
        ans = Math.E * (1 - 2 / x * ans);
        return logb(base, logb(base, logb(base, ans)));
    }

    public static double penee1(double x) {
        double ans = Math.pow(4.804889, x - 12.224969336229359) - 1.5610637239194365;
        for (int i = 1; i <= 14; i++) {
            ans = tetee1(ans);
        }
        return ans;
    }

    public static double SuXexp(double x) {
        if (x > 4) {
            return Double.POSITIVE_INFINITY;
        } else {
            x -= 9.128694199839897;
            double y = Math.log(-x);
            double ans = -1 / x + 0.5 * y / (x * x);
            ans += (-0.25 * y * y + 0.25 * y - 0.2) / (x * x * x);
            ans += (0.125 * y * y * y - 0.3125 * y * y + 0.375 * y - 7 / 48d) / (x * x * x * x);
            for (int i = 1; i <= 8; i++) {
                ans = xexp(ans);
            }
            return ans;
        }
    }

    public static double AuXexp(double x) {
        double W6 = W0(W0(W0(W0(W0(W0(x))))));
        return 7.100549278723377 - 1 / W6 + 0.5 * Math.log(W6);
    }

    public static double xexpEIterated(double x, double iterations) {
        return SuXexp(iterations + AuXexp(x));
    }

    public static double Q_rsqrt(double number) {
        int i;
        float x2, y;
        final float threehalfs = 1.5F;

        x2 = (float) (number) * 0.5F;
        y = (float) (number);
        i = Float.floatToRawIntBits(y);
        i = 0x5f3759df - (i >> 1);
        y = Float.intBitsToFloat(i);
        y = y * (threehalfs - (x2 * y * y));
        // y  = y * ( threehalfs - ( x2 * y * y ) );

        return (double) y;
    }

    public static double Si(double x) {
        double k = sgn(x);
        x = Math.abs(x);
        if (x < 10) {
            double ans = 0;
            for (int n = 0; n <= (x < 1 ? 6 : x < 3 ? 10 : x < 5 ? 13 : 20); n++) {
                ans += Math.pow(-1, n) * Math.pow(x, 2 * n) / (Math.pow(1 + 2 * n, 2) * fact(2 * n));
            }
            return k * x * ans;
        } else if (x == Double.POSITIVE_INFINITY) {
            return Double.POSITIVE_INFINITY * k;
        } else if (Double.isNaN(x)) {
            return Double.NaN;
        } else {
            double ans1 = 0;
            for (int n = 0; n < (x > 20 ? 9 : 4); n++) {
                ans1 += Math.pow(-1, n) * fact(2 * n + 1) / Math.pow(x, 2 * n + 1);
            }
            ans1 *= sinc(x);
            double ans2 = 0;
            for (int n = 0; n < (x > 20 ? 9 : 4); n++) {
                ans2 += Math.pow(-1, n) * fact(2 * n) / Math.pow(x, 2 * n);
            }
            ans2 *= cosc(x);
            return k * (Math.PI / 2 - ans1 - ans2);
        }
    }

    public static double Ci(double x) {
        if (x < 0 || Double.isNaN(x)) {
            return Double.NaN;
        } else if (x == 0) {
            return Double.NEGATIVE_INFINITY;
        } else if (x == Double.POSITIVE_INFINITY) {
            return 0;
        } else if (x<1.5) {
            double ans=EULER_GAMMA+Math.log(x);
            for (int i = 2; i <= 18; i+=2) {
                ans+=Math.pow(x,i)/ci_reciprocal_taylor[i/2-1];
            }
            return ans;
        } else if (x < 10) {
            double a = max(0.5, round(2 * x) / 2);
            double ans = Math.log(x)+ci_values[(int)(2*a-1)]-Math.log(a);
            ans+=(x-a)*(Math.cos(a)-1)/a;
            ans+=(x-a)*(x-a)*(1-a*Math.sin(a)-Math.cos(a))/(2*a*a);
            ans+=Math.pow(x-a,3)*((2*a*Math.sin(a)-(a*a-2)*Math.cos(a))/(6*a*a*a)-1/(3*a*a*a));
            ans+=Math.pow(x-a,4)*(1/(4*a*a*a*a)+(a*(a*a-6)*Math.sin(a)+3*(a*a-1)*Math.cos(a))/(24*a*a*a*a));
            return ans;
        } else {
            double ans1 = 0;
            for (int n = 0; n <= (x>7?3:2); n++) {
                ans1+=Math.pow(-1,n)*fact(2*n)/Math.pow(x,2*n);
            }
            ans1*=sinc(x);
            double ans2 = 0;
            for (int n = 0; n <= (x>7?3:2); n++) {
                ans2+=Math.pow(-1,n)*fact(2*n+1)/Math.pow(x,2*n+1);
            }
            ans2*=cosc(x);
            return ans1-ans2;
        }
    }
}
// erfi
// arcerf
// arcerfc
// arcerfi
// W-1
// Li
// li
// Ei
// E1
// Ein
// Ci
// Cin
// Shi
// Chi
// Fresnel integral
// polylog
// stieljes
// inverf
// invquadratic
// gbarnes
// hyperfact
// gammaincomplete
// Gammaincomople
// gammagerner
// betaincomplete
// betareg
// tet
// slog2
// slog
// pen
// ape2
// apeee1
// ape
// tetit
// p_n
// zetaHurwitz
// polygamma
// Eulers

// javadoc