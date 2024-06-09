package com.tntntn.axone;

public class M {
    public static final double qpi = 0.7853981633974483;
    public static final double hpi = 1.5707963267948966;
    public static final double e = 2.718281828459045;
    public static final double pi = 3.141592653589793;
    public static final double tau = 6.283185307179586;
    private static final double[] zeta_taylor = {-0.5, -0.91893853320467274178, -1.0031782279542924256, -1.000785194477042408, -0.99987929950057116496, -1.000001940896320456, -1.0000013011460139596, -0.99999983138417361078, -1.0000000057646759799, -1.0000000009110164892, -0.9999999998502992406, -1.0000000000094068957, -1.0000000000000409258, -0.999999999999934601, -1.000000000000006544, -0.9999999999999996988};
    public static double floor(double x) {
        return x < 0 ? (int) x - 1 : (int) x;
    }
    public static double round(double x) {
        return floor(x + 0.5);
    }
    public static double ceil(double x) {
        return -floor(-x);
    }
    public static double mod(double a, double b) {
        return b == 0 ? 0 : a - b * floor(a / b);
    }
    public static double abs(double x) {
        return x >= 0 ? x : -x;
    }
    public static double sgn(double x) {
        double ans = Double.NaN;
        if (x < 0) {
            ans = -1;
        } else if (x > 0) {
            ans = 1;
        } else if (x == 0) {
            ans = 0;
        }
        return ans;
    }
    private static double intpow(double x, long n) {
        double ans = 1;
        long b = (long) abs(n);
        while (b > 0) {
            if (b % 2 == 1) {
                ans *= x;
            }
            x *= x;
            b >>= 1;
        }
        return sgn(n) == -1 ? 1 / ans : ans;
    }
    public static double fact(double n) {
        double ans = 1;
        for (int i = 1; i <= n; i++) {
            ans *= i;
        }
        return ans;
    }
    public static double exp(double x) {
        double ans = 0;
        for (int i = 0; i <= 23; i++) {
            ans += intpow(mod(x, 1), i) / fact(i);
        }
        return intpow(e, (int) floor(x)) * ans;
    }
    private static double intln(double x) {
        double ans = 0;
        double m = 1;
        x = 1.6487212707001281 * x;
        if (x > 1) {
            while (true) {
                x /= e;
                ans++;
                if (x <= 1) {
                    break;
                }
            }
            ans--;
        } else if (x > 0) {
            x = 1 / x;
            m = -1;
            while (true) {
                x /= e;
                ans++;
                if (x <= 1) {
                    break;
                }
            }
        } else {
            ans = Double.NaN;
        }
        return m * ans;
    }
    public static double ln(double x) {
        double ans = Double.NaN;
        double x0 = x;
        if (x == 0) {
            ans = Double.NEGATIVE_INFINITY;
        } else if (x > 0) {
            x = x / exp(intln(x));
            ans = 0;
            for (int i = 1; i <= 25; i += 2) {
                ans += 2 * intpow((x - 1) / (x + 1), i) / i;
            }
        }
        return ans + intln(x0);
    }
    public static double log(double x, double base) {
        return ln(x) / ln(base);
    }
    public static double log10(double x) {
        return ln(x) / ln(10);
    }
    public static double log2(double x) {
        return ln(x) / ln(2);
    }
    public static double pow(double x, double y) {
        double ans;
        if (x == 0 && y == 0) {
            ans = 1;
        } else if (x == 0 && y != 0) {
            ans = 0;
        } else if (mod(y, 1) == 0) {
            ans = intpow(x, (long) y);
        } else {
            ans = exp(y * ln(x));
        }
        return ans;
    }
    public static double sqrt(double x) {
        return pow(x, 0.5);
    }
    public static double cbrt(double x) {
        return pow(x, 1 / 3d);
    }
    private static double sin1(double x) {
        double ans = 0;
        for (int i = 0; i <= 10; i++) {
            ans += intpow(-1, i) * intpow(x, 2 * i + 1) / fact(2 * i + 1);
        }
        return ans;
    }
    public static double sin(double x) {
        double ans = 0;
        x = mod(x, tau);
        if (x <= hpi) {
            ans = sin1(x);
        } else if (x <= pi) {
            ans = sin1(pi - x);
        } else if (x <= 3 * hpi) {
            ans = -sin1(x - pi);
        } else {
            ans = -sin1(tau - x);
        }
        return ans;
    }
    public static double cos(double x) {
        return sin(x + hpi);
    }
    public static double tan(double x) {
        return sin(x) / cos(x);
    }
    public static double cot(double x) {
        return cos(x) / sin(x);
    }
    public static double sec(double x) {
        return 1 / cos(x);
    }
    public static double csc(double x) {
        return 1 / sin(x);
    }
    public static double arcsin(double x) {
        double ans=0;
        double c=sqrt(pi);
        for (int i = 0; i <= 43; i++) {
            ans+=(pow(2,-i)*pow(1-abs(x),i)*c)/(fact(i)*(1+2*i)*sqrt(pi));
            c*=i+0.5;
        }
        return (sgn(x)*hpi)-sgn(x)*sqrt(2)*sqrt(1-abs(x))*ans;
    }
    public static double arccos(double x) {
        return hpi - arcsin(x);
    }
    public static double arctan(double x) {
        return arcsin(x / sqrt(x * x + 1));
    }
    public static double arccot(double x) {
        double ans=arctan(1 / x);
        return x < 0 ? ans + pi : ans;
    }
    public static double arcsec(double x) {
        return arccos(1 / x);
    }
    public static double arccsc(double x) {
        return arcsin(1 / x);
    }
    public static double sinh(double x) {
        return (exp(x) - exp(-x)) / 2;
    }
    public static double cosh(double x) {
        return (exp(x) + exp(-x)) / 2;
    }
    public static double tanh(double x) {
        return (exp(2 * x) - 1) / (exp(2 * x) + 1);
    }
    public static double coth(double x) {
        return (exp(2 * x) + 1) / (exp(2 * x) - 1);
    }
    public static double sech(double x) {
        return 2 / (exp(x) + exp(-x));
    }
    public static double csch(double x) {
        return 2 / (exp(x) - exp(-x));
    }
    public static double arcsinh(double x) {
        return ln(x + sqrt(x * x + 1));
    }
    public static double arccosh(double x) {
        return ln(x + sqrt(x * x - 1));
    }
    public static double arctanh(double x) {
        return ln((1 + x) / (1 - x)) / 2;
    }
    public static double arccoth(double x) {
        return ln((x + 1) / (x - 1)) / 2;
    }
    public static double arcsech(double x) {
        return ln((1 + sqrt(1 - x * x)) / x);
    }
    public static double arccsch(double x) {
        return ln((1 + sgn(x) * sqrt(1 + x * x)) / x);
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
    public static double toDegrees(double radians) {
        return radians * 57.29577951308232;
    }
    public static double toRadians(double degrees) {
        return degrees / 57.29577951308232;
    }
    public static double hypot(double... v) {
        double ans = 0;
        for (double value : v) {
            ans += intpow(value, 2);
        }
        return sqrt(ans);
    }
    public static double copySign(double x, double sign) {
        return abs(x) * sgn(sign);
    }
    public static double atan2(double y, double x) {
        double ans = 0;
        if (x > 0) {
            ans = arctan(y / x);
        } else if (x < 0 && y >= 0) {
            ans = arctan(y / x) + pi;
        } else if (x < 0 && y < 0) {
            ans = arctan(y / x) - pi;
        } else if (x == 0 && y > 0) {
            ans = hpi;
        } else if (x == 0 && y < 0) {
            ans = -hpi;
        } else if (x == 0 && y == 0) {
            ans = 0;
        }
        return ans;
    }
    public static double W0(double x) {
        double min = -1;
        double max = 703.2270331047702;
        double mid = 352.1135165523851;
        double ans;
        if (x == Double.POSITIVE_INFINITY) {
            ans = Double.POSITIVE_INFINITY;
        } else if (x < -1 / e) {
            ans = Double.NaN;
        } else if (Double.isNaN(x)) {
            ans = Double.NaN;
        } else {
            while (true) {
                if (mid * exp(mid) < x) {
                    min = mid;
                } else {
                    max = mid;
                }
                mid = (min + max) / 2.0;
                if (abs(max - min) < 2e-13) {
                    ans = mid;
                    break;
                }
            }
        }
        return ans;
    }
    private static double delta(double x) {
        return pow(2, max(ln(abs(x)), -1023) - 35);
    }
    public static double Gamma(double x) {
        double ans = 0;
        if (x <= 10) {
            ans = sqrt(tau * (x + 159)) * pow((x + 159) / e, x + 159) * exp(1 / (12 * (x + 159) + 0.55));
            for (int i = 0; i <= 159; i++) {
                ans /= x + i;
            }
        } else {
            ans = sqrt(tau * (x - 1)) * pow((x - 1) / e, x - 1) * exp(1 / (12 * x - 11.45));
        }
        return ans;
    }
    public static double DiGamma(double x) {
        return (ln(abs(Gamma(x + delta(x)))) - ln(abs(Gamma(x)))) / delta(x);
    }
    public static double invGamma(double x) {
        double min = 1.4616321;
        double max = 171.6243769;
        double mid = 86.5430045;
        double ans;
        if (x == Double.POSITIVE_INFINITY) {
            ans = Double.POSITIVE_INFINITY;
        } else if (x < 0.8856032) {
            ans = Double.NaN;
        } else if (Double.isNaN(x)) {
            ans = Double.NaN;
        } else {
            while (true) {
                if (Gamma(mid) < x) {
                    min = mid;
                } else {
                    max = mid;
                }
                mid = (min + max) / 2.0;
                if (abs(max - min) < 5e-14) {
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
            for (int j = 2; j <= sqrt(i); j++) {
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
            ans = ans * pow(2, x) * pow(pi, x - 1) * sin(hpi * x) * Gamma(1 - x);
        } else {
            for (int i = 0; i <= 15; i++) {
                ans += zeta_taylor[i] * pow(x, i);
            }
            ans += pow(x, 16) / (x - 1);
        }
        return ans;
    }
}
