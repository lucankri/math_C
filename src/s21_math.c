#include <stdio.h>
#include <stdlib.h>
#include "s21_math.h"
#include <string.h>
#include <float.h>

int s21_abs(int x) {
    return x < 0 ? -x : x;
}

long double s21_fabs(double x) {
    return x < 0 ? -x : x;
}

long double s21_ceil(double x) {
  if (!is_fin(x)) {
    return x;
  }
  long double ceil_x = (long long int)x;
  if (s21_fabs(x) > 0. && x != ceil_x) {
    if (x != DBL_MAX) {
      if (x > 0.) {
        ceil_x += 1;
      }
    } else {
      return DBL_MAX;
    }
  }
  return ceil_x;
}

long double s21_floor(double x) {
  double arg = x;
  int num = arg;
  long double res;
  if (x == s21_INFINITY_POSITIVE || x == s21_INFINITY_NEGATIVE) {
    res = x;
  } else if (arg == num) {
    res = num;
  } else if (arg > (long double)num) {
    res = num;
  } else if (arg < (long double)num) {
    res = num - 1;
  }

  return x != x ? s21_NAN : res;
}

long double s21_fmod(double x, double y) {
    long double res;

    if (!is_fin(x) || is_nan(y)) {
        res = s21_NAN;
    } else if (is_inf(x) && is_inf(y)) {
        res = s21_NAN;
    } else if (is_inf(y)) {
        res = x;
    } else if (s21_fabs(y) < 1e-7) {
        res = s21_NAN;
    } else if (s21_fabs(x) < 1e-7) {
        res = 0;
    } else {
        long long int mod = 0;
        mod = x / y;
        res = (long double)x - mod * (long double)y;
    }
    return res;
}


long double s21_exp(double x) {
    long double res = 0;
    if (is_nan(x)) {
        res = s21_NAN;
    } else if (x == s21_INFINITY_POSITIVE) {
        res = s21_INFINITY_POSITIVE;
    } else if (x == s21_INFINITY_NEGATIVE) {
        res = 0;
    } else if (x == 0) {
        res = 1;
    } else {
        long double sum = 1, member = 1, eps = s21_E;
        int negative = 0;
        if (x < 0) {
            x *= -1;
            negative = 1;
        }
        for (int i = 1;; i++) {
            member *= x / i;
            sum += member;
            // printf("sum %.10Lf\n", sum);
            if (member <= eps || sum > DBL_MAX) {
                break;
            }
        }
        if (negative) {
            if (sum > DBL_MAX) {
                res = 0;
            } else {
                res = 1 / sum;
            }
        } else {
            if (sum > DBL_MAX) {
                res = s21_INFINITY_POSITIVE;
            } else {
                res = sum;
            }
        }
    }
    return res;
}



long double s21_pow(double base, double exp) {
    long double res;
    if (base == 1 || exp == 0) {
        // trivial cases
        res = 1;
    } else if (base != base) {
        // nan
        res = base;
    } else if (base == s21_INFINITY_POSITIVE) {
        if (exp == s21_INFINITY_NEGATIVE) {
            res = 0;
        } else if (exp != exp || exp == s21_INFINITY_POSITIVE) {
            res = exp;
        } else {
            // normal exp
            if (exp < 0) {
                res = 0;
            } else {
                res = base;
            }
        }
    } else if (base == s21_INFINITY_NEGATIVE) {
        if (exp == s21_INFINITY_NEGATIVE) {
            res = 0;
        } else if (exp != exp || exp == s21_INFINITY_POSITIVE) {
            res = exp;
        } else {
            // normal exp
            if (exp < 0) {
                if (exp == (int)exp && (int)exp % 2 != 0) {
                    res = -0.0;
                } else {
                    res = 0;
                }
            } else {
                if (exp == (int)exp && (int)exp % 2 != 0) {
                    res = s21_INFINITY_NEGATIVE;
                } else {
                    res = s21_INFINITY_POSITIVE;
                }
            }
        }
    } else {
        // normal base
        if (exp == s21_INFINITY_NEGATIVE) {
            if (base == -1) {
                res = 1;
            } else if (base > -1 && base < 1) {
                res = s21_INFINITY_POSITIVE;
            } else {
                res = 0;
            }
        } else if (exp == s21_INFINITY_POSITIVE) {
            if (base == -1) {
                res = 1;
            } else if (base > -1 && base < 1) {
                res = 0;
            } else {
                res = s21_INFINITY_POSITIVE;
            }
        } else if (base < 0) {
            if (exp == (int)exp) {
                if ((int)exp % 2 == 0) {
                    res = s21_pow(-base, exp);
                } else {
                    res = -s21_pow(-base, exp);
                }
            } else if (exp != exp) {
                res = exp;
            } else {
                res = s21_NAN_NEGATIVE;
            }
        } else if (base == 0) {
            if (exp != exp) {
                res = exp;
            } else if (exp < 0) {
                res = s21_INFINITY_POSITIVE;
            } else {
                res = 0;
            }
        } else {
            if (exp == (int)exp) {
                res = 1;
                for (int i = 0; i < (int)exp; i++) {
                    res *= base;
                }
                if (exp < 0) {
                    res = 1 / res;
                }
            } else {
                res = s21_exp(exp * s21_log(base));
            }
        }
    }

    return res;
}

long double s21_asin(double x) {
    long double res = 0;
    if (x > 1.0 || x < -1.0) {
        res = s21_NAN;
    } else if (x == -1) {
        res = -s21_M_PI / 2;
    } else if (x == 1) {
        res = s21_M_PI / 2;
    } else {
        long double y = 1 + s21_E;
        int n = 0;
        while (s21_fabs(y) > s21_E) {
            if (n == 0) {
                y = x;
            } else {
                // умножаем и делим по очереди, во избежание переполнения числа
                y = y * (2 * n - 1) / (2 * n + 1) * (2 * n - 1) / (2 * n) * x * x;
            }
            res += y;
            n++;
        }
    }
    return res;
}



long double s21_acos(double x) {
    return (s21_M_PI / 2) - s21_asin(x);
}

long double s21_sqrt(double x) {
    return s21_pow(x, 0.5);
}

long double s21_atan(double x) {
    long double res = 0;
    const long double s21_atan_1 = 0.7853981633974480L;
    if (x == s21_INFINITY_POSITIVE || x == s21_INFINITY_NEGATIVE) {
        if (x == s21_INFINITY_POSITIVE) {
            res = s21_M_PI / 2;
        } else if (x == s21_INFINITY_NEGATIVE) {
            res = -s21_M_PI / 2;
        }
    } else if (x == 1) {
        res = s21_atan_1;
    } else if (x == -1) {
        res = -s21_atan_1;
    } else if (x > 1.0 || x < -1.0) {
        res = s21_NAN;
    } else {
        long double y = 1 + s21_E;
        int n = 1;
        int sign = 1;
        while (s21_fabs(y) > s21_E) {
            y = sign * s21_pow(x, 2 * n - 1) / (2 * n - 1);
            res += y;
            n++;
            sign *= -1;
        }
    }
    return res;
}


long double s21_sin(double x) {
    long double res = 0;
    if (x == s21_NAN_NEGATIVE || x == s21_INFINITY_NEGATIVE || x == s21_INFINITY_POSITIVE) {
        res = s21_NAN_NEGATIVE;
    } else {
        double convert_x = convertTrigonometrik(x);
        long double y = s21_E + 1;
        long int n = 0;
        int firstPass = 1;
        while (s21_fabs(y) > s21_E) {
            if (firstPass == 1) {
                y = convert_x;
            } else {
                y = - y * convert_x / (2 * n) * convert_x / (2 * n + 1);
            }
            firstPass = 0;
            res += y;
            n++;
        }
    }

    return res;
}

long double s21_cos(double x) {
    long double res = 0;
    if (x == s21_NAN_NEGATIVE || x == s21_INFINITY_POSITIVE || x == s21_INFINITY_NEGATIVE) {
        res = s21_NAN_NEGATIVE;
    } else {
        double convert_x = convertTrigonometrik(x);
        double y = s21_E + 1;
        int n = 0;
        int firstPass = 1;
        while (s21_fabs(y) > s21_E) {
            if (firstPass == 1) {
                y = 1;
            } else {
                y = - y * convert_x / (2 * n - 1) * convert_x / (2 * n);
            }
            firstPass = 0;
            res += y;
            n++;
        }
    }

    return res;
}

long double s21_tan(double x) {
    return s21_sin(x) / s21_cos(x);
}


long double s21_log(double x) {
    long double res;
    if (x < 0) {
        res = s21_NAN;
    } else if (x == 0) {
        res = s21_INFINITY_NEGATIVE;
    } else if (x == s21_INFINITY_POSITIVE) {
        res = s21_INFINITY_POSITIVE;
    } else if (x != x) {
        res = x;
    } else {
        double t = (x - 1) / (x + 1);
        res = log_mercator(t) - log_mercator(-t);
    }
    return res;
}


/**
 * @param x любое положительное число от 0 до +inf
 * @param x любое отрицательное число от -inf до 0
 * @return y, такой что sin(x) = sin(y) и cos(x) = cos(y), причём y от [-pi, pi]
 **/
double convertTrigonometrik(double x) {
    double res;
    if (x > 0) {
        res = x - ((int)(x / (2 * s21_M_PI)) * (2 * s21_M_PI));
    } else {
        res = x - ((int)(x / (-2 * s21_M_PI)) * (-2 * s21_M_PI));
    }
    return res;
}


/**
 * Вспомогательная функция, вычисляет по ряду Меркатора ln(1+x), где -1 < x < 1.
 **/ 
long double log_mercator(double x) {
    long double res = 0;
    long double y = s21_E + 1;
    int n = 1;
    long double pow = 1;
    int sign = 1;
    while (s21_fabs(y) > s21_E) {
        pow *= x;
        y = sign * pow / n;
        res += y;
        sign *= -1;
        n++;
    }
    return res;
}
