#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_
#define s21_M_E 2.7182818284590452354
#define S21_EPS 1e-9
#define s21_E 0.00000000000000000001
#define s21_M_PI 3.14159265358979323846
#define s21_INFINITY_NEGATIVE -1.0/0.0
#define s21_INFINITY_POSITIVE 1.0/0.0
#define s21_NAN -(0.0/0.0)
#define s21_NAN_NEGATIVE 0.0 / 0.0
#define is_fin(x) __builtin_isfinite(x)
#define is_nan(x) __builtin_isnan(x)
#define is_inf(x) __builtin_isinf(x)


int s21_abs(int x);
long double s21_fabs(double x);
long double s21_ceil(double x);
long double s21_floor(double x);
long double s21_fmod(double x, double y);
long double s21_exp(double x);
long double s21_pow(double base, double exp);
long double s21_sin(double x);
long double s21_sqrt(double x);
long double s21_cos(double x);
long double s21_tan(double x);
long double s21_log(double x);
long double s21_asin(double x);
long double s21_acos(double x);
long double s21_atan(double x);
double convertTrigonometrik(double x);
long double log_mercator(double x);

#endif  // SRC_S21_MATH_H_
