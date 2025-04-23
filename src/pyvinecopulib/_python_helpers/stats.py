import math

import numpy as np


# https://stackoverflow.com/questions/42381244/pure-python-inverse-error-function
# This is a pure Python implementation of the inverse error function from the scipy library.
def polevl(x, coefs, N):
  ans = 0
  power = len(coefs) - 1
  for coef in coefs:
    ans += coef * x**power
    power -= 1
  return ans


def p1evl(x, coefs, N):
  return polevl(x, [1] + coefs, N)


def inv_erf(z):
  if z < -1 or z > 1:
    raise ValueError("`z` must be between -1 and 1 inclusive")

  if z == 0:
    return 0
  if z == 1:
    return math.inf
  if z == -1:
    return -math.inf

  # From scipy special/cephes/ndrti.c
  def ndtri(y):
    # approximation for 0 <= abs(z - 0.5) <= 3/8
    P0 = [
      -5.99633501014107895267e1,
      9.80010754185999661536e1,
      -5.66762857469070293439e1,
      1.39312609387279679503e1,
      -1.23916583867381258016e0,
    ]

    Q0 = [
      1.95448858338141759834e0,
      4.67627912898881538453e0,
      8.63602421390890590575e1,
      -2.25462687854119370527e2,
      2.00260212380060660359e2,
      -8.20372256168333339912e1,
      1.59056225126211695515e1,
      -1.18331621121330003142e0,
    ]

    # Approximation for interval z = sqrt(-2 log y ) between 2 and 8
    # i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
    P1 = [
      4.05544892305962419923e0,
      3.15251094599893866154e1,
      5.71628192246421288162e1,
      4.40805073893200834700e1,
      1.46849561928858024014e1,
      2.18663306850790267539e0,
      -1.40256079171354495875e-1,
      -3.50424626827848203418e-2,
      -8.57456785154685413611e-4,
    ]

    Q1 = [
      1.57799883256466749731e1,
      4.53907635128879210584e1,
      4.13172038254672030440e1,
      1.50425385692907503408e1,
      2.50464946208309415979e0,
      -1.42182922854787788574e-1,
      -3.80806407691578277194e-2,
      -9.33259480895457427372e-4,
    ]

    # Approximation for interval z = sqrt(-2 log y ) between 8 and 64
    # i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
    P2 = [
      3.23774891776946035970e0,
      6.91522889068984211695e0,
      3.93881025292474443415e0,
      1.33303460815807542389e0,
      2.01485389549179081538e-1,
      1.23716634817820021358e-2,
      3.01581553508235416007e-4,
      2.65806974686737550832e-6,
      6.23974539184983293730e-9,
    ]

    Q2 = [
      6.02427039364742014255e0,
      3.67983563856160859403e0,
      1.37702099489081330271e0,
      2.16236993594496635890e-1,
      1.34204006088543189037e-2,
      3.28014464682127739104e-4,
      2.89247864745380683936e-6,
      6.79019408009981274425e-9,
    ]

    s2pi = 2.50662827463100050242
    code = 1

    if y > (1.0 - 0.13533528323661269189):  # 0.135... = exp(-2)
      y = 1.0 - y
      code = 0

    if y > 0.13533528323661269189:
      y = y - 0.5
      y2 = y * y
      x = y + y * (y2 * polevl(y2, P0, 4) / p1evl(y2, Q0, 8))
      x = x * s2pi
      return x

    x = math.sqrt(-2.0 * math.log(y))
    x0 = x - math.log(x) / x

    z = 1.0 / x
    if x < 8.0:  # y > exp(-32) = 1.2664165549e-14
      x1 = z * polevl(z, P1, 8) / p1evl(z, Q1, 8)
    else:
      x1 = z * polevl(z, P2, 8) / p1evl(z, Q2, 8)

    x = x0 - x1
    if code != 0:
      x = -x

    return x

  result = ndtri((z + 1) / 2.0) / math.sqrt(2)

  return result


# Vectorized version of custom inv_erf function
erf_vec = np.vectorize(math.erf)
erfinv_vec = np.vectorize(inv_erf)


# Cumulative Distribution Function (CDF)
def norm_cdf(x, mean=0, std=1):
  return 0.5 * (1 + erf_vec((x - mean) / (std * np.sqrt(2))))


# Percent Point Function (Inverse CDF, PPF)
def norm_ppf(p, mean=0, std=1):
  return mean + std * np.sqrt(2) * erfinv_vec(2 * p - 1)


# Probability Density Function (PDF)
def norm_pdf(x, mean=0, std=1):
  return (1 / (std * np.sqrt(2 * np.pi))) * np.exp(
    -0.5 * ((x - mean) / std) ** 2
  )


# Cumulative Distribution Function (CDF)
def expon_cdf(x, scale=1):
  return 1 - np.exp(-x / scale)


# Percent Point Function (Inverse CDF, PPF)
def expon_ppf(p, scale=1):
  return -scale * np.log(1 - p)


# Probability Density Function (PDF)
def expon_pdf(x, scale=1):
  return (1 / scale) * np.exp(-x / scale)
