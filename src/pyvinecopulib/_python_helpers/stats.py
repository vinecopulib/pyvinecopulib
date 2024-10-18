import math


# Cumulative Distribution Function (CDF)
def norm_cdf(x, mean=0, std=1):
  return 0.5 * (1 + math.erf((x - mean) / (std * math.sqrt(2))))


# Percent Point Function (Inverse CDF, PPF)
def norm_ppf(p, mean=0, std=1):
  return mean + std * math.sqrt(2) * math.erfinv(2 * p - 1)


# Probability Density Function (PDF)
def norm_pdf(x, mean=0, std=1):
  return (1 / (std * math.sqrt(2 * math.pi))) * math.exp(
    -0.5 * ((x - mean) / std) ** 2
  )


# Cumulative Distribution Function (CDF)
def expon_cdf(x, scale=1):
  return 1 - math.exp(-x / scale)


# Percent Point Function (Inverse CDF, PPF)
def expon_ppf(p, scale=1):
  return -scale * math.log(1 - p)


# Probability Density Function (PDF)
def expon_pdf(x, scale=1):
  return (1 / scale) * math.exp(-x / scale)