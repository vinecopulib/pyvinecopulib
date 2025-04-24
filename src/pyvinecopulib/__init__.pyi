import collections
from typing import Any
from numpy.typing import ArrayLike

class Bicop:
  """
  A class for bivariate copula models.
  
  The model is fully characterized by the family, rotation (one of ``0``, ``90``, ``180``, ``270``), a
  matrix of parameters, and variable types (two strings, one for each variable, either ``"c"`` for
  continuous or ``"d"`` for discrete).
  
  Implemented families (see ``BicopFamily``):
  
  
  ::
  
      | type          | full name             | string identifier     |
      |---------------|-----------------------|-----------------------|
      | -             | Independence          | "indep"               |
      | Elliptical    | Gaussian              | "gaussian"            |
      |               | Student t             | "student"             |
      | Archimedean   | Clayton               | "clayton"             |
      |               | Gumbel                | "gumbel"              |
      |               | Frank                 | "frank"               |
      |               | Joe                   | "joe"                 |
      |               | Clayton-Gumbel (BB1)  | "bb1"                 |
      |               | Joe-Gumbel (BB6)      | "bb6"                 |
      |               | Joe-Clayton (BB7)     | "bb7"                 |
      |               | Joe-Frank (BB8)       | "bb8"                 |
      | Extreme-Value | Tawn                  | "tawn"                |
      | Nonparametric | Transformation kernel | "tll"                 |
  """
  def __init__(self, family: "BicopFamily" = ..., rotation: int = 0, parameters: ArrayLike = ..., var_types: collections.abc.Sequence[str] = ['c', 'c']) -> None:
    """
    
    Default constructor for the ``Bicop`` class.
    
    The default constructor uses ``Bicop.from_family()`` to instantiate an
    independent bivariate copula. It can then be used to select a model from data using ``Bicop.select()``. Or if a ``BicopFamily`` is passed to the constructor, then the method ``Bicop.fit()`` can be used to fit the copula to data.
    
    Alternatives to instantiate bivariate copulas are:
    
    - ``Bicop.from_family()``: Instantiate from a family, rotation, parameters, and variable types.
    - ``Bicop.from_data()``: Instantiate from data, as well as optional controls and variable types.
    - ``Bicop.from_file()``: Instantiate from a file.
    - ``Bicop.from_json()``: Instantiate from a JSON string.
    """
    ...

  def aic(self, u: ArrayLike = ...) -> float:
    """
    
    Evaluates the Akaike information criterion (AIC).
    
    The AIC is defined as
    
    .. math:: \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p,
    
    where :math:`\mathrm{loglik}` is the log-liklihood (see ``Bicop.loglik()``) and :math:`p` is the
    (effective) number of parameters of the model. The AIC is a consistent model selection criterion
    even for nonparametric models.
    
    Parameters
    ----------
    u :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    Returns
    -------
    The AIC evaluated at ``u``.
    """
    ...

  def bic(self, u: ArrayLike = ...) -> float:
    """
    
    Evaluates the Bayesian information criterion (BIC).
    
    The BIC is defined as
    
    .. math:: \mathrm{BIC} = -2\, \mathrm{loglik} + \log(n) p,
    
    where :math:`\mathrm{loglik}` is the log-liklihood (see ``Bicop.loglik()``) and :math:`p` is the
    (effective) number of parameters of the model. The BIC is a consistent model selection criterion for
    parametric models.
    
    Parameters
    ----------
    u :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    Returns
    -------
    The BIC evaluated at ``u``.
    """
    ...

  def cdf(self, u: ArrayLike) -> ArrayLike:
    """
    
    Evaluates the copula distribution.
    
    When at least one variable is discrete, more than two columns are required for ``u``: the first
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1), F_{X_2}(x_2))`. The second
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1^-), F_{X_2}(x_2^-))`. The
    minus indicates a left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
    :math:`F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1)`. For continuous variables the left limit and the cdf
    itself coincide. Respective columns can be omitted in the second block.
    
    Parameters
    ----------
    u :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    Returns
    -------
    A length n vector of copula probabilities evaluated at ``u``.
    """
    ...

  @property
  def family(self) -> Any: ...
  def fit(self, data: ArrayLike, controls: "FitControlsBicop" = ...) -> None:
    """
    
    Fits a bivariate copula (with fixed family) to data.
    
    For parametric models, two different methods are available. ``"mle"`` fits the parameters by
    maximum-likelihood. ``"itau"`` uses inversion of Kendall's :math:`\tau`, but is only available for
    one-parameter families and the Student t copula. For the latter, there is a one-to-one
    transformation for the first parameter, the second is found by profile likelihood optimization (with
    accuracy of at least 0.5). Nonparametric families have specialized methods, no specification is
    required.
    
    When at least one variable is discrete, two types of "observations" are required: the first :math:`n
    \times 2` block contains realizations of :math:`F_{X_1}(X_1), F_{X_2}(X_2)`. Let :math:`k` denote
    the number of discrete variables (either one or two). Then the second :math:`n \times k` block
    contains realizations of :math:`F_{X_k}(X_k^-)`. The minus indicates a left-sided limit of the cdf.
    For continuous variables the left limit and the cdf itself coincide. For, e.g., an integer-valued
    variable, it holds :math:`F_{X_k}(X_k^-) = F_{X_k}(X_k - 1)`.
    
    Incomplete observations (i.e., ones with a NaN value) are discarded.
    
    Parameters
    ----------
    data :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    controls :
        The controls (see ``FitControlsBicop``).
    """
    ...


  @staticmethod
  def from_data(data: ArrayLike, controls: "FitControlsBicop" = ..., var_types: collections.abc.Sequence[str] = ['c', 'c']) -> "Bicop":
    """
    
    Instantiates from data.
    
    Equivalent to creating a default ``Bicop()`` and then selecting the model using ``Bicop.select()``.
    
    Parameters
    ----------
    data :
        See ``Bicop.select()``.
    
    controls :
        See ``Bicop.select()``.
    
    var_types :
        Two strings specifying the types of the variables, e.g., ``("c", "d")`` means first variable
        continuous, second discrete.
    """
    ...


  @staticmethod
  def from_family(family: "BicopFamily" = ..., rotation: int = 0, parameters: ArrayLike = ..., var_types: collections.abc.Sequence[str] = ['c', 'c']) -> "Bicop":
    """
    
    Instantiates a specific bivariate copula model.
    
    Parameters
    ----------
    family :
        The copula family.
    
    rotation :
        The rotation of the copula; one of 0, 90, 180, or 270 (for Independence, Gaussian, Student,
        Frank, and nonparametric families, only 0 is allowed).
    
    parameters :
        The copula parameters.
    
    var_types :
        Two strings specifying the types of the variables, e.g., ``("c", "d")`` means first variable
        continuous, second discrete.
    """
    ...


  @staticmethod
  def from_file(filename: str) -> "Bicop":
    """
    
    Instantiates from a JSON file.
    
    The input file contains four attributes: ``"fam"``, ``"rot"``, ``"par"``, ``"vt"`` respectively a
    string for the family name, an integer for the rotation, and a numeric matrix for the parameters,
    and a list of two strings for the variable types.
    
    Parameters
    ----------
    filename :
        The name of the JSON file to read.
    """
    ...


  @staticmethod
  def from_json(json: str) -> "Bicop":
    """
    
    Instantiates from a JSON-like `str` object.
    
    Parameters
    ----------
    input :
        The JSON-like `str` object to convert from (see ``to_json()`` for the structure of the
        input).
    """
    ...

  def hfunc1(self, u: ArrayLike) -> ArrayLike:
    """
    
    Evaluates the first h-function.
    
    The first h-function is :math:`h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 = u_1)`.
    
    When at least one variable is discrete, more than two columns are required for ``u``: the first
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1), F_{X_2}(x_2))`. The second
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1^-), F_{X_2}(x_2^-))`. The
    minus indicates a left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
    :math:`F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1)`. For continuous variables the left limit and the cdf
    itself coincide. Respective columns can be omitted in the second block.
    
    Parameters
    ----------
    u :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    Returns
    -------
    A length n vector of the first h-function evaluated at ``u``.
    """
    ...

  def hfunc2(self, u: ArrayLike) -> ArrayLike:
    """
    
    Evaluates the second h-function.
    
    The second h-function is :math:`h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 = u_2)`.
    
    When at least one variable is discrete, more than two columns are required for ``u``: the first
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1), F_{X_2}(x_2))`. The second
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1^-), F_{X_2}(x_2^-))`. The
    minus indicates a left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
    :math:`F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1)`. For continuous variables the left limit and the cdf
    itself coincide. Respective columns can be omitted in the second block.
    
    Parameters
    ----------
    u :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    Returns
    -------
    A length n vector of the second h-function evaluated at ``u``.
    """
    ...

  def hinv1(self, u: ArrayLike) -> ArrayLike:
    """
    
    Evaluates the inverse of the first h-function.
    
    The first h-function is :math:`h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 = u_1)`. The inverse is calulated
    w.r.t. the second argument.
    
    When at least one variable is discrete, more than two columns are required for ``u``: the first
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1), F_{X_2}(x_2))`. The second
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1^-), F_{X_2}(x_2^-))`. The
    minus indicates a left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
    :math:`F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1)`. For continuous variables the left limit and the cdf
    itself coincide. Respective columns can be omitted in the second block.
    
    Parameters
    ----------
    u :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    Returns
    -------
    A length n vector of the inverse of the first h-function evaluated at ``u``.
    """
    ...

  def hinv2(self, u: ArrayLike) -> ArrayLike:
    """
    
    Evaluates the inverse of the second h-function.
    
    The second h-function is :math:`h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 = u_2)`. The inverse is
    calculated w.r.t. the first argument.
    
    When at least one variable is discrete, more than two columns are required for ``u``: the first
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1), F_{X_2}(x_2))`. The second
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1^-), F_{X_2}(x_2^-))`. The
    minus indicates a left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
    :math:`F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1)`. For continuous variables the left limit and the cdf
    itself coincide. Respective columns can be omitted in the second block.
    
    Parameters
    ----------
    u :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    Returns
    -------
    A length n vector of the inverse of the second h-function evaluated at ``u``.
    """
    ...

  def loglik(self, u: ArrayLike = ...) -> float:
    """
    
    Evaluates the log-likelihood.
    
    The log-likelihood is defined as
    
    .. math:: \mathrm{loglik} = \sum_{i = 1}^n \log c(U_{1, i}, U_{2, i}),
    
    where :math:`c` is the copula density, see ``Bicop.pdf()``.
    
    When at least one variable is discrete, more than two columns are required for ``u``: the first
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1), F_{X_2}(x_2))`. The second
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1^-), F_{X_2}(x_2^-))`. The
    minus indicates a left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
    :math:`F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1)`. For continuous variables the left limit and the cdf
    itself coincide. Respective columns can be omitted in the second block.
    
    Parameters
    ----------
    u :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    Returns
    -------
    The log-likelihood evaluated at ``u``.
    """
    ...

  def mbic(self, u: ArrayLike = ..., psi0: float = 0.9) -> float:
    """
    
    Evaluates the modified Bayesian information criterion (mBIC).
    
    The mBIC is defined as
    
    .. math:: \mathrm{BIC} = -2\, \mathrm{loglik} + p \log(n) - 2 (I \log(\psi_0) + (1 - I) \log(1 -
    \psi_0),
    
    where :math:`\mathrm{loglik}` is the \log-liklihood (see ``Bicop.loglik()``), :math:`p` is the
    (effective) number of parameters of the model, and :math:`\psi_0` is the prior probability of having
    a non-independence copula and :math:`I` is an indicator for the family being non-independence.
    
    Parameters
    ----------
    u :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    psi0 :
        Prior probability of a non-independence copula.
    
    Returns
    -------
    The mBIC evaluated at ``u``.
    """
    ...

  @property
  def nobs(self) -> Any: ...
  @property
  def npars(self) -> Any: ...
  @property
  def parameters(self) -> Any: ...
  @parameters.setter
  def parameters(self, value: Any) -> None: ...
  @property
  def parameters_lower_bounds(self) -> Any: ...
  def parameters_to_tau(self, parameters: ArrayLike) -> float:
    """
    
    Converts the copula parameters to Kendall's :math:`tau`.
    
    Parameters
    ----------
    parameters :
        The parameters (must be a valid parametrization of the current family).
    """
    ...

  @property
  def parameters_upper_bounds(self) -> Any: ...
  def pdf(self, u: ArrayLike) -> ArrayLike:
    """
    
    Evaluates the copula density.
    
    The copula density is defined as joint density divided by marginal densities, irrespective of
    variable types.
    
    When at least one variable is discrete, more than two columns are required for ``u``: the first
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1), F_{X_2}(x_2))`. The second
    :math:`n \times 2` block contains realizations of :math:`(F_{X_1}(x_1^-), F_{X_2}(x_2^-))`. The
    minus indicates a left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
    :math:`F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1)`. For continuous variables the left limit and the cdf
    itself coincide. Respective columns can be omitted in the second block.
    
    Parameters
    ----------
    u :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    Returns
    -------
    A length n vector of copula densities evaluated at ``u``.
    """
    ...

  def plot(self, type: str = 'surface', margin_type: str = 'unif', xylim: object | None = None, grid_size: object | None = None) -> None:
    """
    
    Generates a plot for the Bicop object.
    
    This method generates a contour or surface plot of the copula density. It can be used to visualize the copula density with different types of margins.
    
    Parameters
    ----------
    plot_type : str (default="contour")
        The type of plot to generate. Either `"contour"` or `"surface"`.
    margin_type : str (default="unif")
        The type of margins to use. Either `"unif"`, `"norm"`, or `"exp"`.
    xylim : tuple (default=None)
        The limits for the x and y axes. Automatically set if None.
    grid_size : int (default=None)
        The number of grid points to use. Automatically set if None.
    
    Returns
    -------
    Nothing, the function generates a plot and shows it using matplotlib.
    
    Usage
    -----
    .. code-block:: python
    
        import pyvinecopulib as pv
        import numpy as np
        cop = pv.Bicop(family=pv.BicopFamily.gaussian, parameters=np.array([[0.5]]))
        cop.plot() # surface plot of copula density
        cop.plot(plot_type="contour", margin_type="norm") # contour plot with normal margins
        cop.plot(plot_type="contour", margin_type="unif") # contour plot of copula density
    """
    ...

  @property
  def rotation(self) -> Any: ...
  @rotation.setter
  def rotation(self, value: Any) -> None: ...
  def select(self, data: ArrayLike, controls: "FitControlsBicop" = ...) -> None:
    """
    
    Selects the best fitting model.
    
    The function calls ``Bicop.fit()`` for all families in ``family_set`` and selecting the best
    fitting model by either BIC or AIC, see ``Bicop.bic()`` and ``Bicop.aic()``.
    
    When at least one variable is discrete, two types of "observations" are required: the first :math:`n
    \times 2` block contains realizations of :math:`F_{X_1}(X_1), F_{X_2}(X_2)`. Let :math:`k` denote
    the number of discrete variables (either one or two). Then the second :math:`n \times k` block
    contains realizations of :math:`F_{X_k}(X_k^-)`. The minus indicates a left-sided limit of the cdf.
    For continuous variables the left limit and the cdf itself coincide. For, e.g., an integer-valued
    variable, it holds :math:`F_{X_k}(X_k^-) = F_{X_k}(X_k - 1)`.
    
    Incomplete observations (i.e., ones with a NaN value) are discarded.
    
    Parameters
    ----------
    data :
        An :math:`n \times (2 + k)` matrix of observations contained in :math:`(0, 1)`, where :math:`k`
        is the number of discrete variables.
    
    controls :
        The controls (see ``FitControlsBicop``).
    """
    ...

  def simulate(self, n: int, qrng: bool = False, seeds: collections.abc.Sequence[int] = []) -> ArrayLike:
    """
    
    Simulates from a bivariate copula.
    
    If ``qrng = TRUE``, generalized Halton sequences are used. For more information on Generalized
    Halton sequences, see Faure, H., Lemieux, C. (2009). Generalized Halton Sequences in 2008: A
    Comparative Study. ACM-TOMACS 19(4), Article 15.
    
    Parameters
    ----------
    n :
        Number of observations.
    
    qrng :
        Set to true for quasi-random numbers.
    
    seeds :
        Seeds of the (quasi-)random number generator; if empty (default), the (quasi-)random number
        generator is seeded randomly.
    
    Returns
    -------
    An :math:`n \times 2` matrix of samples from the copula model.
    """
    ...

  @property
  def tau(self) -> Any: ...
  def tau_to_parameters(self, tau: float) -> ArrayLike:
    """
    
    Converts a Kendall's :math:`\tau` into copula parameters for one-parameter families.
    
    Parameters
    ----------
    tau :
        A value in :math:`(-1, 1)`.
    """
    ...

  def to_file(self, filename: str) -> None:
    """
    
    Write the copula object into a JSON file.
    
    The written file contains four attributes: ``"fam"``, ``"rot"``, ``"par"``, ``"vt"``, ``"nobs"``,
    ``"ll"``, ``"npars"`` respectively a string for the family name, an integer for the rotation, and a
    numeric matrix for the parameters, a list of two strings for the variable types, an integer for the
    number of observations (if fitted), a double for the log-likelihood (if fitted), and a double for
    the number of parameters (can be non-integer in nonparametric models).
    
    Parameters
    ----------
    filename :
        The name of the file to write.
    """
    ...

  def to_json(self) -> str:
    """
    
    Convert the copula into a nlohmann::json object.
    
    The JSON-like `str` is contains of three values named ``"fam"``, ``"rot"``, ``"par"``, ``"vt"``,
    respectively a string for the family name, an integer for the rotation, a numeric matrix for the
    parameters and a list of two strings for the variables types.
    
    Returns
    -------
    The JSON-like `str` object containing the copula.
    """
    ...

  @property
  def var_types(self) -> Any: ...
  @var_types.setter
  def var_types(self, value: Any) -> None: ...

class BicopFamily:
  """
  A bivariate copula family identifier.
  
  Contains the following families:
  
    - ``indep``: Independent copula,
    - ``gaussian``: Gaussian copula,
    - ``student``: Student t copula,
    - ``clayton``: Clayton copula,
    - ``gumbel``: Gumbel copula,
    - ``frank``: Frank copula,
    - ``joe``: Joe copula,
    - ``bb1``: BB1 copula,
    - ``bb6``: BB6 copula,
    - ``bb7``: BB7 copula,
    - ``bb8``: BB8 copula,
    - ``tawn``: Tawn copula,
    - ``tll``: Transformation local-likelihood (nonparametric) copula.
  
  The following convenient sets of families are provided:
  
    - ``all`` contains all the families,
    - ``parametric`` contains the parametric families (all except ``tll``),
    - ``nonparametric`` contains the nonparametric families
      (``indep`` and ``tll``),
    - ``one_par`` contains the parametric families with a single parameter
      (``gaussian``, ``clayton``, ``gumbel``, ``frank``, and ``joe``),
    - ``two_par`` contains the parametric families with two parameters
      (``student``, ``bb1``, ``bb6``, ``bb7``, and ``bb8``),
    - ``three_par`` contains the parametric families with three parameters
      (``tawn``),
    - ``elliptical`` contains the elliptical families (``gaussian`` and
      ``student``),
    - ``archimedean`` contains the archimedean families (``clayton``,
      ``gumbel``, ``frank``, ``joe``, ``bb1``, ``bb6``, ``bb7``, and ``bb8``),
    - ``extreme_value`` contains the extreme value families (``tawn`` and
      ``gumbel``),
    - ``bb`` contains the BB families (``bb1``, ``bb6``, ``bb7``, and 
      ``bb8``),
    - ``itau`` families for which estimation by Kendall's tau inversion is
      available (``indep``, ``gaussian``, ``student``, ``clayton``,
      ``gumbel``, ``frank``, ``joe``),
    - ``lt`` contains the families that are lower-tail dependent (``clayton``,
      ``bb1``, ``bb7``, ``tawn``),
    - ``ut`` contains the families that are upper-tail dependent (``gumbel``,
      ``joe``, ``bb1``, ``bb6``, ``bb7``, ``bb8``, ``tawn``),
    - ``rotationless`` contains families that don't have a rotation 
      because they already cover positive and negative dependence (``indep``, 
      ``gaussian``, ``student``, ``frank``, ``tll``).
  """

class CVineStructure:
  """
  A class for C-vine structures.
  
  C-vines are a special class of R-vines where each tree is a star. A C-vine structure is determined
  entirely by the order of variables. For example, if the order is ``{1, 2, 3, 4}``, the first tree in
  the vine connects variable 4 with all others, the second tree connects variable 3 with all others,
  etc.
  
  Note that ``CVineStructure`` objects inherit the methods and attributes of ``RVineStructure``
  objects.
  """
  def __init__(self, order: collections.abc.Sequence[int], trunc_lvl: int = 18446744073709551615) -> None:
    """
    
    Parameters
    ----------
    order :
        The order of variables in the C-vine (diagonal entries in the R-vine array); must be a
        permutation of 1, ..., d.
    
    trunc_lvl :
        The truncation level.
    """
    ...

  @property
  def dim(self) -> Any: ...

  @staticmethod
  def from_dimension(d: int = 1, trunc_lvl: int = 18446744073709551615) -> "RVineStructure":
    """
    
    Instantiates as a D-vine for a given dimension.
    
    Parameters
    ----------
    d :
        The dimension.
    
    trunc_lvl :
        The truncation level. By default, it is dim - 1.
    """
    ...


  @staticmethod
  def from_file(filename: str, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates an RVineStructure from a JSON file.
    
    The file needs to contain two values: ``"array"`` for the structure triangular array and ``"order"``
    for the order vector.
    
    Parameters
    ----------
    filename :
        The name of the JSON file to read.
    
    check :
        Whether to check if the input represents a valid R-vine matrix.
    """
    ...


  @staticmethod
  def from_json(json: str, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates from a JSON-like `str` object.
    
    Parameters
    ----------
    input :
        The JSON-like `str` object to convert from (see ``to_json()`` for the structure of the
        input).
    
    check :
        Whether to check if the input represents a valid R-vine structure.
    """
    ...


  @staticmethod
  def from_matrix(mat: ArrayLike, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates an RVineStructure object from a matrix representing an R-vine array.
    
    The matrix must contain zeros in the lower right triangle and the upper left triangle must be a
    valid R-vine array. Truncated vines can be encoded by putting zeros above the digonal in all rows
    below the truncation level. Example of a 1-truncated matrix:
    
    
    ::
    
        4 4 4 4
        0 0 3 0
        0 2 0 0
        1 0 0 0
    
    Parameters
    ----------
    mat :
        A matrix representing a valid R-vine array.
    
    check :
        Whether ``mat`` shall be checked for validity.
    """
    ...


  @staticmethod
  def from_order(order: collections.abc.Sequence[int], trunc_lvl: int = 18446744073709551615, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates as a D-vine with a given ordering of the variables.
    
    Parameters
    ----------
    order :
        The order of variables in the D-vine (diagonal entries in the R-vine array); must be a
        permutation of 1, ..., d.
    
    trunc_lvl :
        The truncation level. By default, it is d - 1.
    
    check :
        Whether `order shall be checked for validity.
    """
    ...

  @property
  def matrix(self) -> Any: ...
  def min_array(self, tree: int, edge: int) -> int:
    """
    
    Access elements of the minimum array.
    
    Parameters
    ----------
    tree :
        Tree index.
    
    edge :
        Edge index.
    """
    ...

  def needed_hfunc1(self, tree: int, edge: int) -> bool:
    """
    
    Access elements of the needed_hfunc1 array.
    
    Parameters
    ----------
    tree :
        Tree index.
    
    edge :
        Edge index.
    """
    ...

  def needed_hfunc2(self, tree: int, edge: int) -> bool:
    """
    
    Access elements of the needed_hfunc2 array.
    """
    ...

  @property
  def order(self) -> Any: ...

  @staticmethod
  def simulate(d: int, natural_order: bool = False, seeds: collections.abc.Sequence[int] = []) -> "RVineStructure":
    """
    
    Randomly sample a regular vine structure.
    
    Simulates from a uniform distribution over all R-vine structures on d variables
    
    Implementation of Algorithm 13 in Harry Joe's 2014 book (p. 288), but there's a typo: the end of
    line 6 in the book should be 'column j' instead of 'column k'.
    
    Parameters
    ----------
    d :
        The dimension.
    
    natural_order :
        Should the sampled structure be in natural order?
    
    seeds :
        Seeds of the random number generator; if empty (default), the random number generator is seeded
        randomly.
    """
    ...

  def struct_array(self, tree: int, edge: int, natural_order: bool = False) -> int:
    """
    
    Accesses elements of the structure array.
    
    Parameters
    ----------
    tree :
        Tree index.
    
    edge :
        Edge index.
    
    natural_order :
        Whether indices correspond to natural order.
    """
    ...

  def to_file(self, filename: str) -> None:
    """
    
    Write the structure into a JSON file.
    
    The written file contains two values: ``"array"`` for the structure triangular array and ``"order"``
    for the order vector.
    
    Parameters
    ----------
    filename :
        The name of the file to write.
    """
    ...

  def to_json(self) -> str:
    """
    
    Converts the structure into a JSON-like `str` object.
    
    The JSON-like `str` object contains two nodes: ``"array"`` for the structure triangular array and
    ``"order"`` for the order vector.
    
    Returns
    -------
    The JSON-like `str` object containing the structure.
    """
    ...

  @property
  def trunc_lvl(self) -> Any: ...
  def truncate(self, trunc_lvl: int) -> None:
    """
    
    Truncates the R-vine structure.
    
    While a structure of dimension ``d`` contains at most ``d-1`` nested levels, this function extracts
    a sub-structure based on a given truncation level.
    
    If the structure is already truncated at a level less than ``trunc_lvl``, the function does nothing.
    
    Parameters
    ----------
    trunc_lvl :
        The truncation level.
    """
    ...


class DVineStructure:
  """
  A class for D-vine structures.
  
  D-vines are a special class of R-vines where each tree is a path. A D-vine structure is determined
  entirely by the order of variables. For example, if the order is ``(1, 2, 3, 4)``, the first tree in
  the vine is 1-2-3-4 and all further trees are unique due to the proximity condition.
  
  Note that ``DVineStructure`` objects inherit the methods and attributes of ``RVineStructure``
  objects.
  """
  def __init__(self, order: collections.abc.Sequence[int], trunc_lvl: int = 18446744073709551615) -> None:
    """
    
    Parameters
    ----------
    order :
        The order of variables in the D-vine (diagonal entries in the R-vine array); must be a
        permutation of 1, ..., d.
    
    trunc_lvl :
        The truncation level.
    """
    ...

  @property
  def dim(self) -> Any: ...

  @staticmethod
  def from_dimension(d: int = 1, trunc_lvl: int = 18446744073709551615) -> "RVineStructure":
    """
    
    Instantiates as a D-vine for a given dimension.
    
    Parameters
    ----------
    d :
        The dimension.
    
    trunc_lvl :
        The truncation level. By default, it is dim - 1.
    """
    ...


  @staticmethod
  def from_file(filename: str, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates an RVineStructure from a JSON file.
    
    The file needs to contain two values: ``"array"`` for the structure triangular array and ``"order"``
    for the order vector.
    
    Parameters
    ----------
    filename :
        The name of the JSON file to read.
    
    check :
        Whether to check if the input represents a valid R-vine matrix.
    """
    ...


  @staticmethod
  def from_json(json: str, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates from a JSON-like `str` object.
    
    Parameters
    ----------
    input :
        The JSON-like `str` object to convert from (see ``to_json()`` for the structure of the
        input).
    
    check :
        Whether to check if the input represents a valid R-vine structure.
    """
    ...


  @staticmethod
  def from_matrix(mat: ArrayLike, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates an RVineStructure object from a matrix representing an R-vine array.
    
    The matrix must contain zeros in the lower right triangle and the upper left triangle must be a
    valid R-vine array. Truncated vines can be encoded by putting zeros above the digonal in all rows
    below the truncation level. Example of a 1-truncated matrix:
    
    
    ::
    
        4 4 4 4
        0 0 3 0
        0 2 0 0
        1 0 0 0
    
    Parameters
    ----------
    mat :
        A matrix representing a valid R-vine array.
    
    check :
        Whether ``mat`` shall be checked for validity.
    """
    ...


  @staticmethod
  def from_order(order: collections.abc.Sequence[int], trunc_lvl: int = 18446744073709551615, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates as a D-vine with a given ordering of the variables.
    
    Parameters
    ----------
    order :
        The order of variables in the D-vine (diagonal entries in the R-vine array); must be a
        permutation of 1, ..., d.
    
    trunc_lvl :
        The truncation level. By default, it is d - 1.
    
    check :
        Whether `order shall be checked for validity.
    """
    ...

  @property
  def matrix(self) -> Any: ...
  def min_array(self, tree: int, edge: int) -> int:
    """
    
    Access elements of the minimum array.
    
    Parameters
    ----------
    tree :
        Tree index.
    
    edge :
        Edge index.
    """
    ...

  def needed_hfunc1(self, tree: int, edge: int) -> bool:
    """
    
    Access elements of the needed_hfunc1 array.
    
    Parameters
    ----------
    tree :
        Tree index.
    
    edge :
        Edge index.
    """
    ...

  def needed_hfunc2(self, tree: int, edge: int) -> bool:
    """
    
    Access elements of the needed_hfunc2 array.
    """
    ...

  @property
  def order(self) -> Any: ...

  @staticmethod
  def simulate(d: int, natural_order: bool = False, seeds: collections.abc.Sequence[int] = []) -> "RVineStructure":
    """
    
    Randomly sample a regular vine structure.
    
    Simulates from a uniform distribution over all R-vine structures on d variables
    
    Implementation of Algorithm 13 in Harry Joe's 2014 book (p. 288), but there's a typo: the end of
    line 6 in the book should be 'column j' instead of 'column k'.
    
    Parameters
    ----------
    d :
        The dimension.
    
    natural_order :
        Should the sampled structure be in natural order?
    
    seeds :
        Seeds of the random number generator; if empty (default), the random number generator is seeded
        randomly.
    """
    ...

  def struct_array(self, tree: int, edge: int, natural_order: bool = False) -> int:
    """
    
    Accesses elements of the structure array.
    
    Parameters
    ----------
    tree :
        Tree index.
    
    edge :
        Edge index.
    
    natural_order :
        Whether indices correspond to natural order.
    """
    ...

  def to_file(self, filename: str) -> None:
    """
    
    Write the structure into a JSON file.
    
    The written file contains two values: ``"array"`` for the structure triangular array and ``"order"``
    for the order vector.
    
    Parameters
    ----------
    filename :
        The name of the file to write.
    """
    ...

  def to_json(self) -> str:
    """
    
    Converts the structure into a JSON-like `str` object.
    
    The JSON-like `str` object contains two nodes: ``"array"`` for the structure triangular array and
    ``"order"`` for the order vector.
    
    Returns
    -------
    The JSON-like `str` object containing the structure.
    """
    ...

  @property
  def trunc_lvl(self) -> Any: ...
  def truncate(self, trunc_lvl: int) -> None:
    """
    
    Truncates the R-vine structure.
    
    While a structure of dimension ``d`` contains at most ``d-1`` nested levels, this function extracts
    a sub-structure based on a given truncation level.
    
    If the structure is already truncated at a level less than ``trunc_lvl``, the function does nothing.
    
    Parameters
    ----------
    trunc_lvl :
        The truncation level.
    """
    ...


class FitControlsBicop:
  """
  A class for controlling fits of bivariate copula models.
  """
  def __init__(self, family_set: collections.abc.Sequence["BicopFamily"] = ..., parametric_method: str = 'mle', nonparametric_method: str = 'constant', nonparametric_mult: float = 1.0, selection_criterion: str = 'bic', weights: ArrayLike = ..., psi0: float = 0.9, preselect_families: bool = True, allow_rotations: bool = True, num_threads: int = 1) -> None:
    """
    
    Instantiates the controls for fitting bivariate copula models.
    
    Parameters
    ----------
    family_set :
        The set of copula families to consider (if empty, then all families are included).
    
    parametric_method :
        The fit method for parametric families; possible choices: ``"mle"``, ``"itau"``.
    
    nonparametric_method :
        The fit method for the local-likelihood nonparametric family (TLLs); possible choices:
        ``"constant"``, ``"linear"``, ``"quadratic"``.
    
    nonparametric_mult :
        A factor with which the smoothing parameters are multiplied.
    
    selection_criterion :
        The selection criterion (``"loglik"``, ``"aic"`` or ``"bic"``) for the pair copula families.
    
    weights :
        A vector of weights for the observations.
    
    psi0 :
        Only for ``selection_criterion = "mbic"``, the prior probability of non-independence.
    
    preselect_families :
        Whether to exclude families before fitting based on symmetry properties of the data.
    
    allow_rotations :
        Allow rotations for the families when doing model selection (default: true).
    
    num_threads :
        Number of concurrent threads to use while fitting copulas for different families; never uses
        more than the number of concurrent threads supported by the implementation.
    """
    ...

  @property
  def allow_rotations(self) -> Any: ...
  @allow_rotations.setter
  def allow_rotations(self, value: Any) -> None: ...
  @property
  def family_set(self) -> Any: ...
  @family_set.setter
  def family_set(self, value: Any) -> None: ...
  @property
  def nonparametric_method(self) -> Any: ...
  @nonparametric_method.setter
  def nonparametric_method(self, value: Any) -> None: ...
  @property
  def nonparametric_mult(self) -> Any: ...
  @nonparametric_mult.setter
  def nonparametric_mult(self, value: Any) -> None: ...
  @property
  def num_threads(self) -> Any: ...
  @num_threads.setter
  def num_threads(self, value: Any) -> None: ...
  @property
  def parametric_method(self) -> Any: ...
  @parametric_method.setter
  def parametric_method(self, value: Any) -> None: ...
  @property
  def preselect_families(self) -> Any: ...
  @preselect_families.setter
  def preselect_families(self, value: Any) -> None: ...
  @property
  def psi0(self) -> Any: ...
  @psi0.setter
  def psi0(self, value: Any) -> None: ...
  @property
  def selection_criterion(self) -> Any: ...
  @selection_criterion.setter
  def selection_criterion(self, value: Any) -> None: ...
  @property
  def weights(self) -> Any: ...
  @weights.setter
  def weights(self, value: Any) -> None: ...

class FitControlsVinecop:
  """
  A class for controlling fits of vine copula models.
  """
  def __init__(self, family_set: collections.abc.Sequence["BicopFamily"] = ..., parametric_method: str = 'mle', nonparametric_method: str = 'constant', nonparametric_mult: float = 1.0, trunc_lvl: int = 18446744073709551615, tree_criterion: str = 'tau', threshold: float = 0.0, selection_criterion: str = 'bic', weights: ArrayLike = ..., psi0: float = 0.9, preselect_families: bool = True, select_trunc_lvl: bool = False, select_threshold: bool = False, select_families: bool = True, show_trace: bool = False, num_threads: int = 1, tree_algorithm: str = 'mst_prim', allow_rotations: bool = True, seeds: collections.abc.Sequence[int] = []) -> None:
    """
    
    Instantiates custom controls for fitting vine copula models.
    
    Parameters
    ----------
    family_set :
        The set of copula families to consider (if empty, then all families are included).
    
    parametric_method :
        The fit method for parametric families; possible choices: ``"mle"``, ``"itau"``.
    
    nonparametric_method :
        The fit method for the local-likelihood nonparametric family (TLLs); possible choices:
        ``"constant"``, ``"linear"``, ``"quadratic"``.
    
    nonparametric_mult :
        A factor with which the smoothing parameters are multiplied.
    
    trunc_lvl :
        Truncation level for truncated vines.
    
    tree_criterion :
        The criterion for selecting the spanning tree (``"tau"``, ``"hoeffd"``, ``"rho"``, and
        ``"mcor"`` implemented so far) during the tree-wise structure selection.
    
    threshold :
        For thresholded vines (0 = no threshold).
    
    selection_criterion :
        The selection criterion (``"loglik"``, ``"aic"`` or ``"bic"``) for the pair copula families.
    
    weights :
        A vector of weights for the observations.
    
    psi0 :
        Only for ``selection_criterion = "mbic"``, prior probability of non-independence.
    
    preselect_families :
        Whether to exclude families before fitting based on symmetry properties of the data.
    
    select_trunc_lvl :
        Whether the truncation shall be selected automatically.
    
    select_threshold :
        Whether the threshold parameter shall be selected automatically.
    
    select_families :
        Whether the families shall be selected automatically, or should the method simply update the
        parameters for the pair copulas already present in the model.
    
    show_trace :
        Whether to show a trace of the building progress.
    
    num_threads :
        Number of concurrent threads to use while fitting pair copulas within a tree; never uses more
        than the number of concurrent threads supported by the implementation.
    
    tree_algorithm :
        The algorithm for building the spanning tree (``"mst_prim"``, ``"mst_kruskal"``,
        ``"random_weighted"``, or ``"random_unweighted"``) during the tree-wise structure selection.
        ``"mst_prim"`` and ``"mst_kruskal"`` use Prim's and Kruskal's algorithms respectively to select
        the maximum spanning tree, maximizing the sum of the edge weights (i.e., ``tree_criterion``).
        ``"random_weighted"`` and ``"random_unweighted"`` use Wilson's algorithm to generate a random
        spanning tree, either with probability proportional to the product of the edge weights
        (weighted) or uniformly (unweighted).
    
    allow_rotations :
        Allow rotations for the families when doing model selection (default: true).
    
    seeds :
        A vector of random seeds for the random number generator for parts of the algorithm that are
        randomized (e.g., random tree selection).
    """
    ...

  @property
  def allow_rotations(self) -> Any: ...
  @allow_rotations.setter
  def allow_rotations(self, value: Any) -> None: ...
  @property
  def family_set(self) -> Any: ...
  @family_set.setter
  def family_set(self, value: Any) -> None: ...
  @property
  def nonparametric_method(self) -> Any: ...
  @nonparametric_method.setter
  def nonparametric_method(self, value: Any) -> None: ...
  @property
  def nonparametric_mult(self) -> Any: ...
  @nonparametric_mult.setter
  def nonparametric_mult(self, value: Any) -> None: ...
  @property
  def num_threads(self) -> Any: ...
  @num_threads.setter
  def num_threads(self, value: Any) -> None: ...
  @property
  def parametric_method(self) -> Any: ...
  @parametric_method.setter
  def parametric_method(self, value: Any) -> None: ...
  @property
  def preselect_families(self) -> Any: ...
  @preselect_families.setter
  def preselect_families(self, value: Any) -> None: ...
  @property
  def psi0(self) -> Any: ...
  @psi0.setter
  def psi0(self, value: Any) -> None: ...
  @property
  def seeds(self) -> Any: ...
  @seeds.setter
  def seeds(self, value: Any) -> None: ...
  @property
  def select_families(self) -> Any: ...
  @select_families.setter
  def select_families(self, value: Any) -> None: ...
  @property
  def select_threshold(self) -> Any: ...
  @select_threshold.setter
  def select_threshold(self, value: Any) -> None: ...
  @property
  def select_trunc_lvl(self) -> Any: ...
  @select_trunc_lvl.setter
  def select_trunc_lvl(self, value: Any) -> None: ...
  @property
  def selection_criterion(self) -> Any: ...
  @selection_criterion.setter
  def selection_criterion(self, value: Any) -> None: ...
  @property
  def show_trace(self) -> Any: ...
  @show_trace.setter
  def show_trace(self, value: Any) -> None: ...
  @property
  def threshold(self) -> Any: ...
  @threshold.setter
  def threshold(self, value: Any) -> None: ...
  @property
  def tree_algorithm(self) -> Any: ...
  @tree_algorithm.setter
  def tree_algorithm(self, value: Any) -> None: ...
  @property
  def tree_criterion(self) -> Any: ...
  @tree_criterion.setter
  def tree_criterion(self, value: Any) -> None: ...
  @property
  def trunc_lvl(self) -> Any: ...
  @trunc_lvl.setter
  def trunc_lvl(self, value: Any) -> None: ...
  @property
  def weights(self) -> Any: ...
  @weights.setter
  def weights(self, value: Any) -> None: ...

class RVineStructure:
  """
  A class for R-vine structures.
  
  RVineStructure objects encode the tree structure of the vine, i.e. the conditioned/conditioning
  variables of each edge. It is represented by a triangular array. An exemplary array is
  
  
  ::
  
      4 4 4 4
      3 3 3
      2 2
      1
  
  which encodes the following pair-copulas:
  
  
  ::
  
      | tree | edge | pair-copulas |
      |------|------|--------------|
      | 0    | 0    | (1, 4)       |
      |      | 1    | (2, 4)       |
      |      | 2    | (3, 4)       |
      | 1    | 0    | (1, 3; 4)    |
      |      | 1    | (2, 3; 4)    |
      | 2    | 0    | (1, 2; 3, 4) |
  
  Denoting by ``M[i, j]`` the array entry in row ``i`` and column ``j``, the pair-copula index for
  edge ``e`` in tree ``t`` of a ``d`` dimensional vine is ``(M[d - 1 - e, e], M[t, e]; M[t - 1, e],
  ..., M[0, e])``. Less formally,
  
    1. Start with the counter-diagonal element of column ``e``
       (first conditioned variable).
    2. Jump up to the element in row ``t`` (second conditioned variable).
    3. Gather all entries further up in column ``e`` (conditioning set).
  
  Internally, the diagonal is stored separately from the off-diagonal elements, which are stored as a
  triangular array. For instance, the off-diagonal elements off the structure above are stored as
  
  
  ::
  
      4 4 4
      3 3
      2
  
  for the structure above. The reason is that it allows for parsimonious representations of truncated
  models. For instance, the 2-truncated model is represented by the same diagonal and the following
  truncated triangular array:
  
  
  ::
  
      4 4 4
      3 3
  
  A valid R-vine array must satisfy several conditions which are checked when ``RVineStructure()`` is
  called:
  
    1. It only contains numbers between 1 and d.
    2. The diagonal must contain the numbers 1, ..., d.
    3. The diagonal entry of a column must not be contained in any
       column further to the right.
    4. The entries of a column must be contained in all columns to the left.
    5. The proximity condition must hold: For all t = 1, ..., d - 2 and
       e = 0, ..., d - t - 1 there must exist an index j > d, such that
       ``(M[t, e], {M[0, e], ..., M[t-1, e]})`` equals either
       ``(M[d-j-1, j], {M[0, j], ..., M[t-1, j]})`` or
       ``(M[t-1, j], {M[d-j-1, j], M[0, j], ..., M[t-2, j]})``.
  
  An R-vine array is said to be in natural order when the anti-diagonal entries are :math:`1, \dots,
  d` (from left to right). The exemplary arrray above is in natural order. Any R-vine array can be
  characterized by the diagonal entries (called order) and the entries below the diagonal of the
  corresponding R-vine array in natural order. Since most algorithms work with the structure in
  natural order, this is how RVineStructure stores the structure internally.
  """
  def __init__(self, d: int = 1, trunc_lvl: int = 18446744073709551615) -> None:
    """
    
    Default constructor for the ``RVineStructure`` class.
    
    The default constructor uses ``RVineStructure.from_dimension()`` to instantiate
    a default structure of a given dimension and truncation level.
    Alternatives to instantiate structures are:
    
    - ``RVineStructure.from_order()``: Instantiate from an order vector.
    - ``RVineStructure.from_matrix()``: Instantiate from a matrix.
    - ``RVineStructure.from_file()``: Instantiate from a file.
    - ``RVineStructure.from_json()``: Instantiate from a JSON string.
    """
    ...

  @property
  def dim(self) -> Any: ...

  @staticmethod
  def from_dimension(d: int = 1, trunc_lvl: int = 18446744073709551615) -> "RVineStructure":
    """
    
    Instantiates as a D-vine for a given dimension.
    
    Parameters
    ----------
    d :
        The dimension.
    
    trunc_lvl :
        The truncation level. By default, it is dim - 1.
    """
    ...


  @staticmethod
  def from_file(filename: str, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates an RVineStructure from a JSON file.
    
    The file needs to contain two values: ``"array"`` for the structure triangular array and ``"order"``
    for the order vector.
    
    Parameters
    ----------
    filename :
        The name of the JSON file to read.
    
    check :
        Whether to check if the input represents a valid R-vine matrix.
    """
    ...


  @staticmethod
  def from_json(json: str, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates from a JSON-like `str` object.
    
    Parameters
    ----------
    input :
        The JSON-like `str` object to convert from (see ``to_json()`` for the structure of the
        input).
    
    check :
        Whether to check if the input represents a valid R-vine structure.
    """
    ...


  @staticmethod
  def from_matrix(mat: ArrayLike, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates an RVineStructure object from a matrix representing an R-vine array.
    
    The matrix must contain zeros in the lower right triangle and the upper left triangle must be a
    valid R-vine array. Truncated vines can be encoded by putting zeros above the digonal in all rows
    below the truncation level. Example of a 1-truncated matrix:
    
    
    ::
    
        4 4 4 4
        0 0 3 0
        0 2 0 0
        1 0 0 0
    
    Parameters
    ----------
    mat :
        A matrix representing a valid R-vine array.
    
    check :
        Whether ``mat`` shall be checked for validity.
    """
    ...


  @staticmethod
  def from_order(order: collections.abc.Sequence[int], trunc_lvl: int = 18446744073709551615, check: bool = True) -> "RVineStructure":
    """
    
    Instantiates as a D-vine with a given ordering of the variables.
    
    Parameters
    ----------
    order :
        The order of variables in the D-vine (diagonal entries in the R-vine array); must be a
        permutation of 1, ..., d.
    
    trunc_lvl :
        The truncation level. By default, it is d - 1.
    
    check :
        Whether `order shall be checked for validity.
    """
    ...

  @property
  def matrix(self) -> Any: ...
  def min_array(self, tree: int, edge: int) -> int:
    """
    
    Access elements of the minimum array.
    
    Parameters
    ----------
    tree :
        Tree index.
    
    edge :
        Edge index.
    """
    ...

  def needed_hfunc1(self, tree: int, edge: int) -> bool:
    """
    
    Access elements of the needed_hfunc1 array.
    
    Parameters
    ----------
    tree :
        Tree index.
    
    edge :
        Edge index.
    """
    ...

  def needed_hfunc2(self, tree: int, edge: int) -> bool:
    """
    
    Access elements of the needed_hfunc2 array.
    """
    ...

  @property
  def order(self) -> Any: ...

  @staticmethod
  def simulate(d: int, natural_order: bool = False, seeds: collections.abc.Sequence[int] = []) -> "RVineStructure":
    """
    
    Randomly sample a regular vine structure.
    
    Simulates from a uniform distribution over all R-vine structures on d variables
    
    Implementation of Algorithm 13 in Harry Joe's 2014 book (p. 288), but there's a typo: the end of
    line 6 in the book should be 'column j' instead of 'column k'.
    
    Parameters
    ----------
    d :
        The dimension.
    
    natural_order :
        Should the sampled structure be in natural order?
    
    seeds :
        Seeds of the random number generator; if empty (default), the random number generator is seeded
        randomly.
    """
    ...

  def struct_array(self, tree: int, edge: int, natural_order: bool = False) -> int:
    """
    
    Accesses elements of the structure array.
    
    Parameters
    ----------
    tree :
        Tree index.
    
    edge :
        Edge index.
    
    natural_order :
        Whether indices correspond to natural order.
    """
    ...

  def to_file(self, filename: str) -> None:
    """
    
    Write the structure into a JSON file.
    
    The written file contains two values: ``"array"`` for the structure triangular array and ``"order"``
    for the order vector.
    
    Parameters
    ----------
    filename :
        The name of the file to write.
    """
    ...

  def to_json(self) -> str:
    """
    
    Converts the structure into a JSON-like `str` object.
    
    The JSON-like `str` object contains two nodes: ``"array"`` for the structure triangular array and
    ``"order"`` for the order vector.
    
    Returns
    -------
    The JSON-like `str` object containing the structure.
    """
    ...

  @property
  def trunc_lvl(self) -> Any: ...
  def truncate(self, trunc_lvl: int) -> None:
    """
    
    Truncates the R-vine structure.
    
    While a structure of dimension ``d`` contains at most ``d-1`` nested levels, this function extracts
    a sub-structure based on a given truncation level.
    
    If the structure is already truncated at a level less than ``trunc_lvl``, the function does nothing.
    
    Parameters
    ----------
    trunc_lvl :
        The truncation level.
    """
    ...


class Vinecop:
  """
  A class for vine copula models.
  
  A vine copula model is characterized by its structure (see ``RVineStructure`` objects) and the
  pair-copulas (see ``Bicop`` objects).
  """
  def __init__(self, d: int) -> None:
    """
    
    Default constructor for the ``Vinecop`` class.
    
    The default constructor uses ``Vinecop.from_dimension()`` to instantiate an
    empty vine copula of a given dimension. It can then be used to select a model from data using ``Vinecop.select()``. Alternatives to instantiate vine copulas
    are:
    
    - ``Vinecop.from_data()``: Instantiate from data, as well as an optional ``FitControlsVinecop``, an ``RVineStructure`` or matrix, and variable types.
    - ``Vinecop.from_structure()``: Instantiate from an ``RVineStructure`` or matrix, as well as optional pair-copulas and variable types.
    - ``Vinecop.from_file()``: Instantiate from a file.
    - ``Vinecop.from_json()``: Instantiate from a JSON string.
    """
    ...

  def aic(self, u: ArrayLike = ..., num_threads: int = 1) -> float:
    """
    
    Evaluates the Akaike information criterion (AIC).
    
    The AIC is defined as
    
    .. math:: \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p,
    
    where :math:`\mathrm{loglik}` is the log-liklihood (see ``Vinecop.loglik()``) and :math:`p` is the
    (effective) number of parameters of the model. The AIC is a consistent model selection criterion
    even for nonparametric models.
    
    Parameters
    ----------
    u :
        An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of evaluation points, where :math:`k`
        is the number of discrete variables (see ``select()`` or ``Vinecop.pdf()``).
    
    num_threads :
        The number of threads to use for computations; if greater than 1, the function will be applied
        concurrently to ``num_threads`` batches of ``u``.
    
    Returns
    -------
    The AIC as a double.
    """
    ...

  def bic(self, u: ArrayLike = ..., num_threads: int = 1) -> float:
    """
    
    Evaluates the Bayesian information criterion (BIC).
    
    The BIC is defined as
    
    .. math:: \mathrm{BIC} = -2\, \mathrm{loglik} + \log(n) p,
    
    where :math:`\mathrm{loglik}` is the log-liklihood (see ``Vinecop.loglik()``) and :math:`p` is the
    (effective) number of parameters of the model. The BIC is a consistent model selection criterion for
    nonparametric models.
    
    Parameters
    ----------
    u :
        An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of evaluation points, where :math:`k`
        is the number of discrete variables (see ``Vinecop.select()`` or ``Vinecop.pdf()``).
    
    num_threads :
        The number of threads to use for computations; if greater than 1, the function will be applied
        concurrently to ``num_threads`` batches of ``u``.
    
    Returns
    -------
    The BIC as a double.
    """
    ...

  def cdf(self, u: ArrayLike, N: int = 10000, num_threads: int = 1, seeds: collections.abc.Sequence[int] = []) -> ArrayLike:
    """
    
    Evaluates the copula distribution.
    
    Because no closed-form expression is available, the distribution is estimated numerically using
    Monte Carlo integration. The function uses quasi-random numbers from the vine model to do so.
    
    When at least one variable is discrete, two types of "observations" are required in ``u``: the first
    :math:`n \; x \; d` block contains realizations of :math:`F_{X_j}(X_j)`. The second :math:`n \; x \;
    d` block contains realizations of :math:`F_{X_j}(X_j^-)`. The minus indicates a left-sided limit of
    the cdf. For, e.g., an integer-valued variable, it holds :math:`F_{X_j}(X_j^-) = F_{X_j}(X_j - 1)`.
    For continuous variables the left limit and the cdf itself coincide. Respective columns can be
    omitted in the second block.
    
    Parameters
    ----------
    u :
        An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of evaluation points, where :math:`k`
        is the number of discrete variables (see ``Vinecop.select()``).
    
    N :
        Integer for the number of quasi-random numbers to draw to evaluate the distribution (default:
        1e4).
    
    num_threads :
        The number of threads to use for computations; if greater than 1, the function will generate
        ``n`` samples concurrently in ``num_threads`` batches.
    
    seeds :
        Seeds to scramble the quasi-random numbers; if empty (default), the random number
        quasi-generator is seeded randomly.
    
    Returns
    -------
    A vector of length ``n`` containing the copula distribution values.
    """
    ...

  @property
  def dim(self) -> Any: ...
  @property
  def families(self) -> Any: ...
  def fit(self, data: ArrayLike, controls: "FitControlsBicop" = ..., num_threads: int = 1) -> None:
    """
    
    Fits the parameters of a pre-specified vine copula model.
    
    This method fits the pair-copulas of a vine copula model. It is assumed that the structure and
    pair-copula families are already set. The method is equivalent to calling ``fit()`` for each
    pair-copula in the model. The same can be achieved by calling ``select()`` with the same data and a
    ``FitControlsVinecop`` object instantiated with ``select_families = false``.
    
    Parameters
    ----------
    data :
        :math:`n \times (d + k)` or :math:`n \times 2d` matrix of observations, where :math:`k` is the
        number of discrete variables.
    
    controls :
        The controls for each bivariate fit (see ``FitControlsBicop()``).
    
    num_threads :
        The number of threads to use for parallel computation.
    """
    ...

  def format(self, trees: collections.abc.Sequence[int] = []) -> str:
    """
    
    Summarizes the model into a string (can be used for printing).
    
    Parameters
    ----------
    trees :
        A vector of tree indices to summarize; if empty, all trees.
    """
    ...


  @staticmethod
  def from_data(data: ArrayLike, structure: "RVineStructure" | None = None, matrix: ArrayLike | None = None, var_types: collections.abc.Sequence[str] = [], controls: "FitControlsVinecop" = ...) -> "Vinecop":
    """
    
    Factory function to create a Vinecop from data.
    
    Parameters
    ----------
    data :
        Input data matrix.
    
    structure :
        An ``RVineStructure``. Provide either this or `matrix`, but not both.
    
    matrix :
        RVine matrix. Provide either this or `structure`, but not both.
    
    var_types :
        Variable types for each variable (e.g., 'c' for continuous, 'd' for discrete). Defaults to all continuous.
    
    controls :
        Fit controls for the vinecop. Defaults to the default constructor.
    """
    ...


  @staticmethod
  def from_dimension(d: int) -> "Vinecop":
    """
    
    Instantiates a D-vine with all pair-copulas set to independence.
    
    Parameters
    ----------
    d :
        The dimension (= number of variables) of the model.
    """
    ...


  @staticmethod
  def from_file(filename: str, check: bool = True) -> "Vinecop":
    """
    
    Instantiates from a JSON file.
    
    The input file contains 2 attributes : ``"structure"`` for the vine structure, which itself contains
    attributes ``"array"`` for the structure triangular array and ``"order"`` for the order vector, and
    ``"pair copulas"``. ``"pair copulas"`` contains a list of attributes for the trees (``"tree1"``,
    ``"tree2"``, etc), each containing a list of attributes for the edges (``"pc1"``, ``"pc2"``, etc).
    See the corresponding method of ``Bicop`` objects for the encoding of pair-copulas.
    
    Parameters
    ----------
    filename :
        The name of the JSON file to read.
    
    check :
        Whether to check if the ``"structure"`` node of the input represents a valid R-vine structure.
    """
    ...


  @staticmethod
  def from_json(json: str, check: bool = True) -> "Vinecop":
    """
    
    Instantiates from a nlohmann::json object.
    
    Parameters
    ----------
    input :
        The nlohmann::json object to convert from (see ``to_json()`` for the structure of the input).
    
    check :
        Whether to check if the ``"structure"`` node represents a valid R-vine structure.
    """
    ...


  @staticmethod
  def from_structure(structure: "RVineStructure" | None = None, matrix: ArrayLike | None = None, pair_copulas: collections.abc.Sequence[collections.abc.Sequence["Bicop"]] = [], var_types: collections.abc.Sequence[str] = []) -> "Vinecop":
    """
    
    Factory function to create a Vinecop using either a structure or a matrix.
    
    Parameters
    ----------
    structure :
        An ``RVineStructure``. Provide either this or `matrix`, but not both.
    
    matrix :
        Vinecop matrix. Provide either this or `structure`, but not both.
    
    pair_copulas :
        Pairwise copulas for each edge in the vine. Defaults to an empty list.
    
    var_types :
        Variable types for each variable (e.g., 'c' for continuous, 'd' for discrete). Defaults to all continuous.
    """
    ...

  def get_family(self, tree: int, edge: int) -> "BicopFamily":
    """
    
    Gets the family of a pair-copula.
    """
    ...

  def get_pair_copula(self, tree: int, edge: int) -> "Bicop":
    """
    
    Gets a pair-copula.
    """
    ...

  def get_parameters(self, tree: int, edge: int) -> ArrayLike:
    """
    
    Gets the parameters of a pair-copula.
    """
    ...

  def get_rotation(self, tree: int, edge: int) -> int:
    """
    
    Gets the rotation of a pair-copula.
    """
    ...

  def get_tau(self, tree: int, edge: int) -> float:
    """
    
    Gets the kendall's tau of a pair-copula.
    """
    ...

  def inverse_rosenblatt(self, u: ArrayLike, num_threads: int = 1) -> ArrayLike:
    """
    
    Evaluates the inverse Rosenblatt transform.
    
    The inverse Rosenblatt transform can be used for simulation: the function applied to independent
    uniform variates resembles simulated data from the vine copula model.
    
    If the problem is too large, it is split recursively into halves (w.r.t. :math:`n`, the number of
    observations). "Too large" means that the required memory will exceed 1 GB. An examplary
    configuration requiring less than 1 GB is :math:`n = 1000`, :math:`d = 200`.
    
    The Rosenblatt transform (Rosenblatt, 1952) :math:`U = T(V)` of a random vector :math:`V =
    (V_1,\ldots,V_d) ~ F` is defined as
    
    .. math:: U_1= F(V_1), U_{2} = F(V_{2}|V_1), \ldots, U_d =F(V_d|V_1,\ldots,V_{d-1}),
    
    where :math:`F(v_k|v_1,\ldots,v_{k-1})` is the conditional distribution of :math:`V_k` given
    :math:`V_1 \ldots, V_{k-1}, k = 2,\ldots,d`. The vector :math:`U = (U_1, \dots, U_d)` then contains
    independent standard uniform variables. The inverse operation
    
    .. math:: V_1 = F^{-1}(U_1), V_{2} = F^{-1}(U_2|U_1), \ldots, V_d =F^{-1}(U_d|U_1,\ldots,U_{d-1})
    
    can be used to simulate from a distribution. For any copula :math:`F`, if :math:`U` is a vector of
    independent random variables, :math:`V = T^{-1}(U)` has distribution :math:`F`.
    
    The formulas above assume a vine copula model with order :math:`d, \dots, 1`. More generally,
    ``Vinecop.rosenblatt()`` returns the variables
    
    .. math:: U_{M[d - j, j]}= F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0, 0]}),
    
    where :math:`M` is the structure matrix. Similarly, ``Vinecop.inverse_rosenblatt()`` computes
    
    .. math:: V_{M[d - j, j]}= F^{-1}(U_{M[d - j, j]} | U_{M[d - j - 1, j - 1]}, \dots, U_{M[0, 0]}).
    
    Parameters
    ----------
    u :
        An :math:`n \times d` matrix of evaluation points.
    
    num_threads :
        The number of threads to use for computations; if greater than 1, the function will be applied
        concurrently to ``num_threads`` batches of ``u``.
    
    Returns
    -------
    An :math:`n \times d` matrix of evaluations.
    """
    ...

  def loglik(self, u: ArrayLike = ..., num_threads: int = 1) -> float:
    """
    
    Evaluates the log-likelihood.
    
    The log-likelihood is defined as
    
    .. math:: \mathrm{loglik} = \sum_{i = 1}^n \log c(U_{1, i}, ..., U_{d, i}),
    
    where :math:`c` is the copula density, see ``Vinecop.pdf()``.
    
    Parameters
    ----------
    u :
        An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of evaluation points, where :math:`k`
        is the number of discrete variables (see ``select()`` or ``Vinecop.pdf()``).
    
    num_threads :
        The number of threads to use for computations; if greater than 1, the function will be applied
        concurrently to ``num_threads`` batches of ``u``.
    
    Returns
    -------
    The log-likelihood as a double.
    """
    ...

  @property
  def matrix(self) -> Any: ...
  def mbicv(self, u: ArrayLike = ..., psi0: float = 0.9, num_threads: int = 1) -> float:
    """
    
    Evaluates the modified Bayesian information criterion for vines (mBICV).
    
    The mBICV is defined as
    
    .. math:: \mathrm{mBICV} = -2\, \mathrm{loglik} + \log(n) p, - 2 * \sum_{t=1}^(d - 1) \{q_t
    \log(\psi_0^t) - (d - t - q_t) \log(1 -\psi_0^t)\},
    
    where :math:`\mathrm{loglik}` is the log-liklihood, :math:`p` is the (effective) number of
    parameters of the model, :math:`t` is the tree level, :math:`\psi_0` is the prior probability of
    having a non-independence copula in the first tree, and :math:`q_t` is the number of
    non-independence copulas in tree :math:`t`; The vBIC is a consistent model selection criterion for
    parametric sparse vine copula models when :math:`d = o(\sqrt{n \log n})`.
    
    Parameters
    ----------
    u :
        An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of evaluation points, where :math:`k`
        is the number of discrete variables (see ``Vinecop.select()`` or ``Vinecop.pdf()``).
    
    psi0 :
        Baseline prior probability of a non-independence copula.
    
    num_threads :
        The number of threads to use for computations; if greater than 1, the function will be applied
        concurrently to ``num_threads`` batches of ``u``.
    
    Returns
    -------
    The mBICV as a double.
    """
    ...

  @property
  def nobs(self) -> Any: ...
  @property
  def npars(self) -> Any: ...
  @property
  def order(self) -> Any: ...
  @property
  def pair_copulas(self) -> Any: ...
  @property
  def parameters(self) -> Any: ...
  def pdf(self, u: ArrayLike, num_threads: int = 1) -> ArrayLike:
    """
    
    Evaluates the copula density.
    
    The copula density is defined as joint density divided by marginal densities, irrespective of
    variable types.
    
    When at least one variable is discrete, two types of "observations" are required in ``u``: the first
    :math:`n \; x \; d` block contains realizations of :math:`F_{X_j}(X_j)`. The second :math:`n \; x \;
    d` block contains realizations of :math:`F_{X_j}(X_j^-)`. The minus indicates a left-sided limit of
    the cdf. For, e.g., an integer-valued variable, it holds :math:`F_{X_j}(X_j^-) = F_{X_j}(X_j - 1)`.
    For continuous variables the left limit and the cdf itself coincide. Respective columns can be
    omitted in the second block.
    
    Parameters
    ----------
    u :
        An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of evaluation points, where :math:`k`
        is the number of discrete variables (see ``Vinecop.select()``).
    
    num_threads :
        The number of threads to use for computations; if greater than 1, the function will be applied
        concurrently to ``num_threads`` batches of ``u``.
    
    Returns
    -------
    A vector of length ``n`` containing the copula density values.
    """
    ...

  def plot(self, tree: object | None = None, add_edge_labels: bool = True, layout: str = 'graphviz', vars_names: object | None = None) -> None:
    """
    
    Generates a plot for the Vinecop object.
    
    This method generates a plot of the vine copula structure. It can be used to visualize the tree structure of the vine copula.
    
    Parameters
    ----------
    tree: list[int] (default=None)
        The tree indice(s) to plot. If None, all trees are plotted.
    add_edge_labels: bool (default=True)
        Whether to add edge labels to the plot.
    layout: str (default="graphviz")
        The layout to use for plotting. Either "graphviz" or "spring_layout".
    vars_names: list[str] (default=None)
        The names of the variables for the vine model. If None, the indices are used.
    
    Returns
    -------
    Nothing, the function generates a plot and shows it using matplotlib.
    
    Usage
    -----
    .. code-block:: python
    
        import pyvinecopulib as pv
        import numpy as np
        np.random.seed(1234)
        u = np.random.uniform(0, 1, size=(20, 10))
        vc = vc = pv.Vinecop.from_data(u, controls=pv.FitControlsVinecop(family_set=[pv.BicopFamily.indep]))
        vc.plot(tree=[0, 1, 2]) # Plots the first three trees
        vars_names = ["X" + str(i) for i in range(10)]
        vc.plot(vars_names=vars_names) # Using variable names for the plot
    """
    ...

  def rosenblatt(self, u: ArrayLike, num_threads: int = 1, randomize_discrete: bool = True, seeds: collections.abc.Sequence[int] = []) -> ArrayLike:
    """
    
    Evaluates the Rosenblatt transform for a vine copula model.
    
    The Rosenblatt transform converts data from this model into independent uniform variates.
    
    The Rosenblatt transform (Rosenblatt, 1952) :math:`U = T(V)` of a random vector :math:`V =
    (V_1,\ldots,V_d) ~ F` is defined as
    
    .. math:: U_1= F(V_1), U_{2} = F(V_{2}|V_1), \ldots, U_d =F(V_d|V_1,\ldots,V_{d-1}),
    
    where :math:`F(v_k|v_1,\ldots,v_{k-1})` is the conditional distribution of :math:`V_k` given
    :math:`V_1 \ldots, V_{k-1}, k = 2,\ldots,d`. The vector :math:`U = (U_1, \dots, U_d)` then contains
    independent standard uniform variables. The inverse operation
    
    .. math:: V_1 = F^{-1}(U_1), V_{2} = F^{-1}(U_2|U_1), \ldots, V_d =F^{-1}(U_d|U_1,\ldots,U_{d-1})
    
    can be used to simulate from a distribution. For any copula :math:`F`, if :math:`U` is a vector of
    independent random variables, :math:`V = T^{-1}(U)` has distribution :math:`F`.
    
    The formulas above assume a vine copula model with order :math:`d, \dots, 1`. More generally,
    ``Vinecop.rosenblatt()`` returns the variables
    
    .. math:: U_{M[d - j, j]}= F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0, 0]}),
    
    where :math:`M` is the structure matrix. Similarly, ``Vinecop.inverse_rosenblatt()`` computes
    
    .. math:: V_{M[d - j, j]}= F^{-1}(U_{M[d - j, j]} | U_{M[d - j - 1, j - 1]}, \dots, U_{M[0, 0]}).
    
    If some variables have atoms, Brockwell (10.1016/j.spl.2007.02.008) proposed a simple randomization
    scheme to ensure that output is still independent uniform if the model is correct. The
    transformation reads
    
    .. math:: U_{M[d - j, j]}= W_{d - j} F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0,
    0]}) + (1 - W_{d - j}) F^-(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0, 0]}),
    
    where :math:`F^-` is the left limit of the conditional cdf and :math:`W_1, \dots, W_d` are are
    independent standard uniform random variables. This is used by default. If you are interested in the
    conditional probabilities
    
    .. math:: F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0, 0]}),
    
    set ``randomize_discrete = FALSE``.
    
    Parameters
    ----------
    u :
        An :math:`n \times d` matrix of evaluation points.
    
    num_threads :
        The number of threads to use for computations; if greater than 1, the function will be applied
        concurrently to ``num_threads`` batches of ``u``.
    
    randomize_discrete :
        Whether to randomize the transform for discrete variables; see Details.
    
    seeds :
        Seeds to scramble the quasi-random numbers; if empty (default), the random number
        quasi-generator is seeded randomly. Only relevant if there are discrete variables and
        ``randomize_discrete = TRUE``.
    
    Returns
    -------
    An :math:`n \times d` matrix of independent uniform variates.
    """
    ...

  @property
  def rotations(self) -> Any: ...
  def select(self, data: ArrayLike, controls: "FitControlsVinecop" = ...) -> None:
    """
    
    In other words, ``select()`` behaves differently depending on its current truncation level and the
    truncation level specified in the controls, respectively called ``trunc_lvl`` and
    ``controls.trunc_lvl`` in what follows. Essentially, ``controls.trunc_lvl`` defines the object's
    truncation level after calling ``select()``:
    
      - If ``controls.trunc_lvl <= trunc_lvl``, the families and parameters for
        all pairs in trees smaller or equal to ``controls.trunc_lvl``
        are selected, using the current structure.
      - If ``controls.trunc_lvl > trunc_lvl``, ``select()`` behaves as above for
        all trees that are smaller or equal to ``trunc_lvl``, and then it selects
        the structure for higher trees along with the families and parameters.
        This includes the case where ``trunc_lvl = 0``, namely where the
        structure is fully unspecified.
    
    Selection of the structure is performed using the algorithm of Dissmann, J. F., E. C. Brechmann, C.
    Czado, and D. Kurowicka (2013). *Selecting and estimating regular vine copulae and application to
    financial returns.* Computational Statistics & Data Analysis, 59 (1), 52-69. The dependence measure
    used to select trees (default: Kendall's tau) is corrected for ties (see the `wdm
    <https://github.com/tnagler/wdm>`_ library). The dependence measure can be changed using the
    ``controls.tree_criterion``, which can be set to ``"tau"``, ``"rho"`` or ``"hoeffd"``. Both Prim's
    (default: ``"mst_prim"``) and Kruskal's ()``"mst_kruskal"``) algorithms are available through
    ``controls.tree_algorithm`` for the maximum spanning tree selection. An alternative to the maximum
    spanning tree selection is to use random spanning trees, which can be selected using
    ``controls.tree_algorithm`` and come in two flavors, both using Wilson's algorithm loop erased
    random walks:
    
      - "random_weighted"` generates a random spanning tree with probability
        proportional to the product of the weights (i.e., the dependence) of
        the edges in the tree.
      - "random_unweighted"` generates a random spanning tree uniformly over all
        spanning trees satisfying the proximity condition.
    
    If the ``controls`` object has been instantiated with ``select_families = false``, then the method
    simply updates the parameters of the pair-copulas without selecting the families or the structure.
    In this case, this is equivalent to calling ``fit()`` for each pair-copula, albeit potentially in
    parallel if ``num_threads > 1``.
    
    When at least one variable is discrete, two types of "observations" are required: the first :math:`n
    \times d` block contains realizations of :math:`F_Y(Y), F_X(X)`; the second :math:`n \times d` block
    contains realizations of :math:`F_Y(Y^-), F_X(X^-), ...`. The minus indicates a left-sided limit of
    the cdf. For continuous variables the left limit and the cdf itself coincide. For, e.g., an
    integer-valued variable, it holds :math:`F_Y(Y^-) = F_Y(Y - 1)`. Continuous variables in the second
    block can be omitted.
    
    If there are missing data (i.e., NaN entries), incomplete observations are discarded before fitting
    a pair-copula. This is done on a pair-by-pair basis so that the maximal available information is
    used.
    
    Parameters
    ----------
    data :
        :math:`n \times (d + k)` or :math:`n \times 2d` matrix of observations, where :math:`k` is the
        number of discrete variables.
    
    controls :
        The controls to the algorithm (see ``FitControlsVinecop()``).
    """
    ...

  def simulate(self, n: int, qrng: bool = False, num_threads: int = 1, seeds: collections.abc.Sequence[int] = []) -> ArrayLike:
    """
    
    Simulates from a vine copula model, see ``inverse_rosenblatt()``.
    
    Simulated data is always a continous :math:`n \times d` matrix. Sampling from a vine copula model is
    done by first generating :math:`n \times d` uniform random numbers and then applying the inverse
    Rosenblatt transformation.
    
    Parameters
    ----------
    n :
        Number of observations.
    
    qrng :
        Set to true for quasi-random numbers.
    
    num_threads :
        The number of threads to use for computations; if greater than 1, the function will generate
        ``n`` samples concurrently in ``num_threads`` batches.
    
    seeds :
        Seeds of the random number generator; if empty (default), the random number generator is seeded
        randomly.
    
    Returns
    -------
    An :math:`n \times d` matrix of samples from the copula model.
    """
    ...

  @property
  def structure(self) -> Any: ...
  @property
  def taus(self) -> Any: ...
  @property
  def threshold(self) -> Any: ...
  def to_file(self, filename: str) -> None:
    """
    
    Writes the copula object into a JSON file.
    
    The output file contains 2 attributes : ``"structure"`` for the vine structure, which itself
    contains attributes ``"array"`` for the structure triangular array and ``"order"`` for the order
    vector, and ``"pair copulas"``. ``"pair copulas"`` contains a list of attributes for the trees
    (``"tree1"``, ``"tree2"``, etc), each containing a list of attributes for the edges (``"pc1"``,
    ``"pc2"``, etc). See ``Bicop.to_file()`` objects for the encoding of pair-copulas.
    
    Parameters
    ----------
    filename :
        The name of the JSON file to write.
    """
    ...

  def to_json(self) -> str:
    """
    
    Converts the copula into a nlohmann::json object.
    
    The JSON-like `str` object contains two nodes : ``"structure"`` for the vine structure, which
    itself contains nodes ``"array"`` for the structure triangular array and ``"order"`` for the order
    vector, and ``"pair copulas"``. The former two encode the R-Vine structure and the latter is a list
    of child nodes for the trees (``"tree1"``, ``"tree2"``, etc), each containing a list of child nodes
    for the edges (``"pc1"``, ``"pc2"``, etc). See Bicop::to_json() for the encoding of pair-copulas.
    
    Returns
    -------
    the nlohmann::json object containing the copula.
    """
    ...

  @property
  def trunc_lvl(self) -> Any: ...
  def truncate(self, trunc_lvl: int) -> None:
    """
    
    Truncates the vine copula model.
    
    While model for a ``d`` dimensional random vector contains at most ``d-1`` nested trees, this
    function extracts a sub-model based on a given truncation level.
    
    If the model is already truncated at a level less than ``trunc_lvl``, the function does nothing.
    
    Parameters
    ----------
    trunc_lvl :
        The truncation level.
    """
    ...

  @property
  def var_types(self) -> Any: ...
  @var_types.setter
  def var_types(self, value: Any) -> None: ...

__version__: str = ...

all: list[BicopFamily] = ...

archimedean: list[BicopFamily] = ...

bb: list[BicopFamily] = ...

bb1: BicopFamily = ...

bb6: BicopFamily = ...

bb7: BicopFamily = ...

bb8: BicopFamily = ...

def benchmark(data: ArrayLike) -> list[float]:
  ...

clayton: BicopFamily = ...

elliptical: list[BicopFamily] = ...

extreme_value: list[BicopFamily] = ...

frank: BicopFamily = ...

gaussian: BicopFamily = ...

def ghalton(n: int, d: int, seeds: collections.abc.Sequence[int] = []) -> ArrayLike:
  """
  
  Simulates from the multivariate Generalized Halton Sequence.
  
  For more information on Generalized Halton Sequence, see Faure, H., Lemieux, C. (2009). Generalized
  Halton Sequences in 2008: A Comparative Study. ACM-TOMACS 19(4), Article 15.
  
  Parameters
  ----------
  n :
      Number of observations.
  
  d :
      Dimension.
  
  seeds :
      Seeds to scramble the quasi-random numbers; if empty (default), the quasi-random number
      generator is seeded randomly.
  
  Returns
  -------
  An :math:`n \times d` matrix of quasi-random :math:`\mathrm{U}[0, 1]` variables.
  """
  ...

gumbel: BicopFamily = ...

indep: BicopFamily = ...

itau: list[BicopFamily] = ...

joe: BicopFamily = ...

lt: list[BicopFamily] = ...

nonparametric: list[BicopFamily] = ...

one_par: list[BicopFamily] = ...

parametric: list[BicopFamily] = ...

rotationless: list[BicopFamily] = ...

def simulate_uniform(n: int, d: int, qrng: bool = False, seeds: collections.abc.Sequence[int] = []) -> ArrayLike:
  """
  
  Simulates from the multivariate uniform distribution.
  
  If ``qrng = TRUE``, generalized Halton sequences (see ``ghalton()``) are used for :math:`d \leq 300`
  and Sobol sequences otherwise (see ``sobol()``).
  
  Parameters
  ----------
  n :
      Number of observations.
  
  d :
      Dimension.
  
  qrng :
      If true, quasi-numbers are generated.
  
  seeds :
      Seeds of the random number generator; if empty (default), the random number generator is seeded
      randomly.
  
  Returns
  -------
  An :math:`n \times d` matrix of independent :math:`\mathrm{U}[0, 1]` random variables.
  """
  ...

def sobol(n: int, d: int, seeds: collections.abc.Sequence[int] = []) -> ArrayLike:
  """
  
  Simulates from the multivariate Sobol sequence.
  
  For more information on the Sobol sequence, see S. Joe and F. Y. Kuo (2008), constructing Sobol
  sequences with better two-dimensional projections, SIAM J. Sci. Comput. 30, 2635–2654.
  
  Parameters
  ----------
  n :
      Number of observations.
  
  d :
      Dimension.
  
  seeds :
      Seeds to scramble the quasi-random numbers; if empty (default), the quasi-random number
      generator is seeded randomly.
  
  Returns
  -------
  An :math:`n \times d` matrix of quasi-random :math:`\mathrm{U}[0, 1]` variables.
  """
  ...

student: BicopFamily = ...

tawn: BicopFamily = ...

three_par: list[BicopFamily] = ...

tll: BicopFamily = ...

def to_pseudo_obs(x: ArrayLike, ties_method: str = 'average', weights: ArrayLike = ..., seeds: collections.abc.Sequence[int] = []) -> ArrayLike:
  """
  
  Applies the empirical probability integral transform to a data matrix.
  
  Gives pseudo-observations from the copula by applying the empirical distribution function (scaled by
  :math:`n + 1`) to each margin/column.
  
  Parameters
  ----------
  x :
      A matrix of real numbers.
  
  ties_method :
      Indicates how to treat ties; same as in R, see
      https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.
  
  weights :
      Vector of weights for the observations.
  
  Returns
  -------
  Pseudo-observations of the copula, i.e. :math:`F_X(x)` (column-wise).
  """
  ...

two_par: list[BicopFamily] = ...

ut: list[BicopFamily] = ...
