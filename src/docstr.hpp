#pragma once
// GENERATED FILE DO NOT EDIT
// This file contains docstrings for the Python bindings that were
// automatically extracted by mkdoc.py.
#if defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif
// #include "vinecopulib/bicop/class.hpp"
// #include "vinecopulib/bicop/family.hpp"
// #include "vinecopulib/bicop/fit_controls.hpp"
// #include "vinecopulib/bicop/implementation/class.ipp"
// #include "vinecopulib/bicop/implementation/family.ipp"
// #include "vinecopulib/bicop/implementation/fit_controls.ipp"
// #include "vinecopulib/misc/implementation/tools_stats.ipp"
// #include "vinecopulib/misc/tools_stats.hpp"
// #include "vinecopulib/vinecop/class.hpp"
// #include "vinecopulib/vinecop/fit_controls.hpp"
// #include "vinecopulib/vinecop/implementation/class.ipp"
// #include "vinecopulib/vinecop/implementation/fit_controls.ipp"
// #include "vinecopulib/vinecop/implementation/rvine_structure.ipp"
// #include "vinecopulib/vinecop/rvine_structure.hpp"

// Symbol: pyvinecopulib_doc
constexpr struct /* pyvinecopulib_doc */ {
  // Symbol: vinecopulib
  struct /* vinecopulib */ {
    // Symbol: vinecopulib::Bicop
    struct /* Bicop */ {
      // Source: vinecopulib/bicop/class.hpp:22
      const char* doc =
R"""(A class for bivariate copula models.

The copula model is fully characterized by the family, rotation, and
parameters.)""";
      // Symbol: vinecopulib::Bicop::Bicop
      struct /* ctor */ {
        // Source: vinecopulib/bicop/class.hpp:27
        const char* doc_4args_family_rotation_parameters_var_types =
R"""(Instantiates a specific bivariate copula model.

Parameter ``family``:
    The copula family.

Parameter ``rotation``:
    The rotation of the copula; one of 0, 90, 180, or 270 (for
    Independence, Gaussian, Student, Frank, and nonparametric
    families, only 0 is allowed).

Parameter ``parameters``:
    The copula parameters.

Parameter ``var_types``:
    Two strings specifying the types of the variables, e.g., ``("c",
    "d")`` means first variable continuous, second discrete.)""";
        // Source: vinecopulib/bicop/class.hpp:32
        const char* doc_3args_data_controls_var_types =
R"""(Instantiates from data.

Equivalent to creating a default ``Bicop()`` and then selecting the
model using ``select()``.

Parameter ``data``:
    See ``select()``.

Parameter ``controls``:
    See ``select()``.

Parameter ``var_types``:
    Two strings specifying the types of the variables, e.g., ``("c",
    "d")`` means first variable continuous, second discrete.)""";
        // Source: vinecopulib/bicop/class.hpp:36
        const char* doc_copy =
R"""(Copy constructor (deep copy)

Parameter ``other``:
    Bicop object to copy.)""";
        // Source: vinecopulib/bicop/class.hpp:38
        const char* doc_1args_filename =
R"""(Instantiates from a JSON file.

The input file contains four attributes: ``"family"``, ``"rotation"``,
``"parameters"``, ``"var_types"`` respectively a string for the family
name, an integer for the rotation, and a numeric matrix for the
parameters, and a list of two strings for the variable types.

Parameter ``filename``:
    The name of the JSON file to read.)""";
        // Source: vinecopulib/bicop/class.hpp:40
        const char* doc_1args_input =
R"""(Instantiates from a nlohmann::json object.

Parameter ``input``:
    The nlohmann::json object to convert from (see ``to_json()`` for
    the structure of the input).)""";
      } ctor;
      // Symbol: vinecopulib::Bicop::aic
      struct /* aic */ {
        // Source: vinecopulib/bicop/class.hpp:104
        const char* doc =
R"""(Evaluates the Akaike information criterion (AIC).

The AIC is defined as

.. math:: \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p,

where :math:`\mathrm{loglik}` is the log-liklihood (see ``loglik()``)
and :math:`p` is the (effective) number of parameters of the model.
The AIC is a consistent model selection criterion even for
nonparametric models.

Parameter ``u``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.)""";
      } aic;
      // Symbol: vinecopulib::Bicop::as_continuous
      struct /* as_continuous */ {
        // Source: vinecopulib/bicop/class.hpp:124
        const char* doc = R"""()""";
      } as_continuous;
      // Symbol: vinecopulib::Bicop::bic
      struct /* bic */ {
        // Source: vinecopulib/bicop/class.hpp:106
        const char* doc =
R"""(Evaluates the Bayesian information criterion (BIC).

The BIC is defined as

.. math:: \mathrm{BIC} = -2\, \mathrm{loglik} + \log(n) p,

where :math:`\mathrm{loglik}` is the log-liklihood (see ``loglik()``)
and :math:`p` is the (effective) number of parameters of the model.
The BIC is a consistent model selection criterion for parametric
models.

Parameter ``u``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.)""";
      } bic;
      // Symbol: vinecopulib::Bicop::cdf
      struct /* cdf */ {
        // Source: vinecopulib/bicop/class.hpp:79
        const char* doc =
R"""(Evaluates the copula distribution.

Parameter ``u``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.

Returns:
    The copula distribution evaluated at ``u``.)""";
      } cdf;
      // Symbol: vinecopulib::Bicop::check_data
      struct /* check_data */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:579
        const char* doc = R"""()""";
      } check_data;
      // Symbol: vinecopulib::Bicop::check_data_dim
      struct /* check_data_dim */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:586
        const char* doc = R"""()""";
      } check_data_dim;
      // Symbol: vinecopulib::Bicop::check_fitted
      struct /* check_fitted */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:971
        const char* doc =
R"""(Checks whether the Bicop object was fitted to data.)""";
      } check_fitted;
      // Symbol: vinecopulib::Bicop::check_rotation
      struct /* check_rotation */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:944
        const char* doc =
R"""(Checks whether the supplied rotation is valid (only 0, 90, 180, 270
allowd).)""";
      } check_rotation;
      // Symbol: vinecopulib::Bicop::check_var_types
      struct /* check_var_types */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:982
        const char* doc =
R"""(Checks whether var_types have the correct length and are either "c" or
"d".)""";
      } check_var_types;
      // Symbol: vinecopulib::Bicop::check_weights_size
      struct /* check_weights_size */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:961
        const char* doc =
R"""(Checks whether weights and data have matching sizes.)""";
      } check_weights_size;
      // Symbol: vinecopulib::Bicop::compute_mbic_penalty
      struct /* compute_mbic_penalty */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:550
        const char* doc = R"""()""";
      } compute_mbic_penalty;
      // Symbol: vinecopulib::Bicop::fit
      struct /* fit */ {
        // Source: vinecopulib/bicop/class.hpp:95
        const char* doc =
R"""(Fits a bivariate copula (with fixed family) to data.

For parametric models, two different methods are available. ``"mle"``
fits the parameters by maximum-likelihood. ``"itau"`` uses inversion
of Kendall's :math:`\tau`, but is only available for one-parameter
families and the Student t copula. For the latter, there is a
one-to-one transformation for the first parameter, the second is found
by profile likelihood optimization (with accuracy of at least 0.5).
Nonparametric families have specialized methods, no specification is
required.

When at least one variable is discrete, two types of "observations"
are required: the first :math:`n \times 2` block contains realizations
of :math:`F_{X_1}(X_1), F_{X_2}(X_2)`. Let :math:`k` denote the number
of discrete variables (either one or two). Then the second :math:`n
\times k` block contains realizations of :math:`F_{X_k}(X_k^-)`. The
minus indicates a left-sided limit of the cdf. For continuous
variables the left limit and the cdf itself coincide. For, e.g., an
integer-valued variable, it holds :math:`F_{X_k}(X_k^-) = F_{X_k}(X_k
- 1)`.

Incomplete observations (i.e., ones with a NaN value) are discarded.

Parameter ``data``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.

Parameter ``controls``:
    The controls (see FitControlsBicop).)""";
      } fit;
      // Symbol: vinecopulib::Bicop::flip
      struct /* flip */ {
        // Source: vinecopulib/bicop/class.hpp:118
        const char* doc =
R"""(Adjusts the copula model to a change in the variable order.)""";
      } flip;
      // Symbol: vinecopulib::Bicop::flip_abstract_var_types
      struct /* flip_abstract_var_types */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:608
        const char* doc = R"""()""";
      } flip_abstract_var_types;
      // Symbol: vinecopulib::Bicop::format_data
      struct /* format_data */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:868
        const char* doc =
R"""(Adds an additional column if there's only one discrete variable;
removes superfluous columns for continuous variables. (continuous
models only require two columns, discrete models always four))""";
      } format_data;
      // Symbol: vinecopulib::Bicop::get_aic
      struct /* get_aic */ {
        // Source: vinecopulib/bicop/class.hpp:64
        const char* doc = R"""(Gets the aic (only for fitted objects).)""";
      } get_aic;
      // Symbol: vinecopulib::Bicop::get_bic
      struct /* get_bic */ {
        // Source: vinecopulib/bicop/class.hpp:65
        const char* doc = R"""(Gets the bic (only for fitted objects).)""";
      } get_bic;
      // Symbol: vinecopulib::Bicop::get_bicop
      struct /* get_bicop */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:704
        const char* doc = R"""()""";
      } get_bicop;
      // Symbol: vinecopulib::Bicop::get_family
      struct /* get_family */ {
        // Source: vinecopulib/bicop/class.hpp:50
        const char* doc = R"""(Gets the copula family.)""";
      } get_family;
      // Symbol: vinecopulib::Bicop::get_family_name
      struct /* get_family_name */ {
        // Source: vinecopulib/bicop/class.hpp:52
        const char* doc = R"""(Gets the copula family as a string.)""";
      } get_family_name;
      // Symbol: vinecopulib::Bicop::get_loglik
      struct /* get_loglik */ {
        // Source: vinecopulib/bicop/class.hpp:62
        const char* doc =
R"""(Gets the log-likelihood (only for fitted objects).)""";
      } get_loglik;
      // Symbol: vinecopulib::Bicop::get_mbic
      struct /* get_mbic */ {
        // Source: vinecopulib/bicop/class.hpp:66
        const char* doc =
R"""(Gets the modified bic (only for fitted objects).)""";
      } get_mbic;
      // Symbol: vinecopulib::Bicop::get_n_discrete
      struct /* get_n_discrete */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:996
        const char* doc =
R"""(Returns the number of discrete variables.)""";
      } get_n_discrete;
      // Symbol: vinecopulib::Bicop::get_nobs
      struct /* get_nobs */ {
        // Source: vinecopulib/bicop/class.hpp:63
        const char* doc =
R"""(Gets the number of observations (only for fitted objects).)""";
      } get_nobs;
      // Symbol: vinecopulib::Bicop::get_npars
      struct /* get_npars */ {
        // Source: vinecopulib/bicop/class.hpp:60
        const char* doc =
R"""(Returns the actual number of parameters for parameteric families.

For nonparametric families, there is a conceptually similar definition
in the sense that it can be used in the calculation of fit statistics.)""";
      } get_npars;
      // Symbol: vinecopulib::Bicop::get_parameters
      struct /* get_parameters */ {
        // Source: vinecopulib/bicop/class.hpp:56
        const char* doc = R"""(Gets the parameters.)""";
      } get_parameters;
      // Symbol: vinecopulib::Bicop::get_parameters_lower_bounds
      struct /* get_parameters_lower_bounds */ {
        // Source: vinecopulib/bicop/class.hpp:120
        const char* doc =
R"""(Gets lower bounds for copula parameters.)""";
      } get_parameters_lower_bounds;
      // Symbol: vinecopulib::Bicop::get_parameters_upper_bounds
      struct /* get_parameters_upper_bounds */ {
        // Source: vinecopulib/bicop/class.hpp:122
        const char* doc =
R"""(Gets upper bounds for copula parameters.)""";
      } get_parameters_upper_bounds;
      // Symbol: vinecopulib::Bicop::get_rotation
      struct /* get_rotation */ {
        // Source: vinecopulib/bicop/class.hpp:54
        const char* doc = R"""(Gets the rotation.)""";
      } get_rotation;
      // Symbol: vinecopulib::Bicop::get_tau
      struct /* get_tau */ {
        // Source: vinecopulib/bicop/class.hpp:58
        const char* doc = R"""(Gets the Kendall's tau.)""";
      } get_tau;
      // Symbol: vinecopulib::Bicop::get_var_types
      struct /* get_var_types */ {
        // Source: vinecopulib/bicop/class.hpp:74
        const char* doc = R"""(Gets variable types.)""";
      } get_var_types;
      // Symbol: vinecopulib::Bicop::hfunc1
      struct /* hfunc1 */ {
        // Source: vinecopulib/bicop/class.hpp:81
        const char* doc =
R"""(Evaluates the first h-function.

The first h-function is :math:`h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 =
u_1)`.

Parameter ``u``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.)""";
      } hfunc1;
      // Symbol: vinecopulib::Bicop::hfunc2
      struct /* hfunc2 */ {
        // Source: vinecopulib/bicop/class.hpp:83
        const char* doc =
R"""(Evaluates the second h-function.

The second h-function is :math:`h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 =
u_2)`.

Parameter ``u``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.)""";
      } hfunc2;
      // Symbol: vinecopulib::Bicop::hinv1
      struct /* hinv1 */ {
        // Source: vinecopulib/bicop/class.hpp:85
        const char* doc =
R"""(Evaluates the inverse of the first h-function.

The first h-function is :math:`h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 =
u_1)`. The inverse is calulated w.r.t. the second argument.

Parameter ``u``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.)""";
      } hinv1;
      // Symbol: vinecopulib::Bicop::hinv2
      struct /* hinv2 */ {
        // Source: vinecopulib/bicop/class.hpp:87
        const char* doc =
R"""(Evaluates the inverse of the second h-function.

The second h-function is :math:`h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 =
u_2)`. The inverse is calculated w.r.t. the first argument.

Parameter ``u``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.)""";
      } hinv2;
      // Symbol: vinecopulib::Bicop::loglik
      struct /* loglik */ {
        // Source: vinecopulib/bicop/class.hpp:102
        const char* doc =
R"""(Evaluates the log-likelihood.

The log-likelihood is defined as

.. math:: \mathrm{loglik} = \sum_{i = 1}^n \log c(U_{1, i}, U_{2, i}),

where :math:`c` is the copula density, see ``Bicop::pdf()``.

Parameter ``u``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.)""";
      } loglik;
      // Symbol: vinecopulib::Bicop::mbic
      struct /* mbic */ {
        // Source: vinecopulib/bicop/class.hpp:108
        const char* doc =
R"""(Evaluates the modified Bayesian information criterion (mBIC).

The mBIC is defined as

.. math:: \mathrm{BIC} = -2\, \mathrm{loglik} + p \log(n) - 2 (I
\log(\psi_0) + (1 - I) \log(1 - \psi_0),

where :math:`\mathrm{loglik}` is the \log-liklihood (see
``loglik()``), :math:`p` is the (effective) number of parameters of
the model, and :math:`\psi_0` is the prior probability of having a
non-independence copula and :math:`I` is an indicator for the family
being non-independence.

Parameter ``u``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.

Parameter ``psi0``:
    Prior probability of a non-independence copula.)""";
      } mbic;
      // Symbol: vinecopulib::Bicop::parameters_to_tau
      struct /* parameters_to_tau */ {
        // Source: vinecopulib/bicop/class.hpp:114
        const char* doc =
R"""(Converts the copula parameters to Kendall's :math:`tau`.

Parameter ``parameters``:
    The parameters (must be a valid parametrization of the current
    family).)""";
      } parameters_to_tau;
      // Symbol: vinecopulib::Bicop::pdf
      struct /* pdf */ {
        // Source: vinecopulib/bicop/class.hpp:77
        const char* doc =
R"""(Evaluates the copula density.

The copula density is defined as joint density divided by marginal
densities, irrespective of variable types.

Parameter ``u``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.

Returns:
    The copula density evaluated at ``u``.)""";
      } pdf;
      // Symbol: vinecopulib::Bicop::prep_for_abstract
      struct /* prep_for_abstract */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:933
        const char* doc =
R"""(Prepares data for use with the ``AbstractBicop`` class: - add an
additional column if there's only one discrete variable. - trim the
data to the interval [1e-10, 1 - 1e-10] for numerical stability. -
rotate the data appropriately (``AbstractBicop`` is always
0deg-rotation).)""";
      } prep_for_abstract;
      // Symbol: vinecopulib::Bicop::rotate_data
      struct /* rotate_data */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:897
        const char* doc =
R"""(Rotates the data corresponding to the models rotation.

Parameter ``u``:
    An ``n x 2`` matrix.)""";
      } rotate_data;
      // Symbol: vinecopulib::Bicop::select
      struct /* select */ {
        // Source: vinecopulib/bicop/class.hpp:98
        const char* doc =
R"""(Selects the best fitting model.

The function calls ``fit()`` for all families in ``family_set`` and
selecting the best fitting model by either BIC or AIC, see ``bic()``
and ``aic()``.

When at least one variable is discrete, two types of "observations"
are required: the first :math:`n \times 2` block contains realizations
of :math:`F_{X_1}(X_1), F_{X_2}(X_2)`. Let :math:`k` denote the number
of discrete variables (either one or two). Then the second :math:`n
\times k` block contains realizations of :math:`F_{X_k}(X_k^-)`. The
minus indicates a left-sided limit of the cdf. For continuous
variables the left limit and the cdf itself coincide. For, e.g., an
integer-valued variable, it holds :math:`F_{X_k}(X_k^-) = F_{X_k}(X_k
- 1)`.

Incomplete observations (i.e., ones with a NaN value) are discarded.

Parameter ``data``:
    An :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)`, where :math:`k` is the number of discrete
    variables.

Parameter ``controls``:
    The controls (see FitControlsBicop).)""";
      } select;
      // Symbol: vinecopulib::Bicop::set_parameters
      struct /* set_parameters */ {
        // Source: vinecopulib/bicop/class.hpp:70
        const char* doc = R"""()""";
      } set_parameters;
      // Symbol: vinecopulib::Bicop::set_rotation
      struct /* set_rotation */ {
        // Source: vinecopulib/bicop/class.hpp:68
        const char* doc = R"""(Sets the rotation.)""";
      } set_rotation;
      // Symbol: vinecopulib::Bicop::set_var_types
      struct /* set_var_types */ {
        // Source: vinecopulib/bicop/class.hpp:72
        const char* doc =
R"""(Sets variable types.

Parameter ``var_types``:
    A vector of size two specifying the types of the variables, e.g.,
    ``{"c", "d"}`` means first variable continuous, second discrete.)""";
      } set_var_types;
      // Symbol: vinecopulib::Bicop::simulate
      struct /* simulate */ {
        // Source: vinecopulib/bicop/class.hpp:89
        const char* doc =
R"""(Simulates from a bivariate copula.

If ``qrng = TRUE``, generalized Halton sequences are used. For more
information on Generalized Halton sequences, see Faure, H., Lemieux,
C. (2009). Generalized Halton Sequences in 2008: A Comparative Study.
ACM-TOMACS 19(4), Article 15.

Parameter ``n``:
    Number of observations.

Parameter ``qrng``:
    Set to true for quasi-random numbers.

Parameter ``seeds``:
    Seeds of the (quasi-)random number generator; if empty (default),
    the (quasi-)random number generator is seeded randomly.

Returns:
    An :math:`n \times 2` matrix of samples from the copula model.)""";
      } simulate;
      // Symbol: vinecopulib::Bicop::str
      struct /* str */ {
        // Source: vinecopulib/bicop/class.hpp:112
        const char* doc =
R"""(Summarizes the model into a string (can be used for printing).)""";
      } str;
      // Symbol: vinecopulib::Bicop::tau_to_parameters
      struct /* tau_to_parameters */ {
        // Source: vinecopulib/bicop/class.hpp:116
        const char* doc =
R"""(Converts a Kendall's :math:`\tau` into copula parameters.

It only works for one-parameter families.

Parameter ``tau``:
    A value in :math:`(-1, 1)`.)""";
      } tau_to_parameters;
      // Symbol: vinecopulib::Bicop::to_file
      struct /* to_file */ {
        // Source: vinecopulib/bicop/class.hpp:47
        const char* doc =
R"""(Write the copula object into a JSON file.

The written file contains four attributes: ``"family"``,
``"rotation"``, ``"parameters"``, ``"var_types"`` respectively a
string for the family name, an integer for the rotation, and a numeric
matrix for the parameters, and a list of two strings for the variable
types.

Parameter ``filename``:
    The name of the file to write.)""";
      } to_file;
      // Symbol: vinecopulib::Bicop::to_json
      struct /* to_json */ {
        // Source: vinecopulib/bicop/class.hpp:45
        const char* doc =
R"""(Convert the copula into a nlohmann::json object.

The nlohmann::json is contains of three values named ``"family"``,
``"rotation"``, ``"parameters"``, ``"var_types"``, respectively a
string for the family name, an integer for the rotation, a numeric
matrix for the parameters and a list of two strings for the variables
types.

Returns:
    the nlohmann::json object containing the copula.)""";
      } to_json;
    } Bicop;
    // Symbol: vinecopulib::BicopFamily
    struct /* BicopFamily */ {
      // Source: vinecopulib/bicop/family.hpp:15
      const char* doc = R"""(A bivariate copula family identifier.)""";
      // Symbol: vinecopulib::BicopFamily::bb1
      struct /* bb1 */ {
        // Source: vinecopulib/bicop/family.hpp:24
        const char* doc = R"""(BB1 copula)""";
      } bb1;
      // Symbol: vinecopulib::BicopFamily::bb6
      struct /* bb6 */ {
        // Source: vinecopulib/bicop/family.hpp:25
        const char* doc = R"""(BB6 copula)""";
      } bb6;
      // Symbol: vinecopulib::BicopFamily::bb7
      struct /* bb7 */ {
        // Source: vinecopulib/bicop/family.hpp:26
        const char* doc = R"""(BB7 copula)""";
      } bb7;
      // Symbol: vinecopulib::BicopFamily::bb8
      struct /* bb8 */ {
        // Source: vinecopulib/bicop/family.hpp:27
        const char* doc = R"""(BB8 copula)""";
      } bb8;
      // Symbol: vinecopulib::BicopFamily::clayton
      struct /* clayton */ {
        // Source: vinecopulib/bicop/family.hpp:20
        const char* doc = R"""(Clayton copula)""";
      } clayton;
      // Symbol: vinecopulib::BicopFamily::frank
      struct /* frank */ {
        // Source: vinecopulib/bicop/family.hpp:22
        const char* doc = R"""(Frank copula)""";
      } frank;
      // Symbol: vinecopulib::BicopFamily::gaussian
      struct /* gaussian */ {
        // Source: vinecopulib/bicop/family.hpp:18
        const char* doc = R"""(Gaussian copula)""";
      } gaussian;
      // Symbol: vinecopulib::BicopFamily::gumbel
      struct /* gumbel */ {
        // Source: vinecopulib/bicop/family.hpp:21
        const char* doc = R"""(Gumbel copula)""";
      } gumbel;
      // Symbol: vinecopulib::BicopFamily::indep
      struct /* indep */ {
        // Source: vinecopulib/bicop/family.hpp:17
        const char* doc = R"""(Independence copula)""";
      } indep;
      // Symbol: vinecopulib::BicopFamily::joe
      struct /* joe */ {
        // Source: vinecopulib/bicop/family.hpp:23
        const char* doc = R"""(Joe copula)""";
      } joe;
      // Symbol: vinecopulib::BicopFamily::student
      struct /* student */ {
        // Source: vinecopulib/bicop/family.hpp:19
        const char* doc = R"""(Student t copula)""";
      } student;
      // Symbol: vinecopulib::BicopFamily::tll
      struct /* tll */ {
        // Source: vinecopulib/bicop/family.hpp:28
        const char* doc =
R"""(Transformation local likelihood kernel estimator)""";
      } tll;
    } BicopFamily;
    // Symbol: vinecopulib::BicopPtr
    struct /* BicopPtr */ {
      // Source: vinecopulib/bicop/class.hpp:16
      const char* doc =
R"""(A shared pointer to an object of class AbstracBicop.)""";
    } BicopPtr;
    // Symbol: vinecopulib::CVineStructure
    struct /* CVineStructure */ {
      // Source: vinecopulib/vinecop/rvine_structure.hpp:178
      const char* doc =
R"""(A class for C-vine structures.

C-vines are a special class of R-vines where each tree is a star. A
C-vine structure is determined entirely by the order of variables. For
example, if the order is ``{1, 2, 3, 4}``, the first tree in the vine
connects variable 4 with all others, the second tree connects variable
3 with all others, etc.

Note that ``CVineStructure`` objects inherit the methods and
attributes of ``RVineStructure`` objects.)""";
      // Symbol: vinecopulib::CVineStructure::CVineStructure
      struct /* ctor */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:743
        const char* doc_1args =
R"""(Parameter ``order``:
    The order of variables in the C-vine (diagonal entries in the
    R-vine array); must be a permutation of 1, ..., d.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:753
        const char* doc_2args =
R"""(Parameter ``order``:
    The order of variables in the C-vine (diagonal entries in the
    R-vine array); must be a permutation of 1, ..., d.

Parameter ``trunc_lvl``:
    The truncation level.)""";
      } ctor;
    } CVineStructure;
    // Symbol: vinecopulib::DVineStructure
    struct /* DVineStructure */ {
      // Source: vinecopulib/vinecop/rvine_structure.hpp:161
      const char* doc =
R"""(A class for D-vine structures.

D-vines are a special class of R-vines where each tree is a path. A
D-vine structure is determined entirely by the order of variables. For
example, if the order is ``(1, 2, 3, 4)``, the first tree in the vine
is 1-2-3-4 and all further trees are unique due to the proximity
condition.

Note that ``DVineStructure`` objects inherit the methods and
attributes of ``RVineStructure`` objects.)""";
      // Symbol: vinecopulib::DVineStructure::DVineStructure
      struct /* ctor */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:723
        const char* doc_1args =
R"""(Parameter ``order``:
    The order of variables in the D-vine (diagonal entries in the
    R-vine array); must be a permutation of 1, ..., d.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:733
        const char* doc_2args =
R"""(Parameter ``order``:
    The order of variables in the D-vine (diagonal entries in the
    R-vine array); must be a permutation of 1, ..., d.

Parameter ``trunc_lvl``:
    The truncation level.)""";
      } ctor;
    } DVineStructure;
    // Symbol: vinecopulib::FitControlsBicop
    struct /* FitControlsBicop */ {
      // Source: vinecopulib/bicop/fit_controls.hpp:16
      const char* doc =
R"""(A class for controlling fits of bivariate copula models.)""";
      // Symbol: vinecopulib::FitControlsBicop::FitControlsBicop
      struct /* ctor */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:20
        const char* doc_9args =
R"""(Instantiates the controls for fitting bivariate copula models.

Parameter ``family_set``:
    The set of copula families to consider (if empty, then all
    families are included).

Parameter ``parametric_method``:
    The fit method for parametric families; possible choices:
    ``"mle"``, ``"itau"``.

Parameter ``nonparametric_method``:
    The fit method for the local-likelihood nonparametric family
    (TLLs); possible choices: ``"constant"``, ``"linear"``,
    ``"quadratic"``.

Parameter ``nonparametric_mult``:
    A factor with which the smoothing parameters are multiplied.

Parameter ``selection_criterion``:
    The selection criterion (``"loglik"``, ``"aic"`` or ``"bic"``).

Parameter ``weights``:
    A vector of weights for the observations.

Parameter ``psi0``:
    Only for `selection_criterion = "mbic", the prior probability of
    non-independence.

Parameter ``preselect_families``:
    Whether to exclude families before fitting based on symmetry
    properties of the data.

Parameter ``num_threads``:
    Number of concurrent threads to use while fitting copulas for
    different families; never uses more than the number of concurrent
    threads supported by the implementation.)""";
        // Source: vinecopulib/bicop/fit_controls.hpp:30
        const char* doc_1args =
R"""(Instantiates default controls except for the parameteric method.

Parameter ``parametric_method``:
    The fit method for parametric families; possible choices:
    ``"mle"``, ``"itau"``.)""";
        // Source: vinecopulib/bicop/fit_controls.hpp:32
        const char* doc_2args =
R"""(Instantiates default controls except for the nonparametric method.

Parameter ``nonparametric_method``:
    The fit method for the local-likelihood nonparametric family
    (TLLs); possible choices: ``"constant"``, ``"linear"``,
    ``"quadratic"``.

Parameter ``nonparametric_mult``:
    A factor with which the smoothing parameters are multiplied.)""";
      } ctor;
      // Symbol: vinecopulib::FitControlsBicop::check_nonparametric_method
      struct /* check_nonparametric_method */ {
        // Source: vinecopulib/bicop/implementation/fit_controls.ipp:89
        const char* doc = R"""()""";
      } check_nonparametric_method;
      // Symbol: vinecopulib::FitControlsBicop::check_nonparametric_mult
      struct /* check_nonparametric_mult */ {
        // Source: vinecopulib/bicop/implementation/fit_controls.ipp:99
        const char* doc = R"""()""";
      } check_nonparametric_mult;
      // Symbol: vinecopulib::FitControlsBicop::check_parametric_method
      struct /* check_parametric_method */ {
        // Source: vinecopulib/bicop/implementation/fit_controls.ipp:81
        const char* doc = R"""(@name Sanity checks)""";
      } check_parametric_method;
      // Symbol: vinecopulib::FitControlsBicop::check_psi0
      struct /* check_psi0 */ {
        // Source: vinecopulib/bicop/implementation/fit_controls.ipp:119
        const char* doc = R"""()""";
      } check_psi0;
      // Symbol: vinecopulib::FitControlsBicop::check_selection_criterion
      struct /* check_selection_criterion */ {
        // Source: vinecopulib/bicop/implementation/fit_controls.ipp:107
        const char* doc = R"""()""";
      } check_selection_criterion;
      // Symbol: vinecopulib::FitControlsBicop::get_family_set
      struct /* get_family_set */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:36
        const char* doc = R"""(returns the family set.)""";
      } get_family_set;
      // Symbol: vinecopulib::FitControlsBicop::get_nonparametric_method
      struct /* get_nonparametric_method */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:40
        const char* doc = R"""(returns the nonparametric method.)""";
      } get_nonparametric_method;
      // Symbol: vinecopulib::FitControlsBicop::get_nonparametric_mult
      struct /* get_nonparametric_mult */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:42
        const char* doc =
R"""(returns the nonparametric bandwidth multiplier.)""";
      } get_nonparametric_mult;
      // Symbol: vinecopulib::FitControlsBicop::get_num_threads
      struct /* get_num_threads */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:52
        const char* doc = R"""(returns the number of threads.)""";
      } get_num_threads;
      // Symbol: vinecopulib::FitControlsBicop::get_parametric_method
      struct /* get_parametric_method */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:38
        const char* doc = R"""(returns the parametric method.)""";
      } get_parametric_method;
      // Symbol: vinecopulib::FitControlsBicop::get_preselect_families
      struct /* get_preselect_families */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:48
        const char* doc = R"""(returns whether to preselect families.)""";
      } get_preselect_families;
      // Symbol: vinecopulib::FitControlsBicop::get_psi0
      struct /* get_psi0 */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:50
        const char* doc =
R"""(returns the baseline probability for mBIC selection.)""";
      } get_psi0;
      // Symbol: vinecopulib::FitControlsBicop::get_selection_criterion
      struct /* get_selection_criterion */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:44
        const char* doc = R"""()""";
      } get_selection_criterion;
      // Symbol: vinecopulib::FitControlsBicop::get_weights
      struct /* get_weights */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:46
        const char* doc = R"""(returns the observation weights.)""";
      } get_weights;
      // Symbol: vinecopulib::FitControlsBicop::process_num_threads
      struct /* process_num_threads */ {
        // Source: vinecopulib/bicop/implementation/fit_controls.ipp:262
        const char* doc = R"""()""";
      } process_num_threads;
      // Symbol: vinecopulib::FitControlsBicop::set_family_set
      struct /* set_family_set */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:55
        const char* doc = R"""(Sets the family set.)""";
      } set_family_set;
      // Symbol: vinecopulib::FitControlsBicop::set_nonparametric_method
      struct /* set_nonparametric_method */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:59
        const char* doc = R"""(Sets the nonparmetric method.)""";
      } set_nonparametric_method;
      // Symbol: vinecopulib::FitControlsBicop::set_nonparametric_mult
      struct /* set_nonparametric_mult */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:61
        const char* doc = R"""(Sets the nonparametric multiplier.)""";
      } set_nonparametric_mult;
      // Symbol: vinecopulib::FitControlsBicop::set_num_threads
      struct /* set_num_threads */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:71
        const char* doc = R"""(Sets the number of threads.)""";
      } set_num_threads;
      // Symbol: vinecopulib::FitControlsBicop::set_parametric_method
      struct /* set_parametric_method */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:57
        const char* doc = R"""(Sets the parametric method.)""";
      } set_parametric_method;
      // Symbol: vinecopulib::FitControlsBicop::set_preselect_families
      struct /* set_preselect_families */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:67
        const char* doc = R"""(Sets whether to preselect the families.)""";
      } set_preselect_families;
      // Symbol: vinecopulib::FitControlsBicop::set_psi0
      struct /* set_psi0 */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:69
        const char* doc = R"""(Sets the prior probability for mBIC.)""";
      } set_psi0;
      // Symbol: vinecopulib::FitControlsBicop::set_selection_criterion
      struct /* set_selection_criterion */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:63
        const char* doc = R"""(Sets the selection criterion.)""";
      } set_selection_criterion;
      // Symbol: vinecopulib::FitControlsBicop::set_weights
      struct /* set_weights */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:65
        const char* doc = R"""(Sets the observation weights.)""";
      } set_weights;
      // Symbol: vinecopulib::FitControlsBicop::str
      struct /* str */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:74
        const char* doc =
R"""(Summarizes the controls into a string (can be used for printing).)""";
      } str;
      // Symbol: vinecopulib::FitControlsBicop::str_internal
      struct /* str_internal */ {
        // Source: vinecopulib/bicop/fit_controls.hpp:77
        const char* doc = R"""()""";
      } str_internal;
    } FitControlsBicop;
    // Symbol: vinecopulib::FitControlsVinecop
    struct /* FitControlsVinecop */ {
      // Source: vinecopulib/vinecop/fit_controls.hpp:24
      const char* doc =
R"""(A class for controlling fits of vine copula models.)""";
      // Symbol: vinecopulib::FitControlsVinecop::FitControlsVinecop
      struct /* ctor */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:28
        const char* doc_0args =
R"""(Instantiates default controls for fitting vine copula models.)""";
        // Source: vinecopulib/vinecop/fit_controls.hpp:47
        const char* doc_8args =
R"""(Instantiates custom controls for fitting vine copula models.

Parameter ``trunc_lvl``:
    Truncation level for truncated vines.

Parameter ``tree_criterion``:
    The criterion for selecting the maximum spanning tree ("tau",
    "hoeffd" and "rho" implemented so far).

Parameter ``threshold``:
    For thresholded vines (0 = no threshold).

Parameter ``show_trace``:
    Whether to show a trace of the building progress.

Parameter ``select_trunc_lvl``:
    Whether the truncation shall be selected automatically.

Parameter ``select_threshold``:
    Whether the threshold parameter shall be selected automatically.

Parameter ``controls``:
    See FitControlsBicop.

Parameter ``num_threads``:
    Number of concurrent threads to use while fitting pair copulas
    within a tree; never uses more than the number returned by
    `std::thread::hardware_concurrency()``.)""";
        // Source: vinecopulib/vinecop/implementation/fit_controls.ipp:56
        const char* doc_15args =
R"""(Instantiates custom controls for fitting vine copula models.

Parameter ``family_set``:
    The set of copula families to consider (if empty, then all
    families are included).

Parameter ``parametric_method``:
    The fit method for parametric families; possible choices:
    ``"mle"``, ``"itau"``.

Parameter ``nonparametric_method``:
    The fit method for the local-likelihood nonparametric family
    (TLLs); possible choices: ``"constant"``, ``"linear"``,
    ``"quadratic"``.

Parameter ``nonparametric_mult``:
    A factor with which the smoothing parameters are multiplied.

Parameter ``trunc_lvl``:
    Truncation level for truncated vines.

Parameter ``tree_criterion``:
    The criterion for selecting the maximum spanning tree ("tau",
    "hoeffd", "rho", and "mcor" implemented so far).

Parameter ``threshold``:
    For thresholded vines (0 = no threshold).

Parameter ``selection_criterion``:
    The selection criterion (``"loglik"``, ``"aic"`` or ``"bic"``).

Parameter ``weights``:
    A vector of weights for the observations.

Parameter ``psi0``:
    Only for `selection_criterion = "mbic", prior probability of
    non-independence.

Parameter ``preselect_families``:
    Whether to exclude families before fitting based on symmetry
    properties of the data.

Parameter ``select_trunc_lvl``:
    Whether the truncation shall be selected automatically.

Parameter ``select_threshold``:
    Whether the threshold parameter shall be selected automatically.

Parameter ``show_trace``:
    Whether to show a trace of the building progress.

Parameter ``num_threads``:
    Number of concurrent threads to use while fitting pair copulas
    within a tree; never uses more than the number of concurrent
    threads supported by the implementation.)""";
      } ctor;
      // Symbol: vinecopulib::FitControlsVinecop::check_threshold
      struct /* check_threshold */ {
        // Source: vinecopulib/vinecop/implementation/fit_controls.ipp:136
        const char* doc = R"""()""";
      } check_threshold;
      // Symbol: vinecopulib::FitControlsVinecop::check_tree_criterion
      struct /* check_tree_criterion */ {
        // Source: vinecopulib/vinecop/implementation/fit_controls.ipp:126
        const char* doc = R"""(@name Sanity checks)""";
      } check_tree_criterion;
      // Symbol: vinecopulib::FitControlsVinecop::get_fit_controls_bicop
      struct /* get_fit_controls_bicop */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:74
        const char* doc =
R"""(Returns the fit controls for bivariate fitting.)""";
      } get_fit_controls_bicop;
      // Symbol: vinecopulib::FitControlsVinecop::get_select_threshold
      struct /* get_select_threshold */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:70
        const char* doc =
R"""(returns whether to select the threshold automatically.)""";
      } get_select_threshold;
      // Symbol: vinecopulib::FitControlsVinecop::get_select_trunc_lvl
      struct /* get_select_trunc_lvl */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:68
        const char* doc =
R"""(returns whether to select the truncation level automatically.)""";
      } get_select_trunc_lvl;
      // Symbol: vinecopulib::FitControlsVinecop::get_select_truncation_level
      struct /* get_select_truncation_level */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:67
        const char* doc = R"""()""";
      } get_select_truncation_level;
      // Symbol: vinecopulib::FitControlsVinecop::get_show_trace
      struct /* get_show_trace */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:65
        const char* doc =
R"""(returns whether to show a trace is during fitting.)""";
      } get_show_trace;
      // Symbol: vinecopulib::FitControlsVinecop::get_threshold
      struct /* get_threshold */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:63
        const char* doc = R"""(returns the threshold parameter.)""";
      } get_threshold;
      // Symbol: vinecopulib::FitControlsVinecop::get_tree_criterion
      struct /* get_tree_criterion */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:61
        const char* doc =
R"""(returns the criterion for tree selection.)""";
      } get_tree_criterion;
      // Symbol: vinecopulib::FitControlsVinecop::get_trunc_lvl
      struct /* get_trunc_lvl */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:59
        const char* doc = R"""(returns the truncation level.)""";
      } get_trunc_lvl;
      // Symbol: vinecopulib::FitControlsVinecop::get_truncation_level
      struct /* get_truncation_level */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:58
        const char* doc = R"""()""";
      } get_truncation_level;
      // Symbol: vinecopulib::FitControlsVinecop::needs_sparse_select
      struct /* needs_sparse_select */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:72
        const char* doc = R"""()""";
      } needs_sparse_select;
      // Symbol: vinecopulib::FitControlsVinecop::set_fit_controls_bicop
      struct /* set_fit_controls_bicop */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:91
        const char* doc =
R"""(Sets the fit controls for bivariate fitting.)""";
      } set_fit_controls_bicop;
      // Symbol: vinecopulib::FitControlsVinecop::set_select_threshold
      struct /* set_select_threshold */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:89
        const char* doc =
R"""(Sets whether to select the threshold automatically.)""";
      } set_select_threshold;
      // Symbol: vinecopulib::FitControlsVinecop::set_select_trunc_lvl
      struct /* set_select_trunc_lvl */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:87
        const char* doc =
R"""(Sets whether to select the truncation level automatically.)""";
      } set_select_trunc_lvl;
      // Symbol: vinecopulib::FitControlsVinecop::set_select_truncation_level
      struct /* set_select_truncation_level */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:86
        const char* doc = R"""()""";
      } set_select_truncation_level;
      // Symbol: vinecopulib::FitControlsVinecop::set_show_trace
      struct /* set_show_trace */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:84
        const char* doc =
R"""(returns whether to show a trace is during fitting.)""";
      } set_show_trace;
      // Symbol: vinecopulib::FitControlsVinecop::set_threshold
      struct /* set_threshold */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:82
        const char* doc = R"""(Sets the threshold parameter.)""";
      } set_threshold;
      // Symbol: vinecopulib::FitControlsVinecop::set_tree_criterion
      struct /* set_tree_criterion */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:80
        const char* doc = R"""(Sets the criterion for tree selection.)""";
      } set_tree_criterion;
      // Symbol: vinecopulib::FitControlsVinecop::set_trunc_lvl
      struct /* set_trunc_lvl */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:78
        const char* doc = R"""(Sets the truncation level.)""";
      } set_trunc_lvl;
      // Symbol: vinecopulib::FitControlsVinecop::set_truncation_level
      struct /* set_truncation_level */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:77
        const char* doc = R"""()""";
      } set_truncation_level;
      // Symbol: vinecopulib::FitControlsVinecop::str
      struct /* str */ {
        // Source: vinecopulib/vinecop/fit_controls.hpp:94
        const char* doc =
R"""(Summarizes the controls into a string (can be used for printing).)""";
      } str;
    } FitControlsVinecop;
    // Symbol: vinecopulib::RVineStructure
    struct /* RVineStructure */ {
      // Source: vinecopulib/vinecop/rvine_structure.hpp:67
      const char* doc =
R"""(A class for R-vine structures.

RVineStructure objects encode the tree structure of the vine, i.e. the
conditioned/conditioning variables of each edge. It is represented by
a triangular array. An exemplary array is


::

    4 4 4 4
    3 3 3
    2 2
    1

which encodes the following pair-copulas:


::

    | tree | edge | pair-copulas   |
    |------|------|----------------|
    | 0    | 0    | ``(1, 4)``       |
    |      | 1    | ``(2, 4)``       |
    |      | 2    | ``(3, 4)``       |
    | 1    | 0    | ``(1, 3; 4)``    |
    |      | 1    | ``(2, 3; 4)``    |
    | 2    | 0    | ``(1, 2; 3, 4)`` |

Denoting by ``M[i, j]`` the array entry in row ``i`` and column ``j``,
the pair-copula index for edge ``e`` in tree ``t`` of a ``d``
dimensional vine is ``(M[d - 1 - e, e], M[t, e]; M[t - 1, e], ...,
M[0, e])``. Less formally, 1. Start with the counter-diagonal element
of column ``e`` (first conditioned variable). 2. Jump up to the
element in row ``t`` (second conditioned variable). 3. Gather all
entries further up in column ``e`` (conditioning set).

A valid R-vine array must satisfy several conditions which are checked
when ``RVineStructure()`` is called: 1. It only contains numbers
between 1 and d. 2. The diagonal must contain the numbers 1, ..., d.
3. The diagonal entry of a column must not be contained in any column
further to the right. 4. The entries of a column must be contained in
all columns to the left. 5. The proximity condition must hold: For all
t = 1, ..., d - 2 and e = 0, ..., d - t - 1 there must exist an index
j > d, such that ``(M[t, e], {M[0, e], ..., M[t-1, e]})`` equals
either ``(M[d-j-1, j], {M[0, j], ..., M[t-1, j]})`` or ``(M[t-1, j],
{M[d-j-1, j], M[0, j], ..., M[t-2, j]})``.

An R-vine array is said to be in natural order when the anti-diagonal
entries are :math:`1, \dots, d` (from left to right). The exemplary
arrray above is in natural order. Any R-vine array can be
characterized by the diagonal entries (called order) and the entries
below the diagonal of the corresponding R-vine array in natural order.
Since most algorithms work with the structure in natural order, this
is how RVineStructure stores the structure internally.)""";
      // Symbol: vinecopulib::RVineStructure::RVineStructure
      struct /* ctor */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:28
        const char* doc_2args_mat_check =
R"""(Instantiates an RVineStructure object from a matrix representing an
R-vine array.

The matrix must contain zeros in the lower right triangle and the
upper left triangle must be a valid R-vine array. Truncated vines can
be encoded by putting zeros above the digonal in all rows below the
truncation level. Example of a 1-truncated matrix:


::

    4 4 4 4
    0 0 3 0
    0 2 0 0
    1 0 0 0

Parameter ``mat``:
    A matrix representing a valid R-vine array.

Parameter ``check``:
    Whether ``mat`` shall be checked for validity.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:64
        const char* doc_2args_d_trunc_lvl =
R"""(Instantiates an RVineStructure object to a D-vine for a given
dimension.

Parameter ``d``:
    The dimension.

Parameter ``trunc_lvl``:
    The truncation level. By default, it is dim - 1.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:74
        const char* doc_3args_order_trunc_lvl_check =
R"""(Instantiates an RVineStructure object to a D-vine with a given
ordering of the variables.

Parameter ``order``:
    The order of variables in the D-vine (diagonal entries in the
    R-vine array); must be a permutation of 1, ..., d.

Parameter ``trunc_lvl``:
    The truncation level. By default, it is d - 1.

Parameter ``check``:
    Whether `order shall be checked for validity.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:100
        const char* doc_4args_order_struct_array_natural_order_check =
R"""(Instantiates an RVineStructure object from the variable order
(diagonal elements of the R-vine array) and a triangular structure
array (all elements above the diagonal).

Parameter ``order``:
    The order of variables (diagonal entries in the R-vine array);
    must be a permutation of 1, ..., d.

Parameter ``struct_array``:
    The structure array (all elements above the diagonal in the R-vine
    array). For truncated vines, all rows below the truncation level
    are omitted.

Parameter ``natural_order``:
    Whether ``struct_array`` is already in natural order.

Parameter ``check``:
    Whether ``order`` and ``struct_array`` shall be checked for
    validity.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:145
        const char* doc_2args_input_check =
R"""(Instantiates from a nlohmann::json object.

Parameter ``input``:
    The nlohmann::json object to convert from (see ``to_json()`` for
    the structure of the input).

Parameter ``check``:
    Whether to check if the input represents a valid R-vine structure.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:161
        const char* doc_2args_filename_check =
R"""(Instantiates an RVineStructure from a JSON file.

The file needs to contain two values: ``"array"`` for the structure
triangular array and ``"order"`` for the order vector.

Parameter ``filename``:
    The name of the JSON file to read.

Parameter ``check``:
    Whether to check if the input represents a valid R-vine matrix.)""";
      } ctor;
      // Symbol: vinecopulib::RVineStructure::check_antidiagonal
      struct /* check_antidiagonal */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:668
        const char* doc = R"""()""";
      } check_antidiagonal;
      // Symbol: vinecopulib::RVineStructure::check_columns
      struct /* check_columns */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:641
        const char* doc = R"""()""";
      } check_columns;
      // Symbol: vinecopulib::RVineStructure::check_if_quadratic
      struct /* check_if_quadratic */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:601
        const char* doc = R"""()""";
      } check_if_quadratic;
      // Symbol: vinecopulib::RVineStructure::check_lower_tri
      struct /* check_lower_tri */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:611
        const char* doc = R"""()""";
      } check_lower_tri;
      // Symbol: vinecopulib::RVineStructure::check_proximity_condition
      struct /* check_proximity_condition */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:679
        const char* doc = R"""()""";
      } check_proximity_condition;
      // Symbol: vinecopulib::RVineStructure::check_upper_tri
      struct /* check_upper_tri */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:625
        const char* doc = R"""()""";
      } check_upper_tri;
      // Symbol: vinecopulib::RVineStructure::compute_min_array
      struct /* compute_min_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:551
        const char* doc = R"""()""";
      } compute_min_array;
      // Symbol: vinecopulib::RVineStructure::compute_needed_hfunc1
      struct /* compute_needed_hfunc1 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:564
        const char* doc = R"""()""";
      } compute_needed_hfunc1;
      // Symbol: vinecopulib::RVineStructure::compute_needed_hfunc2
      struct /* compute_needed_hfunc2 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:582
        const char* doc = R"""()""";
      } compute_needed_hfunc2;
      // Symbol: vinecopulib::RVineStructure::d_
      struct /* d_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:140
        const char* doc = R"""()""";
      } d_;
      // Symbol: vinecopulib::RVineStructure::find_trunc_lvl
      struct /* find_trunc_lvl */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:448
        const char* doc =
R"""(Find the truncation level in an R-vine array.

The truncation level is determined by the first row (starting from the
bottom) that contains only zeros above the diagonal.

Parameter ``mat``:
    An array representing the R-vine array.)""";
      } find_trunc_lvl;
      // Symbol: vinecopulib::RVineStructure::get_dim
      struct /* get_dim */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:199
        const char* doc = R"""(Gets the dimension of the vine.)""";
      } get_dim;
      // Symbol: vinecopulib::RVineStructure::get_matrix
      struct /* get_matrix */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:425
        const char* doc = R"""(Gets the R-vine matrix representation.)""";
      } get_matrix;
      // Symbol: vinecopulib::RVineStructure::get_min_array
      struct /* get_min_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:247
        const char* doc =
R"""(Gets the minimum array.

The minimum array is derived from an R-vine array by iteratively
computing the (elementwise) minimum of two subsequent rows (starting
from the top). It is used in estimation and evaluation algorithms to
find the two edges in the previous tree that are joined by the current
edge.)""";
      } get_min_array;
      // Symbol: vinecopulib::RVineStructure::get_needed_hfunc1
      struct /* get_needed_hfunc1 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:257
        const char* doc =
R"""(Gets an array indicating which of the first h-functions are needed.

It is usually not necessary to compute both h-functions for each
pair-copula.)""";
      } get_needed_hfunc1;
      // Symbol: vinecopulib::RVineStructure::get_needed_hfunc2
      struct /* get_needed_hfunc2 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:267
        const char* doc =
R"""(Gets an array indicating which of the second h-functions are needed.

It is usually not necessary to compute both h-functions for each
pair-copula.)""";
      } get_needed_hfunc2;
      // Symbol: vinecopulib::RVineStructure::get_order
      struct /* get_order */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:214
        const char* doc_0args =
R"""(Extract the order of variables in the vine (diagonal entries in the
R-vine array).)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:476
        const char* doc_1args =
R"""(Find the order of an R-vine array.

The order is contained in the counter-diagonal of the R-vine array.

Parameter ``mat``:
    A matrix representing the R-vine array.)""";
      } get_order;
      // Symbol: vinecopulib::RVineStructure::get_struct_array
      struct /* get_struct_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:223
        const char* doc =
R"""(Extract structure array (all elements above the diagonal in the R-vine
array).

Parameter ``natural_order``:
    Whether indices correspond to natural order.)""";
      } get_struct_array;
      // Symbol: vinecopulib::RVineStructure::get_trunc_lvl
      struct /* get_trunc_lvl */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:206
        const char* doc = R"""(Gets the truncation level of the vine.)""";
      } get_trunc_lvl;
      // Symbol: vinecopulib::RVineStructure::make_cvine_struct_array
      struct /* make_cvine_struct_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:538
        const char* doc =
R"""(Creates a structure array corresponding to a D-vine in natural order.)""";
      } make_cvine_struct_array;
      // Symbol: vinecopulib::RVineStructure::make_dvine_struct_array
      struct /* make_dvine_struct_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:524
        const char* doc =
R"""(Creates a structure array corresponding to a D-vine in natural order.)""";
      } make_dvine_struct_array;
      // Symbol: vinecopulib::RVineStructure::min_array
      struct /* min_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:290
        const char* doc =
R"""(Access elements of the minimum array.

Parameter ``tree``:
    Tree index.

Parameter ``edge``:
    Edge index.)""";
      } min_array;
      // Symbol: vinecopulib::RVineStructure::min_array_
      struct /* min_array_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:143
        const char* doc = R"""()""";
      } min_array_;
      // Symbol: vinecopulib::RVineStructure::needed_hfunc1
      struct /* needed_hfunc1 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:299
        const char* doc =
R"""(Access elements of the needed_hfunc1 array.

Parameter ``tree``:
    Tree index.

Parameter ``edge``:
    Edge index.)""";
      } needed_hfunc1;
      // Symbol: vinecopulib::RVineStructure::needed_hfunc1_
      struct /* needed_hfunc1_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:145
        const char* doc = R"""()""";
      } needed_hfunc1_;
      // Symbol: vinecopulib::RVineStructure::needed_hfunc2
      struct /* needed_hfunc2 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:306
        const char* doc =
R"""(Access elements of the needed_hfunc2 array.)""";
      } needed_hfunc2;
      // Symbol: vinecopulib::RVineStructure::needed_hfunc2_
      struct /* needed_hfunc2_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:146
        const char* doc = R"""()""";
      } needed_hfunc2_;
      // Symbol: vinecopulib::RVineStructure::order_
      struct /* order_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:139
        const char* doc = R"""()""";
      } order_;
      // Symbol: vinecopulib::RVineStructure::simulate
      struct /* simulate */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:358
        const char* doc =
R"""(Randomly sample a regular vine structure.

Parameter ``d``:
    The dimension.

Parameter ``natural_order``:
    Should the sampled structure be in natural order?

Parameter ``seeds``:
    Seeds of the random number generator; if empty (default), the
    random number generator is seeded randomly.

Note:
    Implementation of Algorithm 13 in Harry Joe's 2014 book (p. 288),
    but there's a typo: the end of line 6 in the book should be
    'column j' instead of 'column k'.)""";
      } simulate;
      // Symbol: vinecopulib::RVineStructure::str
      struct /* str */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:331
        const char* doc =
R"""(Converts the structure to a string representation (most useful for
printing).)""";
      } str;
      // Symbol: vinecopulib::RVineStructure::struct_array
      struct /* struct_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:277
        const char* doc =
R"""(Accesses elements of the structure array.

Parameter ``tree``:
    Tree index.

Parameter ``edge``:
    Edge index.

Parameter ``natural_order``:
    Whether indices correspond to natural order.)""";
      } struct_array;
      // Symbol: vinecopulib::RVineStructure::struct_array_
      struct /* struct_array_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:142
        const char* doc = R"""()""";
      } struct_array_;
      // Symbol: vinecopulib::RVineStructure::to_file
      struct /* to_file */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:192
        const char* doc =
R"""(Write the structure into a JSON file.

The written file contains two values: ``"array"`` for the structure
triangular array and ``"order"`` for the order vector.

Parameter ``filename``:
    The name of the file to write.)""";
      } to_file;
      // Symbol: vinecopulib::RVineStructure::to_json
      struct /* to_json */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:173
        const char* doc =
R"""(Converts the structure into a nlohmann::json object.

The ``nlohmann::json`` object contains two nodes: ``"array"`` for the
structure triangular array and ``"order"`` for the order vector.

Returns:
    the nlohmann::json object containing the structure.)""";
      } to_json;
      // Symbol: vinecopulib::RVineStructure::to_natural_order
      struct /* to_natural_order */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:506
        const char* doc =
R"""(Converts ``struct_array_`` to natural order.)""";
      } to_natural_order;
      // Symbol: vinecopulib::RVineStructure::to_rvine_array
      struct /* to_rvine_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:490
        const char* doc =
R"""(Gets the structure array (entries above the diagonal in R-vine.
array).

Parameter ``mat``:
    A array representing the R-vine array.)""";
      } to_rvine_array;
      // Symbol: vinecopulib::RVineStructure::trunc_lvl_
      struct /* trunc_lvl_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:141
        const char* doc = R"""()""";
      } trunc_lvl_;
      // Symbol: vinecopulib::RVineStructure::truncate
      struct /* truncate */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:317
        const char* doc =
R"""(Truncates the R-vine structure.

Parameter ``trunc_lvl``:
    The truncation level.

If the structure is already truncated at a level less than
``trunc_lvl``, the function does nothing.)""";
      } truncate;
    } RVineStructure;
    // Symbol: vinecopulib::Vinecop
    struct /* Vinecop */ {
      // Source: vinecopulib/vinecop/class.hpp:25
      const char* doc =
R"""(A class for vine copula models.

A vine copula model is characterized by its structure (see
``RVineStructure`` objects) and the pair-copulas.)""";
      // Symbol: vinecopulib::Vinecop::Vinecop
      struct /* ctor */ {
        // Source: vinecopulib/vinecop/class.hpp:29
        const char* doc_0args = R"""()""";
        // Source: vinecopulib/vinecop/class.hpp:31
        const char* doc_1args_d =
R"""(Instantiates a D-vine with all pair-copulas set to independence.

Parameter ``d``:
    The dimension (= number of variables) of the model.)""";
        // Source: vinecopulib/vinecop/class.hpp:57
        const char* doc_2args_filename_check =
R"""(Instantiates from a JSON file.

The input file contains 2 attributes : ``"structure"`` for the vine
structure, which itself contains attributes ``"array"`` for the
structure triangular array and ``"order"`` for the order vector, and
``"pair copulas"``. ``"pair copulas"`` contains a list of attributes
for the trees (``"tree1"``, ``"tree2"``, etc), each containing a list
of attributes for the edges (``"pc1"``, ``"pc2"``, etc). See the
corresponding method of ``Bicop`` objects for the encoding of
pair-copulas.

Parameter ``filename``:
    The name of the JSON file to read.

Parameter ``check``:
    Whether to check if the ``"structure"`` node of the input
    represents a valid R-vine structure.)""";
        // Source: vinecopulib/vinecop/class.hpp:58
        const char* doc_2args_input_check =
R"""(Instantiates from a nlohmann::json object.

Parameter ``input``:
    The nlohmann::json object to convert from (see ``to_json()`` for
    the structure of the input).

Parameter ``check``:
    Whether to check if the ``"structure"`` node represents a valid
    R-vine structure.)""";
        // Source: vinecopulib/vinecop/implementation/class.ipp:32
        const char* doc_3args_structure_pair_copulas_var_types =
R"""(Instantiates an arbitrary vine copula model.

Parameter ``structure``:
    An RVineStructure object specifying the vine structure.

Parameter ``pair_copulas``:
    Bicop objects specifying the pair-copulas, namely a nested list
    such that ``pc_store[t][e]`` contains a ``Bicop`` object for the
    pair copula corresponding to tree ``t`` and edge ``e``.

Parameter ``var_types``:
    Strings specifying the types of the variables, e.g., ``("c",
    "d")`` means first variable continuous, second discrete. If empty,
    then all variables are set as continuous.)""";
        // Source: vinecopulib/vinecop/implementation/class.ipp:58
        const char* doc_3args_matrix_pair_copulas_var_types =
R"""(Instantiates an arbitrary vine copula model.

Parameter ``matrix``:
    An R-vine matrix specifying the vine structure.

Parameter ``pair_copulas``:
    Bicop objects specifying the pair-copulas, namely a nested list
    such that ``pc_store[t][e]`` contains a ``Bicop`` object for the
    pair copula corresponding to tree ``t`` and edge ``e``.

Parameter ``var_types``:
    Strings specifying the types of the variables, e.g., ``("c",
    "d")`` means first variable continuous, second discrete. If empty,
    then all variables are set as continuous.)""";
        // Source: vinecopulib/vinecop/implementation/class.ipp:77
        const char* doc_4args_data_structure_var_types_controls =
R"""(Instantiates from data.

Equivalent to creating a default ``Vinecop()`` and then selecting the
model using ``select()``.

Parameter ``data``:
    An :math:`n \times d` matrix of observations.

Parameter ``structure``:
    An RVineStructure object specifying the vine structure. If empty,
    then it is selected as part of the fit.

Parameter ``var_types``:
    Strings specifying the types of the variables, e.g., ``("c",
    "d")`` means first variable continuous, second discrete. If empty,
    then all variables are set as continuous.

Parameter ``controls``:
    See ``FitControlsVinecop()``.)""";
        // Source: vinecopulib/vinecop/implementation/class.ipp:116
        const char* doc_4args_data_matrix_var_types_controls =
R"""(Instantiates from data.

Equivalent to creating a default ``Vinecop()`` and then selecting the
model using ``select()``.

Parameter ``data``:
    An :math:`n \times d` matrix of observations.

Parameter ``matrix``:
    Either an empty matrix (default) or an R-vine structure matrix,
    see ``select()``. If empty, then it is selected as part of the
    fit.

Parameter ``var_types``:
    Strings specifying the types of the variables, e.g., ``("c",
    "d")`` means first variable continuous, second discrete. If empty,
    then all variables are set as continuous.

Parameter ``controls``:
    See ``FitControlsVinecop()``.)""";
      } ctor;
      // Symbol: vinecopulib::Vinecop::aic
      struct /* aic */ {
        // Source: vinecopulib/vinecop/class.hpp:150
        const char* doc =
R"""(Evaluates the Akaike information criterion (AIC).

The AIC is defined as

.. math:: \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p,

where :math:`\mathrm{loglik}` is the log-liklihood (see ``loglik()``)
and :math:`p` is the (effective) number of parameters of the model.
The AIC is a consistent model selection criterion even for
nonparametric models.

Parameter ``u``:
    An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``select()``).

Parameter ``num_threads``:
    The number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } aic;
      // Symbol: vinecopulib::Vinecop::bic
      struct /* bic */ {
        // Source: vinecopulib/vinecop/class.hpp:153
        const char* doc =
R"""(Evaluates the Bayesian information criterion (BIC).

The BIC is defined as

.. math:: \mathrm{BIC} = -2\, \mathrm{loglik} + \log(n) p,

where :math:`\mathrm{loglik}` is the log-liklihood (see ``loglik()``)
and :math:`p` is the (effective) number of parameters of the model.
The BIC is a consistent model selection criterion for nonparametric
models.

Parameter ``u``:
    An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``select()``).

Parameter ``num_threads``:
    The number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } bic;
      // Symbol: vinecopulib::Vinecop::calculate_mbicv_penalty
      struct /* calculate_mbicv_penalty */ {
        // Source: vinecopulib/vinecop/class.hpp:181
        const char* doc = R"""(Computes the penalty term for mBICV.)""";
      } calculate_mbicv_penalty;
      // Symbol: vinecopulib::Vinecop::cdf
      struct /* cdf */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:859
        const char* doc =
R"""(Evaluates the copula distribution.

Because no closed-form expression is available, the distribution is
estimated numerically using Monte Carlo integration.

Parameter ``u``:
    An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``select()``).

Parameter ``N``:
    Integer for the number of quasi-random numbers to draw to evaluate
    the distribution (default: 1e4).

Parameter ``num_threads``:
    The number of threads to use for computations; if greater than 1,
    the function will generate ``n`` samples concurrently in
    ``num_threads`` batches.

Parameter ``seeds``:
    Seeds to scramble the quasi-random numbers; if empty (default),
    the random number quasi-generator is seeded randomly.)""";
      } cdf;
      // Symbol: vinecopulib::Vinecop::check_data
      struct /* check_data */ {
        // Source: vinecopulib/vinecop/class.hpp:178
        const char* doc =
R"""(Checks if dimension d of the data matches the dimension of the vine.)""";
      } check_data;
      // Symbol: vinecopulib::Vinecop::check_data_dim
      struct /* check_data_dim */ {
        // Source: vinecopulib/vinecop/class.hpp:177
        const char* doc =
R"""(Checks if dimension d of the data matches the dimension of the vine.)""";
      } check_data_dim;
      // Symbol: vinecopulib::Vinecop::check_enough_data
      struct /* check_enough_data */ {
        // Source: vinecopulib/vinecop/class.hpp:185
        const char* doc = R"""(Checks if data size is large enough.)""";
      } check_enough_data;
      // Symbol: vinecopulib::Vinecop::check_fitted
      struct /* check_fitted */ {
        // Source: vinecopulib/vinecop/class.hpp:186
        const char* doc = R"""()""";
      } check_fitted;
      // Symbol: vinecopulib::Vinecop::check_indices
      struct /* check_indices */ {
        // Source: vinecopulib/vinecop/class.hpp:187
        const char* doc = R"""()""";
      } check_indices;
      // Symbol: vinecopulib::Vinecop::check_pair_copulas_rvine_structure
      struct /* check_pair_copulas_rvine_structure */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:1246
        const char* doc =
R"""(Checks if pair copulas are compatible with the R-vine structure.)""";
      } check_pair_copulas_rvine_structure;
      // Symbol: vinecopulib::Vinecop::check_var_types
      struct /* check_var_types */ {
        // Source: vinecopulib/vinecop/class.hpp:188
        const char* doc = R"""()""";
      } check_var_types;
      // Symbol: vinecopulib::Vinecop::check_weights_size
      struct /* check_weights_size */ {
        // Source: vinecopulib/vinecop/class.hpp:183
        const char* doc =
R"""(Checks if weights are compatible with the data.)""";
      } check_weights_size;
      // Symbol: vinecopulib::Vinecop::collapse_data
      struct /* collapse_data */ {
        // Source: vinecopulib/vinecop/class.hpp:192
        const char* doc =
R"""(Removes superfluous columns for continuous data.)""";
      } collapse_data;
      // Symbol: vinecopulib::Vinecop::d_
      struct /* d_ */ {
        // Source: vinecopulib/vinecop/class.hpp:169
        const char* doc = R"""()""";
      } d_;
      // Symbol: vinecopulib::Vinecop::finalize_fit
      struct /* finalize_fit */ {
        // Source: vinecopulib/vinecop/class.hpp:182
        const char* doc = R"""()""";
      } finalize_fit;
      // Symbol: vinecopulib::Vinecop::get_aic
      struct /* get_aic */ {
        // Source: vinecopulib/vinecop/class.hpp:115
        const char* doc =
R"""(Gets the AIC.

The function throws an error if model has not been fitted to data.)""";
      } get_aic;
      // Symbol: vinecopulib::Vinecop::get_all_families
      struct /* get_all_families */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:404
        const char* doc =
R"""(Gets the families of all pair copulas.

Returns:
    a nested std::vector with entry ``[t][e]`` corresponding to edge
    ``e`` in tree ``t``.)""";
      } get_all_families;
      // Symbol: vinecopulib::Vinecop::get_all_pair_copulas
      struct /* get_all_pair_copulas */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:380
        const char* doc =
R"""(Gets all pair copulas.

Returns:
    a nested std::vector with entry ``[t][e]`` corresponding to edge
    ``e`` in tree ``t``.)""";
      } get_all_pair_copulas;
      // Symbol: vinecopulib::Vinecop::get_all_parameters
      struct /* get_all_parameters */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:488
        const char* doc =
R"""(Gets the parameters of all pair copulas.

Returns:
    a nested std::vector with entry ``[t][e]`` corresponding to edge
    ``e`` in tree ``t``.)""";
      } get_all_parameters;
      // Symbol: vinecopulib::Vinecop::get_all_rotations
      struct /* get_all_rotations */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:436
        const char* doc =
R"""(Gets the rotations of all pair copulas.

Returns:
    a nested std::vector with entry ``[t][e]`` corresponding to edge
    ``e`` in tree ``t``.)""";
      } get_all_rotations;
      // Symbol: vinecopulib::Vinecop::get_all_taus
      struct /* get_all_taus */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:506
        const char* doc =
R"""(Gets the Kendall's :math:`tau`s of all pair copulas.

Returns:
    a nested std::vector with entry ``[t][e]`` corresponding to edge
    ``e`` in tree ``t``.)""";
      } get_all_taus;
      // Symbol: vinecopulib::Vinecop::get_bic
      struct /* get_bic */ {
        // Source: vinecopulib/vinecop/class.hpp:116
        const char* doc =
R"""(Gets the BIC.

The function throws an error if model has not been fitted to data.)""";
      } get_bic;
      // Symbol: vinecopulib::Vinecop::get_dim
      struct /* get_dim */ {
        // Source: vinecopulib/vinecop/class.hpp:101
        const char* doc =
R"""(Gets the dimension of the vine copula model.)""";
      } get_dim;
      // Symbol: vinecopulib::Vinecop::get_family
      struct /* get_family */ {
        // Source: vinecopulib/vinecop/class.hpp:79
        const char* doc =
R"""(Gets the family of a pair copula.

Parameter ``tree``:
    Tree index (starting with 0).

Parameter ``edge``:
    Edge index (starting with 0).)""";
      } get_family;
      // Symbol: vinecopulib::Vinecop::get_loglik
      struct /* get_loglik */ {
        // Source: vinecopulib/vinecop/class.hpp:113
        const char* doc =
R"""(Gets the log-likelihood (throws an error if model has not been. fitted
to data).)""";
      } get_loglik;
      // Symbol: vinecopulib::Vinecop::get_matrix
      struct /* get_matrix */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:542
        const char* doc =
R"""(Gets the structure matrix of the vine copula model.)""";
      } get_matrix;
      // Symbol: vinecopulib::Vinecop::get_mbicv
      struct /* get_mbicv */ {
        // Source: vinecopulib/vinecop/class.hpp:117
        const char* doc =
R"""(Gets the log-likelihood.

The function throws an error if model has not been fitted to data.)""";
      } get_mbicv;
      // Symbol: vinecopulib::Vinecop::get_n_discrete
      struct /* get_n_discrete */ {
        // Source: vinecopulib/vinecop/class.hpp:191
        const char* doc =
R"""(Returns the number of discrete variables.)""";
      } get_n_discrete;
      // Symbol: vinecopulib::Vinecop::get_nobs
      struct /* get_nobs */ {
        // Source: vinecopulib/vinecop/class.hpp:114
        const char* doc =
R"""(Gets the number of observations used for the fit.

The function throws an error if model has not been fitted to data.)""";
      } get_nobs;
      // Symbol: vinecopulib::Vinecop::get_npars
      struct /* get_npars */ {
        // Source: vinecopulib/vinecop/class.hpp:145
        const char* doc =
R"""(Returns sum of the number of parameters for all pair copulas (see.
Bicop::get_npars()).)""";
      } get_npars;
      // Symbol: vinecopulib::Vinecop::get_order
      struct /* get_order */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:528
        const char* doc =
R"""(Gets the order vector of the vine copula model.)""";
      } get_order;
      // Symbol: vinecopulib::Vinecop::get_pair_copula
      struct /* get_pair_copula */ {
        // Source: vinecopulib/vinecop/class.hpp:77
        const char* doc =
R"""(Gets a pair copula.

Parameter ``tree``:
    Tree index (starting with 0).

Parameter ``edge``:
    Edge index (starting with 0).)""";
      } get_pair_copula;
      // Symbol: vinecopulib::Vinecop::get_parameters
      struct /* get_parameters */ {
        // Source: vinecopulib/vinecop/class.hpp:83
        const char* doc =
R"""(Gets the parameters of a pair copula.

Parameter ``tree``:
    Tree index (starting with 0).

Parameter ``edge``:
    Edge index (starting with 0).)""";
      } get_parameters;
      // Symbol: vinecopulib::Vinecop::get_rotation
      struct /* get_rotation */ {
        // Source: vinecopulib/vinecop/class.hpp:81
        const char* doc =
R"""(Gets the rotation of a pair copula.

Parameter ``tree``:
    Tree index (starting with 0).

Parameter ``edge``:
    Edge index (starting with 0).)""";
      } get_rotation;
      // Symbol: vinecopulib::Vinecop::get_rvine_structure
      struct /* get_rvine_structure */ {
        // Source: vinecopulib/vinecop/class.hpp:105
        const char* doc =
R"""(Gets the structure matrix of the vine copula model.)""";
      } get_rvine_structure;
      // Symbol: vinecopulib::Vinecop::get_struct_array
      struct /* get_struct_array */ {
        // Source: vinecopulib/vinecop/class.hpp:109
        const char* doc =
R"""(Gets the above diagonal coefficients of the vine copula model.

Parameter ``natural_order``:
    Whether indices correspond to natural order.)""";
      } get_struct_array;
      // Symbol: vinecopulib::Vinecop::get_tau
      struct /* get_tau */ {
        // Source: vinecopulib/vinecop/class.hpp:85
        const char* doc =
R"""(Gets the Kendall's :math:`tau` of a pair copula.

Parameter ``tree``:
    Tree index (starting with 0).

Parameter ``edge``:
    Edge index (starting with 0).)""";
      } get_tau;
      // Symbol: vinecopulib::Vinecop::get_threshold
      struct /* get_threshold */ {
        // Source: vinecopulib/vinecop/class.hpp:112
        const char* doc =
R"""(Gets the threshold.

Usually zero except ``select_threshold == TRUE`` in
``FitControlsVinecop()``).)""";
      } get_threshold;
      // Symbol: vinecopulib::Vinecop::get_trunc_lvl
      struct /* get_trunc_lvl */ {
        // Source: vinecopulib/vinecop/class.hpp:87
        const char* doc = R"""()""";
      } get_trunc_lvl;
      // Symbol: vinecopulib::Vinecop::get_var_types
      struct /* get_var_types */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:729
        const char* doc = R"""(Gets the variable types.)""";
      } get_var_types;
      // Symbol: vinecopulib::Vinecop::inverse_rosenblatt
      struct /* inverse_rosenblatt */ {
        // Source: vinecopulib/vinecop/class.hpp:135
        const char* doc =
R"""(Evaluates the inverse Rosenblatt transform.

The inverse Rosenblatt transform can be used for simulation: the
function applied to independent uniform variates resembles simulated
data from the vine copula model.

If the problem is too large, it is split recursively into halves
(w.r.t. :math:`n`, the number of observations). "Too large" means that
the required memory will exceed 1 GB. An examplary configuration
requiring less than 1 GB is :math:`n = 1000`, :math:`d = 200`.

Only works for continous models.

Parameter ``u``:
    An :math:`n \times d` matrix of evaluation points.

Parameter ``num_threads``:
    The number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } inverse_rosenblatt;
      // Symbol: vinecopulib::Vinecop::loglik
      struct /* loglik */ {
        // Source: vinecopulib/vinecop/class.hpp:147
        const char* doc =
R"""(Evaluates the log-likelihood.

The log-likelihood is defined as

.. math:: \mathrm{loglik} = \sum_{i = 1}^n \log c(U_{1, i}, ..., U_{d,
i}),

where :math:`c` is the copula density, see ``Vinecop::pdf()``.

Parameter ``u``:
    An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``select()``).

Parameter ``num_threads``:
    The number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } loglik;
      // Symbol: vinecopulib::Vinecop::loglik_
      struct /* loglik_ */ {
        // Source: vinecopulib/vinecop/class.hpp:173
        const char* doc = R"""()""";
      } loglik_;
      // Symbol: vinecopulib::Vinecop::make_pair_copula_store
      struct /* make_pair_copula_store */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:243
        const char* doc =
R"""(Initializes object for storing pair copulas.

Parameter ``d``:
    Dimension of the vine copula.

Parameter ``trunc_lvl``:
    A truncation level (optional).

Returns:
    A nested list such that ``pc_store[t][e]`` contains a Bicop.
    object for the pair copula corresponding to tree ``t`` and edge
    ``e``.)""";
      } make_pair_copula_store;
      // Symbol: vinecopulib::Vinecop::mbicv
      struct /* mbicv */ {
        // Source: vinecopulib/vinecop/class.hpp:156
        const char* doc =
R"""(Evaluates the modified Bayesian information criterion for vines
(mBICV).

The mBICV is defined as

.. math:: \mathrm{mBICV} = -2\, \mathrm{loglik} + \log(n) p, - 2 *
\sum_{t=1}^(d - 1) \{q_t \log(\psi_0^t) - (d - t - q_t) \log(1
-\psi_0^t)\},

where :math:`\mathrm{loglik}` is the log-liklihood, :math:`p` is the
(effective) number of parameters of the model, :math:`t` is the tree
level, :math:`\psi_0` is the prior probability of having a
non-independence copula in the first tree, and :math:`q_t` is the
number of non-independence copulas in tree :math:`t`; The vBIC is a
consistent model selection criterion for parametric sparse vine copula
models when :math:`d = o(\sqrt{n \log n})`.

Parameter ``u``:
    An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``select()``).

Parameter ``psi0``:
    Baseline prior probability of a non-independence copula.

Parameter ``num_threads``:
    The number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } mbicv;
      // Symbol: vinecopulib::Vinecop::nobs_
      struct /* nobs_ */ {
        // Source: vinecopulib/vinecop/class.hpp:174
        const char* doc = R"""()""";
      } nobs_;
      // Symbol: vinecopulib::Vinecop::pair_copulas_
      struct /* pair_copulas_ */ {
        // Source: vinecopulib/vinecop/class.hpp:171
        const char* doc = R"""()""";
      } pair_copulas_;
      // Symbol: vinecopulib::Vinecop::pdf
      struct /* pdf */ {
        // Source: vinecopulib/vinecop/class.hpp:120
        const char* doc =
R"""(Evaluates the copula density.

The copula density is defined as joint density divided by marginal
densities, irrespective of variable types.

Parameter ``u``:
    An :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``select()``).

Parameter ``num_threads``:
    The number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } pdf;
      // Symbol: vinecopulib::Vinecop::rosenblatt
      struct /* rosenblatt */ {
        // Source: vinecopulib/vinecop/class.hpp:133
        const char* doc =
R"""(Evaluates the Rosenblatt transform for a vine copula model.

The Rosenblatt transform converts data from this model into
independent uniform variates. Only works for continuous data.

Parameter ``u``:
    An :math:`n \times d` or :math:`n \times 2d` matrix of evaluation
    points.

Parameter ``num_threads``:
    The number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } rosenblatt;
      // Symbol: vinecopulib::Vinecop::rvine_structure_
      struct /* rvine_structure_ */ {
        // Source: vinecopulib/vinecop/class.hpp:170
        const char* doc = R"""()""";
      } rvine_structure_;
      // Symbol: vinecopulib::Vinecop::select
      struct /* select */ {
        // Source: vinecopulib/vinecop/class.hpp:65
        const char* doc =
R"""(Automatically fits and selects a vine copula model.

``select()`` behaves differently depending on its current truncation
level and the truncation level specified in the controls, respectively
called ``trunc_lvl`` and ``controls.trunc_lvl`` in what follows.
Essentially, ``controls.trunc_lvl`` defines the object's truncation
level after calling ``select()``: - If ``controls.trunc_lvl <=
trunc_lvl``, the families and parameters for all pairs in trees
smaller or equal to ``controls.trunc_lvl`` are selected, using the
current structure. - If ``controls.trunc_lvl > trunc_lvl``,
``select()`` behaves as above for all trees that are smaller or equal
to ``trunc_lvl``, and then it selects the structure for higher trees
along with the families and parameters. This includes the case where
``trunc_lvl = 0``, namely where the structure is fully unspecified.

Selection of the structure is performed using the algorithm of
Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
*Selecting and estimating regular vine copulae and application to
financial returns.* Computational Statistics & Data Analysis, 59 (1),
52-69. The dependence measure used to select trees (default: Kendall's
tau) is corrected for ties (see the wdm library).

When at least one variable is discrete, two types of "observations"
are required: the first :math:`n \times d` block contains realizations
of :math:`F_Y(Y), F_X(X)`; the second :math:`n \times d` block
contains realizations of :math:`F_Y(Y^-), F_X(X^-), ...`. The minus
indicates a left-sided limit of the cdf. For continuous variables the
left limit and the cdf itself coincide. For, e.g., an integer-valued
variable, it holds :math:`F_Y(Y^-) = F_Y(Y - 1)`. Continuous variables
in the second block can be omitted.

If there are missing data (i.e., NaN entries), incomplete observations are 
discarded before fitting a pair-copula. This is done on a pair-by-pair basis
so that the maximal available information is used.

Parameter ``data``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    observations, where :math:`k` is the number of discrete variables.

Parameter ``controls``:
    The controls to the algorithm (see FitControlsVinecop).)""";
      } select;
      // Symbol: vinecopulib::Vinecop::select_all
      struct /* select_all */ {
        // Source: vinecopulib/vinecop/class.hpp:68
        const char* doc =
R"""(Automatically fits and selects a vine copula model.

Selection of the structure is performed using the algorithm of
Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
*Selecting and estimating regular vine copulae and application to
financial returns.* Computational Statistics & Data Analysis, 59 (1),
52-69.

When at least one variable is discrete, two types of "observations"
are required: the first :math:`n \times d` block contains realizations
of :math:`F_Y(Y), F_X(X)`; the second :math:`n \times d` block
contains realizations of :math:`F_Y(Y^-), F_X(X^-), ...`. The minus
indicates a left-sided limit of the cdf. For continuous variables the
left limit and the cdf itself coincide. For, e.g., an integer-valued
variable, it holds :math:`F_Y(Y^-) = F_Y(Y - 1)`. Continuous variables
in the second block can be omitted.

If there are missing data (i.e., NaN entries), incomplete observations are 
discarded before fitting a pair-copula. This is done on a pair-by-pair basis
so that the maximal available information is used.

Parameter ``data``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    observations, where :math:`k` is the number of discrete variables.

Parameter ``controls``:
    The controls to the algorithm (see FitControlsVinecop).)""";
      } select_all;
      // Symbol: vinecopulib::Vinecop::select_families
      struct /* select_families */ {
        // Source: vinecopulib/vinecop/class.hpp:72
        const char* doc =
R"""(Automatically selects all pair-copula families and fits all.
parameters.

When at least one variable is discrete, two types of "observations"
are required: the first :math:`n \times d` block contains realizations
of :math:`F_Y(Y), F_X(X)`; the second :math:`n \times d` block
contains realizations of :math:`F_Y(Y^-), F_X(X^-), ...`. The minus
indicates a left-sided limit of the cdf. For continuous variables the
left limit and the cdf itself coincide. For, e.g., an integer-valued
variable, it holds :math:`F_Y(Y^-) = F_Y(Y - 1)`. Continuous variables
in the second block can be omitted.

If there are missing data (i.e., NaN entries), incomplete observations are 
discarded before fitting a pair-copula. This is done on a pair-by-pair basis
so that the maximal available information is used.

Parameter ``data``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    observations, where :math:`k` is the number of discrete variables.

Parameter ``controls``:
    The controls to the algorithm (see FitControlsVinecop).)""";
      } select_families;
      // Symbol: vinecopulib::Vinecop::set_all_pair_copulas
      struct /* set_all_pair_copulas */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:663
        const char* doc =
R"""(Sets all pair-copulas.

Parameter ``pair_copulas``:
    A vector of pair-copulas that has to be consistent with the
    current structure (see ``Vinecop()``).)""";
      } set_all_pair_copulas;
      // Symbol: vinecopulib::Vinecop::set_continuous_var_types
      struct /* set_continuous_var_types */ {
        // Source: vinecopulib/vinecop/class.hpp:189
        const char* doc =
R"""(Sets all variable types to continuous. The function can be const,
because var_types_ is mutable.)""";
      } set_continuous_var_types;
      // Symbol: vinecopulib::Vinecop::set_var_types
      struct /* set_var_types */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:653
        const char* doc =
R"""(Sets variable types.

Parameter ``var_types``:
    A vector specifying the types of the variables, e.g., ``{"c",
    "d"}`` means first varible continuous, second discrete.)""";
      } set_var_types;
      // Symbol: vinecopulib::Vinecop::set_var_types_internal
      struct /* set_var_types_internal */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:693
        const char* doc =
R"""(Sets variable types.

Parameter ``var_types``:
    A vector specifying the types of the variables, e.g., ``{"c",
    "d"}`` means first varible continuous, second discrete.)""";
      } set_var_types_internal;
      // Symbol: vinecopulib::Vinecop::simulate
      struct /* simulate */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:902
        const char* doc =
R"""(Simulates from a vine copula model, see ``inverse_rosenblatt()``.

Simulated data is always a continous :math:`n \times d` matrix.

Parameter ``n``:
    Number of observations.

Parameter ``qrng``:
    Set to true for quasi-random numbers.

Parameter ``num_threads``:
    The number of threads to use for computations; if greater than 1,
    the function will generate ``n`` samples concurrently in
    ``num_threads`` batches.

Parameter ``seeds``:
    Seeds of the random number generator; if empty (default), the
    random number generator is seeded randomly.

Returns:
    An :math:`n \times d` matrix of samples from the copula model.)""";
      } simulate;
      // Symbol: vinecopulib::Vinecop::str
      struct /* str */ {
        // Source: vinecopulib/vinecop/class.hpp:166
        const char* doc =
R"""(Summarizes the model into a string (can be used for printing).)""";
      } str;
      // Symbol: vinecopulib::Vinecop::threshold_
      struct /* threshold_ */ {
        // Source: vinecopulib/vinecop/class.hpp:172
        const char* doc = R"""()""";
      } threshold_;
      // Symbol: vinecopulib::Vinecop::to_file
      struct /* to_file */ {
        // Source: vinecopulib/vinecop/class.hpp:62
        const char* doc =
R"""(Writes the copula object into a JSON file.

The output file contains 2 attributes : ``"structure"`` for the vine
structure, which itself contains attributes ``"array"`` for the
structure triangular array and ``"order"`` for the order vector, and
``"pair copulas"``. ``"pair copulas"`` contains a list of attributes
for the trees (``"tree1"``, ``"tree2"``, etc), each containing a list
of attributes for the edges (``"pc1"``, ``"pc2"``, etc). See the
corresponding method of ``Bicop`` objects for the encoding of
pair-copulas.

Parameter ``filename``:
    The name of the JSON file to write.)""";
      } to_file;
      // Symbol: vinecopulib::Vinecop::to_json
      struct /* to_json */ {
        // Source: vinecopulib/vinecop/class.hpp:61
        const char* doc =
R"""(Converts the copula into a nlohmann::json object.

The ``nlohmann::json`` object contains two nodes : ``"structure"`` for
the vine structure, which itself contains nodes ``"array"`` for the
structure triangular array and ``"order"`` for the order vector, and
``"pair copulas"``. The former two encode the R-Vine structure and the
latter is a list of child nodes for the trees (``"tree1"``,
``"tree2"``, etc), each containing a list of child nodes for the edges
(``"pc1"``, ``"pc2"``, etc). See Bicop::to_json() for the encoding of
pair-copulas.

Returns:
    the nlohmann::json object containing the copula.)""";
      } to_json;
      // Symbol: vinecopulib::Vinecop::truncate
      struct /* truncate */ {
        // Source: vinecopulib/vinecop/class.hpp:164
        const char* doc =
R"""(Truncates the vine copula model.

If the model is already truncated at a level less than ``trunc_lvl``,
the function does nothing.

Parameter ``trunc_lvl``:
    The truncation level.)""";
      } truncate;
      // Symbol: vinecopulib::Vinecop::var_types_
      struct /* var_types_ */ {
        // Source: vinecopulib/vinecop/class.hpp:175
        const char* doc = R"""()""";
      } var_types_;
    } Vinecop;
    // Symbol: vinecopulib::bicop_families
    struct /* bicop_families */ {
    } bicop_families;
    // Symbol: vinecopulib::family_bimap
    struct /* family_bimap */ {
      // Source: vinecopulib/bicop/implementation/family.ipp:12
      const char* doc = R"""()""";
    } family_bimap;
    // Symbol: vinecopulib::get_family_enum
    struct /* get_family_enum */ {
      // Source: vinecopulib/bicop/family.hpp:35
      const char* doc =
R"""(Converts a string name into a BicopFamily.

Parameter ``family``:
    The family name.)""";
    } get_family_enum;
    // Symbol: vinecopulib::get_family_name
    struct /* get_family_name */ {
      // Source: vinecopulib/bicop/family.hpp:32
      const char* doc =
R"""(Converts a BicopFamily into a string with its name.

Parameter ``family``:
    The family.)""";
    } get_family_name;
    // Symbol: vinecopulib::tools_select
    struct /* tools_select */ {
    } tools_select;
    // Symbol: vinecopulib::tools_stats
    struct /* tools_stats */ {
      // Symbol: vinecopulib::tools_stats::ace
      struct /* ace */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:189
        const char* doc =
R"""(alternating conditional expectation algorithm)""";
      } ace;
      // Symbol: vinecopulib::tools_stats::cef
      struct /* cef */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:167
        const char* doc =
R"""(helper routine for ace (In R, this would be win(x[ind], wl)[ranks]))""";
      } cef;
      // Symbol: vinecopulib::tools_stats::dependence_matrix
      struct /* dependence_matrix */ {
        // Source: vinecopulib/misc/tools_stats.hpp:114
        const char* doc = R"""()""";
      } dependence_matrix;
      // Symbol: vinecopulib::tools_stats::dnorm
      struct /* dnorm */ {
        // Source: vinecopulib/misc/tools_stats.hpp:22
        const char* doc =
R"""(Density function of the Standard normal distribution.

Parameter ``x``:
    Evaluation points.

Returns:
    An :math:`n \times d` matrix of evaluated densities.)""";
      } dnorm;
      // Symbol: vinecopulib::tools_stats::dt
      struct /* dt */ {
        // Source: vinecopulib/misc/tools_stats.hpp:62
        const char* doc =
R"""(Density function of the Student t distribution.

Parameter ``x``:
    Evaluation points.

Parameter ``nu``:
    Degrees of freedom parameter.

Returns:
    An :math:`n \times d` matrix of evaluated densities.)""";
      } dt;
      // Symbol: vinecopulib::tools_stats::ghalton
      struct /* ghalton */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:316
        const char* doc =
R"""(Simulates from the multivariate Generalized Halton Sequence.

For more information on Generalized Halton Sequence, see Faure, H.,
Lemieux, C. (2009). Generalized Halton Sequences in 2008: A
Comparative Study. ACM-TOMACS 19(4), Article 15.

Parameter ``n``:
    Number of observations.

Parameter ``d``:
    Dimension.

Parameter ``seeds``:
    Seeds to scramble the quasi-random numbers; if empty (default),
    the quasi-random number generator is seeded randomly.

Returns:
    An :math:`n \times d` matrix of quasi-random :math:`\mathrm{U}[0,
    1]` variables.)""";
      } ghalton;
      // Symbol: vinecopulib::tools_stats::pairwise_mcor
      struct /* pairwise_mcor */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:295
        const char* doc =
R"""(calculates the pairwise maximum correlation coefficient.)""";
      } pairwise_mcor;
      // Symbol: vinecopulib::tools_stats::pbvnorm
      struct /* pbvnorm */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:593
        const char* doc =
R"""(Compute bivariate normal probabilities.

A function for computing bivariate normal probabilities; developed
using Drezner, Z. and Wesolowsky, G. O. (1989), On the Computation of
the Bivariate Normal Integral, J. Stat. Comput. Simul.. 35 pp.
101-107. with extensive modications for double precisions by Alan Genz
and Yihong Ge. Translated from the Fortran routines of Alan Genz
(www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f).

Parameter ``z``:
    An :math:`n \times 2` matrix of evaluation points.

Parameter ``rho``:
    Correlation.

Returns:
    An :math:`n \times 1` vector of probabilities.)""";
      } pbvnorm;
      // Symbol: vinecopulib::tools_stats::pbvt
      struct /* pbvt */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:468
        const char* doc =
R"""(Computes bivariate t probabilities.

Based on the method described by Dunnett, C.W. and M. Sobel, (1954), A
bivariate generalization of Student's t-distribution with tables for
certain special cases, Biometrika 41, pp. 153-169. Translated from the
Fortran routines of Alan Genz
(www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f).

Parameter ``z``:
    An :math:`n \times 2` matrix of evaluation points.

Parameter ``nu``:
    Number of degrees of freedom.

Parameter ``rho``:
    Correlation.

Returns:
    An :math:`n \times 1` vector of probabilities.)""";
      } pbvt;
      // Symbol: vinecopulib::tools_stats::pnorm
      struct /* pnorm */ {
        // Source: vinecopulib/misc/tools_stats.hpp:35
        const char* doc =
R"""(Distribution function of the Standard normal distribution.

Parameter ``x``:
    Evaluation points.

Returns:
    An :math:`n \times d` matrix of evaluated probabilities.)""";
      } pnorm;
      // Symbol: vinecopulib::tools_stats::pt
      struct /* pt */ {
        // Source: vinecopulib/misc/tools_stats.hpp:76
        const char* doc =
R"""(Distribution function of the Student t distribution.

Parameter ``x``:
    Evaluation points.

Parameter ``nu``:
    Degrees of freedom parameter.

Returns:
    An :math:`n \times d` matrix of evaluated probabilities.)""";
      } pt;
      // Symbol: vinecopulib::tools_stats::qnorm
      struct /* qnorm */ {
        // Source: vinecopulib/misc/tools_stats.hpp:48
        const char* doc =
R"""(Quantile function of the Standard normal distribution.

Parameter ``x``:
    Evaluation points.

Returns:
    An :math:`n \times d` matrix of evaluated quantiles.)""";
      } qnorm;
      // Symbol: vinecopulib::tools_stats::qt
      struct /* qt */ {
        // Source: vinecopulib/misc/tools_stats.hpp:90
        const char* doc =
R"""(Quantile function of the Student t distribution.

Parameter ``x``:
    Evaluation points.

Parameter ``nu``:
    Degrees of freedom parameter.

Returns:
    An :math:`n \times d` matrix of evaluated quantiles.)""";
      } qt;
      // Symbol: vinecopulib::tools_stats::simulate_uniform
      struct /* simulate_uniform */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:33
        const char* doc =
R"""(Simulates from the multivariate uniform distribution.

Parameter ``n``:
    Number of observations.

Parameter ``d``:
    Dimension.

Parameter ``qrng``:
    If true, quasi-numbers are generated.

Parameter ``seeds``:
    Seeds of the random number generator; if empty (default), the
    random number generator is seeded randomly.

If ``qrng = TRUE``, generalized Halton sequences (see ``ghalton()``)
are used for :math:`d \leq 300` and Sobol sequences otherwise (see
``sobol()``).

Returns:
    An :math:`n \times d` matrix of independent :math:`\mathrm{U}[0,
    1]` random variables.)""";
      } simulate_uniform;
      // Symbol: vinecopulib::tools_stats::sobol
      struct /* sobol */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:377
        const char* doc =
R"""(Simulates from the multivariate Sobol sequence.

For more information on the Sobol sequence, see S. Joe and F. Y. Kuo
(2008), constructing Sobol sequences with better two-dimensional
projections, SIAM J. Sci. Comput. 30, 2635–2654.

Parameter ``n``:
    Number of observations.

Parameter ``d``:
    Dimension.

Parameter ``seeds``:
    Seeds to scramble the quasi-random numbers; if empty (default),
    the quasi-random number generator is seeded randomly.

Returns:
    An :math:`n \times d` matrix of quasi-random :math:`\mathrm{U}[0,
    1]` variables.)""";
      } sobol;
      // Symbol: vinecopulib::tools_stats::to_pseudo_obs
      struct /* to_pseudo_obs */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:77
        const char* doc =
R"""(Applies the empirical probability integral transform to a data matrix.

Gives pseudo-observations from the copula by applying the empirical
distribution function (scaled by :math:`n + 1`) to each margin/column.

Parameter ``x``:
    A matrix of real numbers.

Parameter ``ties_method``:
    Indicates how to treat ties; same as in R, see
    https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.

Returns:
    Pseudo-observations of the copula, i.e. :math:`F_X(x)`
    (column-wise).)""";
      } to_pseudo_obs;
      // Symbol: vinecopulib::tools_stats::to_pseudo_obs_1d
      struct /* to_pseudo_obs_1d */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:97
        const char* doc =
R"""(Applies the empirical probability integral transform to a data vector.

Gives pseudo-observations from the copula by applying the empirical
distribution function (scaled by :math:`n + 1`) to each margin/column.

Parameter ``x``:
    A vector of real numbers.

Parameter ``ties_method``:
    Indicates how to treat ties; same as in R, see
    https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.

Returns:
    Pseudo-observations of the copula, i.e. :math:`F_X(x)`.)""";
      } to_pseudo_obs_1d;
      // Symbol: vinecopulib::tools_stats::win
      struct /* win */ {
        // Source: vinecopulib/misc/implementation/tools_stats.ipp:143
        const char* doc = R"""(window smoother)""";
      } win;
    } tools_stats;
  } vinecopulib;
} pyvinecopulib_doc;

#if defined(__GNUG__)
#pragma GCC diagnostic pop
#endif
