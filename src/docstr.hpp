#pragma once
// GENERATED FILE DO NOT EDIT
// This file contains docstrings for the Python bindings that were
// automatically extracted by mkdoc.py.
#if defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif
// #include "vinecopulib/bicop/class.hpp"
// #include "vinecopulib/bicop/implementation/class.ipp"
// #include "vinecopulib/vinecop/class.hpp"
// #include "vinecopulib/vinecop/implementation/class.ipp"
// #include "vinecopulib/vinecop/implementation/rvine_structure.ipp"
// #include "vinecopulib/vinecop/rvine_structure.hpp"

// Symbol: mkdoc_doc
constexpr struct /* mkdoc_doc */ {
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
        // Source: vinecopulib/bicop/class.hpp:36
        const char* doc_1args_filename =
R"""(creates from a JSON file

Parameter ``filename``:
    the name of the JSON file to read (see to_ptree() for the
    structure of the file).)""";
        // Source: vinecopulib/bicop/class.hpp:38
        const char* doc_1args_input =
R"""(creates from a boost::property_tree::ptree object

Parameter ``input``:
    the boost::property_tree::ptree object to convert from (see
    to_ptree() for the structure of the input).)""";
        // Source: vinecopulib/bicop/implementation/class.ipp:26
        const char* doc_4args_family_rotation_parameters_var_types =
R"""(creates a specific bivariate copula model.

Parameter ``family``:
    the copula family.

Parameter ``rotation``:
    the rotation of the copula; one of 0, 90, 180, or 270 (for
    Independence, Gaussian, Student, Frank, and nonparametric
    families, only 0 is allowed).

Parameter ``parameters``:
    the copula parameters.

Parameter ``var_types``:
    a vector of size two specifying the types of the variables, e.g.,
    ``{"c", "d"}`` means first varible continuous, second discrete.)""";
        // Source: vinecopulib/bicop/implementation/class.ipp:49
        const char* doc_3args_data_controls_var_types =
R"""(create a copula model from the data, equivalent to ``Bicop cop;
cop.select(data, controls)``.

Parameter ``data``:
    see select().

Parameter ``controls``:
    see select().

Parameter ``var_types``:
    a vector of size two specifying the types of the variables, e.g.,
    ``{"c", "d"}`` means first variable continuous, second discrete.)""";
      } ctor;
      // Symbol: vinecopulib::Bicop::aic
      struct /* aic */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:333
        const char* doc =
R"""(calculates the Akaike information criterion (AIC).

The AIC is defined as

.. math:: \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p,

where :math:`\mathrm{loglik}` is the log-liklihood and :math:`p` is
the (effective) number of parameters of the model, see loglik() and
get_npars(). The AIC is a consistent model selection criterion for
nonparametric models.

Parameter ``u``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.)""";
      } aic;
      // Symbol: vinecopulib::Bicop::as_continuous
      struct /* as_continuous */ {
        // Source: vinecopulib/bicop/class.hpp:120
        const char* doc = R"""()""";
      } as_continuous;
      // Symbol: vinecopulib::Bicop::bic
      struct /* bic */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:350
        const char* doc =
R"""(calculates the Bayesian information criterion (BIC).

The BIC is defined as

.. math:: \mathrm{BIC} = -2\, \mathrm{loglik} + \ln(n) p,

where :math:`\mathrm{loglik}` is the log-liklihood and :math:`p` is
the (effective) number of parameters of the model, see loglik() and
get_npars(). The BIC is a consistent model selection criterion for
parametric models.

Parameter ``u``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.)""";
      } bic;
      // Symbol: vinecopulib::Bicop::cdf
      struct /* cdf */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:138
        const char* doc =
R"""(evaluates the copula distribution.

Parameter ``u``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.

Returns:
    The copula distribution evaluated at ``u``.)""";
      } cdf;
      // Symbol: vinecopulib::Bicop::check_data
      struct /* check_data */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:529
        const char* doc = R"""()""";
      } check_data;
      // Symbol: vinecopulib::Bicop::check_data_dim
      struct /* check_data_dim */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:536
        const char* doc = R"""()""";
      } check_data_dim;
      // Symbol: vinecopulib::Bicop::check_fitted
      struct /* check_fitted */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:908
        const char* doc =
R"""(checks whether the Bicop object was fitted to data.)""";
      } check_fitted;
      // Symbol: vinecopulib::Bicop::check_rotation
      struct /* check_rotation */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:881
        const char* doc =
R"""(checks whether the supplied rotation is valid (only 0, 90, 180, 270
allowd).)""";
      } check_rotation;
      // Symbol: vinecopulib::Bicop::check_var_types
      struct /* check_var_types */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:918
        const char* doc =
R"""(checks whether var_types have the correct length and are either "c" or
"d".)""";
      } check_var_types;
      // Symbol: vinecopulib::Bicop::check_weights_size
      struct /* check_weights_size */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:898
        const char* doc =
R"""(checks whether weights and data have matching sizes.)""";
      } check_weights_size;
      // Symbol: vinecopulib::Bicop::compute_mbic_penalty
      struct /* compute_mbic_penalty */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:500
        const char* doc = R"""()""";
      } compute_mbic_penalty;
      // Symbol: vinecopulib::Bicop::fit
      struct /* fit */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:690
        const char* doc =
R"""(fits a bivariate copula (with fixed family) to data.

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

Parameter ``data``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.

Parameter ``controls``:
    the controls (see FitControlsBicop).)""";
      } fit;
      // Symbol: vinecopulib::Bicop::flip
      struct /* flip */ {
        // Source: vinecopulib/bicop/class.hpp:114
        const char* doc =
R"""(adjust's the copula model to a change in the variable order.)""";
      } flip;
      // Symbol: vinecopulib::Bicop::flip_var_types
      struct /* flip_var_types */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:558
        const char* doc = R"""()""";
      } flip_var_types;
      // Symbol: vinecopulib::Bicop::format_data
      struct /* format_data */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:815
        const char* doc =
R"""(adds an additional column if there's only one discrete variable;
removes superfluous columns for continuous variables. (continuous
models only require two columns, discrete models always four))""";
      } format_data;
      // Symbol: vinecopulib::Bicop::get_aic
      struct /* get_aic */ {
        // Source: vinecopulib/bicop/class.hpp:60
        const char* doc = R"""(get the aic (only for fitted objects))""";
      } get_aic;
      // Symbol: vinecopulib::Bicop::get_bic
      struct /* get_bic */ {
        // Source: vinecopulib/bicop/class.hpp:61
        const char* doc = R"""(get the bic (only for fitted objects))""";
      } get_bic;
      // Symbol: vinecopulib::Bicop::get_bicop
      struct /* get_bicop */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:651
        const char* doc = R"""()""";
      } get_bicop;
      // Symbol: vinecopulib::Bicop::get_family
      struct /* get_family */ {
        // Source: vinecopulib/bicop/class.hpp:46
        const char* doc = R"""(get the copula family)""";
      } get_family;
      // Symbol: vinecopulib::Bicop::get_family_name
      struct /* get_family_name */ {
        // Source: vinecopulib/bicop/class.hpp:48
        const char* doc = R"""(get the copula family as a string)""";
      } get_family_name;
      // Symbol: vinecopulib::Bicop::get_loglik
      struct /* get_loglik */ {
        // Source: vinecopulib/bicop/class.hpp:58
        const char* doc =
R"""(get the log-likelihood (only for fitted objects))""";
      } get_loglik;
      // Symbol: vinecopulib::Bicop::get_mbic
      struct /* get_mbic */ {
        // Source: vinecopulib/bicop/class.hpp:62
        const char* doc =
R"""(get the modified bic (only for fitted objects))""";
      } get_mbic;
      // Symbol: vinecopulib::Bicop::get_n_discrete
      struct /* get_n_discrete */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:932
        const char* doc =
R"""(returns the number of discrete variables.)""";
      } get_n_discrete;
      // Symbol: vinecopulib::Bicop::get_nobs
      struct /* get_nobs */ {
        // Source: vinecopulib/bicop/class.hpp:59
        const char* doc =
R"""(get the number of observations (only for fitted objects))""";
      } get_nobs;
      // Symbol: vinecopulib::Bicop::get_npars
      struct /* get_npars */ {
        // Source: vinecopulib/bicop/class.hpp:56
        const char* doc =
R"""(returns the actual number of parameters for parameteric families.

For nonparametric families, there is a conceptually similar definition
in the sense that it can be used in the calculation of fit statistics.)""";
      } get_npars;
      // Symbol: vinecopulib::Bicop::get_parameters
      struct /* get_parameters */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:453
        const char* doc = R"""(get the parameters)""";
      } get_parameters;
      // Symbol: vinecopulib::Bicop::get_parameters_lower_bounds
      struct /* get_parameters_lower_bounds */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:636
        const char* doc =
R"""(extract lower bounds for copula parameters.)""";
      } get_parameters_lower_bounds;
      // Symbol: vinecopulib::Bicop::get_parameters_upper_bounds
      struct /* get_parameters_upper_bounds */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:643
        const char* doc =
R"""(extract upper bounds for copula parameters.)""";
      } get_parameters_upper_bounds;
      // Symbol: vinecopulib::Bicop::get_rotation
      struct /* get_rotation */ {
        // Source: vinecopulib/bicop/class.hpp:50
        const char* doc = R"""(get the rotation)""";
      } get_rotation;
      // Symbol: vinecopulib::Bicop::get_tau
      struct /* get_tau */ {
        // Source: vinecopulib/bicop/class.hpp:54
        const char* doc = R"""(get the Kendall's tau)""";
      } get_tau;
      // Symbol: vinecopulib::Bicop::get_var_types
      struct /* get_var_types */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:588
        const char* doc = R"""(extracts variable types.)""";
      } get_var_types;
      // Symbol: vinecopulib::Bicop::hfunc1
      struct /* hfunc1 */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:164
        const char* doc =
R"""(calculates the first h-function.

The first h-function is :math:`h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 =
u_1)`.

Parameter ``u``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.)""";
      } hfunc1;
      // Symbol: vinecopulib::Bicop::hfunc2
      struct /* hfunc2 */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:196
        const char* doc =
R"""(calculates the second h-function.

The second h-function is :math:`h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 =
u_2)`.

Parameter ``u``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.)""";
      } hfunc2;
      // Symbol: vinecopulib::Bicop::hinv1
      struct /* hinv1 */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:226
        const char* doc =
R"""(calculates the inverse of :math:`h_1` (see hfunc1()) w.r.t. the second
argument.

Parameter ``u``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.)""";
      } hinv1;
      // Symbol: vinecopulib::Bicop::hinv2
      struct /* hinv2 */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:256
        const char* doc =
R"""(calculates the inverse of :math:`h_2` (see hfunc2()) w.r.t. the first
argument.

Parameter ``u``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.)""";
      } hinv2;
      // Symbol: vinecopulib::Bicop::loglik
      struct /* loglik */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:311
        const char* doc =
R"""(calculates the log-likelihood.

The log-likelihood is defined as

.. math:: \mathrm{loglik} = \sum_{i = 1}^n \ln c(U_{1, i}, U_{2, i}),

where :math:`c` is the copula density pdf().

Parameter ``u``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.)""";
      } loglik;
      // Symbol: vinecopulib::Bicop::mbic
      struct /* mbic */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:376
        const char* doc =
R"""(calculates the modified Bayesian information criterion (mBIC).

The mBIC is defined as

.. math:: \mathrm{BIC} = -2\, \mathrm{loglik} + \nu \ln(n) - 2 (I
log(\psi_0) + (1 - I) log(1 - \psi_0)

where :math:`\mathrm{loglik}` is the log-liklihood and :math:`\nu` is
the (effective) number of parameters of the model, :math:`\psi_0` is
the prior probability of having a non-independence copula and
:math:`I` is an indicator for the family being non-independence; see
loglik() and get_npars().

Parameter ``u``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.

Parameter ``psi0``:
    prior probability of a non-independence copula.)""";
      } mbic;
      // Symbol: vinecopulib::Bicop::parameters_to_tau
      struct /* parameters_to_tau */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:417
        const char* doc =
R"""(converts the parameters to the Kendall's :math:`tau` for the current
family.

Parameter ``parameters``:
    the parameters (must be a valid parametrization of the current
    family).)""";
      } parameters_to_tau;
      // Symbol: vinecopulib::Bicop::pdf
      struct /* pdf */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:126
        const char* doc =
R"""(evaluates the copula density.

Parameter ``u``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.

Returns:
    The copula density evaluated at ``u``.)""";
      } pdf;
      // Symbol: vinecopulib::Bicop::prep_for_abstract
      struct /* prep_for_abstract */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:871
        const char* doc =
R"""(prepares data for use with the ``AbstractBicop`` class: - add an
additional column if there's only one discrete variable. - trim the
data to the interval [1e-10, 1 - 1e-10] for numerical stability. -
rotate the data appropriately (``AbstractBicop`` is always
0deg-rotation).)""";
      } prep_for_abstract;
      // Symbol: vinecopulib::Bicop::rotate_data
      struct /* rotate_data */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:835
        const char* doc =
R"""(rotates the data corresponding to the models rotation.

Parameter ``u``:
    an ``n x 2`` matrix.)""";
      } rotate_data;
      // Symbol: vinecopulib::Bicop::select
      struct /* select */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:733
        const char* doc =
R"""(selects the best fitting model.

The function calls fit() for all families in ``family_set`` and
selecting the best fitting model by either BIC or AIC, see bic() and
aic().

When at least one variable is discrete, two types of "observations"
are required: the first :math:`n \times 2` block contains realizations
of :math:`F_{X_1}(X_1), F_{X_2}(X_2)`. Let :math:`k` denote the number
of discrete variables (either one or two). Then the second :math:`n
\times k` block contains realizations of :math:`F_{X_k}(X_k^-)`. The
minus indicates a left-sided limit of the cdf. For continuous
variables the left limit and the cdf itself coincide. For, e.g., an
integer-valued variable, it holds :math:`F_{X_k}(X_k^-) = F_{X_k}(X_k
- 1)`.

Parameter ``data``:
    an :math:`n \times (2 + k)` matrix of observations contained in
    :math:`(0, 1)^2`, where :math:`k` is the number of discrete
    variables.

Parameter ``controls``:
    the controls (see FitControlsBicop).)""";
      } select;
      // Symbol: vinecopulib::Bicop::set_parameters
      struct /* set_parameters */ {
        // Source: vinecopulib/bicop/class.hpp:66
        const char* doc = R"""()""";
      } set_parameters;
      // Symbol: vinecopulib::Bicop::set_rotation
      struct /* set_rotation */ {
        // Source: vinecopulib/bicop/class.hpp:64
        const char* doc = R"""(set the rotation)""";
      } set_rotation;
      // Symbol: vinecopulib::Bicop::set_var_types
      struct /* set_var_types */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:574
        const char* doc =
R"""(sets variable types.

Parameter ``var_types``:
    a vector of size two specifying the types of the variables, e.g.,
    ``{"c", "d"}`` means first variable continuous, second discrete.)""";
      } set_var_types;
      // Symbol: vinecopulib::Bicop::simulate
      struct /* simulate */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:291
        const char* doc =
R"""(simulates from a bivariate copula.

Parameter ``n``:
    number of observations.

Parameter ``qrng``:
    set to true for quasi-random numbers.

Parameter ``seeds``:
    seeds of the (quasi-)random number generator; if empty (default),
    the (quasi-)random number generator is seeded randomly.

Returns:
    An :math:`n \times 2` matrix of samples from the copula model.)""";
      } simulate;
      // Symbol: vinecopulib::Bicop::str
      struct /* str */ {
        // Source: vinecopulib/bicop/class.hpp:108
        const char* doc =
R"""(summarizes the model into a string (can be used for printing).)""";
      } str;
      // Symbol: vinecopulib::Bicop::tau_to_parameters
      struct /* tau_to_parameters */ {
        // Source: vinecopulib/bicop/implementation/class.ipp:406
        const char* doc =
R"""(converts a Kendall's :math:`\tau` to the copula parameters of the
current family

(only works for one-parameter families)

Parameter ``tau``:
    a value in :math:`(-1, 1)`.)""";
      } tau_to_parameters;
      // Symbol: vinecopulib::Bicop::to_json
      struct /* to_json */ {
        // Source: vinecopulib/bicop/class.hpp:43
        const char* doc =
R"""(Write the copula object into a JSON file

See to_ptree() for the structure of the file.

Parameter ``filename``:
    the name of the file to write.)""";
      } to_json;
      // Symbol: vinecopulib::Bicop::to_ptree
      struct /* to_ptree */ {
        // Source: vinecopulib/bicop/class.hpp:41
        const char* doc =
R"""(Convert the copula into a boost::property_tree::ptree object

The boost::property_tree::ptree is contains of three values named
``"family"``, `"rotation"`, ``"parameters"``, respectively a string
for the family name, an integer for the rotation, and an
Eigen::MatrixXd for the parameters.

Returns:
    the boost::property_tree::ptree object containing the copula.)""";
      } to_ptree;
    } Bicop;
    // Symbol: vinecopulib::BicopPtr
    struct /* BicopPtr */ {
      // Source: vinecopulib/bicop/class.hpp:16
      const char* doc =
R"""(A shared pointer to an object of class AbstracBicop.)""";
    } BicopPtr;
    // Symbol: vinecopulib::CVineStructure
    struct /* CVineStructure */ {
      // Source: vinecopulib/vinecop/rvine_structure.hpp:171
      const char* doc =
R"""(C-vine structures

C-vines are a special class of R-vines where each tree is a star. A
C-vine structure is determined entirely by the order of variables. For
example, if the order is ``{1, 2, 3, 4}``, the first tree in the vine
connects variable 4 with all others, the second tree connects variable
3 with all others, etc.)""";
      // Symbol: vinecopulib::CVineStructure::CVineStructure
      struct /* ctor */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:735
        const char* doc_1args =
R"""(Parameter ``order``:
    the order of variables in the C-vine (diagonal entries in the
    R-vine array); must be a permutation of 1, ..., d.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:745
        const char* doc_2args =
R"""(Parameter ``order``:
    the order of variables in the C-vine (diagonal entries in the
    R-vine array); must be a permutation of 1, ..., d.

Parameter ``trunc_lvl``:
    the truncation level.)""";
      } ctor;
    } CVineStructure;
    // Symbol: vinecopulib::DVineStructure
    struct /* DVineStructure */ {
      // Source: vinecopulib/vinecop/rvine_structure.hpp:157
      const char* doc =
R"""(D-vine structures

D-vines are a special class of R-vines where each tree is a path. A
D-vine structure is determined entirely by the order of variables. For
example, if the order is ``{1, 2, 3, 4}``, the first tree in the vine
is 1-2-3-4 and all further trees are unique due to the proximity
condition.)""";
      // Symbol: vinecopulib::DVineStructure::DVineStructure
      struct /* ctor */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:715
        const char* doc_1args =
R"""(Parameter ``order``:
    the order of variables in the D-vine (diagonal entries in the
    R-vine array); must be a permutation of 1, ..., d.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:725
        const char* doc_2args =
R"""(Parameter ``order``:
    the order of variables in the D-vine (diagonal entries in the
    R-vine array); must be a permutation of 1, ..., d.

Parameter ``trunc_lvl``:
    the truncation level.)""";
      } ctor;
    } DVineStructure;
    // Symbol: vinecopulib::RVineStructure
    struct /* RVineStructure */ {
      // Source: vinecopulib/vinecop/rvine_structure.hpp:67
      const char* doc =
R"""(R-vine structures

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
R"""(instantiates an RVineStructure object from a matrix representing an
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
    a matrix representing a valid R-vine array.

Parameter ``check``:
    whether ``mat`` shall be checked for validity.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:64
        const char* doc_2args_d_trunc_lvl =
R"""(instantiates an RVineStructure object to a D-vine for a given
dimension

Parameter ``d``:
    the dimension.

Parameter ``trunc_lvl``:
    the truncation level. By default, it is dim - 1.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:74
        const char* doc_3args_order_trunc_lvl_check =
R"""(instantiates an RVineStructure object to a D-vine with given ordering
of variables.

Parameter ``order``:
    the order of variables in the D-vine (diagonal entries in the
    R-vine array); must be a permutation of 1, ..., d.

Parameter ``trunc_lvl``:
    the truncation level. By default, it is d - 1.

Parameter ``check``:
    whether `order shall be checked for validity.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:100
        const char* doc_4args_order_struct_array_natural_order_check =
R"""(instantiates an RVineStructure object from the variable order
(diagonal elements of the R-vine array) and a triangular structure
array (all elements above the diagonal).

Parameter ``order``:
    the order of variables (diagonal entries in the R-vine array);
    must be a permutation of 1, ..., d.

Parameter ``struct_array``:
    the structure array (all elements above the diagonal in the R-vine
    array). For truncated vines, all rows below the truncation level
    are omitted.

Parameter ``natural_order``:
    whether ``struct_array`` is already in natural order.

Parameter ``check``:
    whether ``order`` and ``struct_array`` shall be checked for
    validity.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:145
        const char* doc_2args_input_check =
R"""(creates from a boost::property_tree::ptree object

Parameter ``input``:
    the boost::property_tree::ptree object to convert from (see
    to_ptree() for the structure of the input).

Parameter ``check``:
    whether to check if the input represents a valid R-vine structure.)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:160
        const char* doc_2args_filename_check =
R"""(creates from a JSON file.

Parameter ``filename``:
    the name of the JSON file to read (see to_ptree() for the
    structure of the file).

Parameter ``check``:
    whether to check if the input represents a valid R-vine matrix.)""";
      } ctor;
      // Symbol: vinecopulib::RVineStructure::check_antidiagonal
      struct /* check_antidiagonal */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:660
        const char* doc = R"""()""";
      } check_antidiagonal;
      // Symbol: vinecopulib::RVineStructure::check_columns
      struct /* check_columns */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:633
        const char* doc = R"""()""";
      } check_columns;
      // Symbol: vinecopulib::RVineStructure::check_if_quadratic
      struct /* check_if_quadratic */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:593
        const char* doc = R"""()""";
      } check_if_quadratic;
      // Symbol: vinecopulib::RVineStructure::check_lower_tri
      struct /* check_lower_tri */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:603
        const char* doc = R"""()""";
      } check_lower_tri;
      // Symbol: vinecopulib::RVineStructure::check_proximity_condition
      struct /* check_proximity_condition */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:671
        const char* doc = R"""()""";
      } check_proximity_condition;
      // Symbol: vinecopulib::RVineStructure::check_upper_tri
      struct /* check_upper_tri */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:617
        const char* doc = R"""()""";
      } check_upper_tri;
      // Symbol: vinecopulib::RVineStructure::compute_min_array
      struct /* compute_min_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:549
        const char* doc = R"""()""";
      } compute_min_array;
      // Symbol: vinecopulib::RVineStructure::compute_needed_hfunc1
      struct /* compute_needed_hfunc1 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:562
        const char* doc = R"""()""";
      } compute_needed_hfunc1;
      // Symbol: vinecopulib::RVineStructure::compute_needed_hfunc2
      struct /* compute_needed_hfunc2 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:577
        const char* doc = R"""()""";
      } compute_needed_hfunc2;
      // Symbol: vinecopulib::RVineStructure::d_
      struct /* d_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:139
        const char* doc = R"""()""";
      } d_;
      // Symbol: vinecopulib::RVineStructure::find_trunc_lvl
      struct /* find_trunc_lvl */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:446
        const char* doc =
R"""(find the truncation level in an R-vine array.

The truncation level is determined by the first row (starting from the
bottom) that contains only zeros above the diagonal.

Parameter ``mat``:
    an array representing the R-vine array.)""";
      } find_trunc_lvl;
      // Symbol: vinecopulib::RVineStructure::get_dim
      struct /* get_dim */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:196
        const char* doc = R"""(extract the dimension of the vine.)""";
      } get_dim;
      // Symbol: vinecopulib::RVineStructure::get_matrix
      struct /* get_matrix */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:423
        const char* doc =
R"""(extract the R-vine matrix representation.)""";
      } get_matrix;
      // Symbol: vinecopulib::RVineStructure::get_min_array
      struct /* get_min_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:244
        const char* doc =
R"""(extracts the minimum array.

The minimum array is derived from an R-vine array by iteratively
computing the (elementwise) minimum of two subsequent rows (starting
from the top). It is used in estimation and evaluation algorithms to
find the two edges in the previous tree that are joined by the current
edge.)""";
      } get_min_array;
      // Symbol: vinecopulib::RVineStructure::get_needed_hfunc1
      struct /* get_needed_hfunc1 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:255
        const char* doc =
R"""(extracts an array indicating which of the first h-functions are
needed.

(it is usually not necessary to compute both h-functions for each
pair-copula).)""";
      } get_needed_hfunc1;
      // Symbol: vinecopulib::RVineStructure::get_needed_hfunc2
      struct /* get_needed_hfunc2 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:266
        const char* doc =
R"""(extracts an array indicating which of the second h-functions are
needed.

(it is usually not necessary to compute both h-functions for each
pair-copula).)""";
      } get_needed_hfunc2;
      // Symbol: vinecopulib::RVineStructure::get_order
      struct /* get_order */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:211
        const char* doc_0args =
R"""(extract the order of variables in the vine (diagonal entries in the
R-vine array).)""";
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:474
        const char* doc_1args =
R"""(find the order of an R-vine array.

The order is contained in the counter-diagonal of the R-vine array.

Parameter ``mat``:
    a matrix representing the R-vine array.)""";
      } get_order;
      // Symbol: vinecopulib::RVineStructure::get_struct_array
      struct /* get_struct_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:220
        const char* doc =
R"""(extract structure array (all elements above the diagonal in the R-vine
array).

Parameter ``natural_order``:
    whether indices correspond to natural order.)""";
      } get_struct_array;
      // Symbol: vinecopulib::RVineStructure::get_trunc_lvl
      struct /* get_trunc_lvl */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:203
        const char* doc =
R"""(extract the truncation level of the vine.)""";
      } get_trunc_lvl;
      // Symbol: vinecopulib::RVineStructure::make_cvine_struct_array
      struct /* make_cvine_struct_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:536
        const char* doc =
R"""(creates a structure array corresponding to a D-vine (in natural
order).)""";
      } make_cvine_struct_array;
      // Symbol: vinecopulib::RVineStructure::make_dvine_struct_array
      struct /* make_dvine_struct_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:522
        const char* doc =
R"""(creates a structure array corresponding to a D-vine (in natural
order).)""";
      } make_dvine_struct_array;
      // Symbol: vinecopulib::RVineStructure::min_array
      struct /* min_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:289
        const char* doc =
R"""(access elements of the minimum array.

Parameter ``tree``:
    tree index.

Parameter ``edge``:
    edge index.)""";
      } min_array;
      // Symbol: vinecopulib::RVineStructure::min_array_
      struct /* min_array_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:142
        const char* doc = R"""()""";
      } min_array_;
      // Symbol: vinecopulib::RVineStructure::needed_hfunc1
      struct /* needed_hfunc1 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:298
        const char* doc =
R"""(access elements of the needed_hfunc1 array.

Parameter ``tree``:
    tree index.

Parameter ``edge``:
    edge index.)""";
      } needed_hfunc1;
      // Symbol: vinecopulib::RVineStructure::needed_hfunc1_
      struct /* needed_hfunc1_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:144
        const char* doc = R"""()""";
      } needed_hfunc1_;
      // Symbol: vinecopulib::RVineStructure::needed_hfunc2
      struct /* needed_hfunc2 */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:305
        const char* doc =
R"""(access elements of the needed_hfunc2 array.)""";
      } needed_hfunc2;
      // Symbol: vinecopulib::RVineStructure::needed_hfunc2_
      struct /* needed_hfunc2_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:145
        const char* doc = R"""()""";
      } needed_hfunc2_;
      // Symbol: vinecopulib::RVineStructure::order_
      struct /* order_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:138
        const char* doc = R"""()""";
      } order_;
      // Symbol: vinecopulib::RVineStructure::simulate
      struct /* simulate */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:357
        const char* doc =
R"""(randomly sample a regular vine structure.

Parameter ``d``:
    the dimension.

Parameter ``natural_order``:
    should the sampled structure be in natural order?

Parameter ``seeds``:
    seeds of the random number generator; if empty (default), the
    random number generator is seeded randomly.

Note:
    Implementation of Algorithm 13 in Harry Joe's 2014 book (p. 288),
    but there's a typo: the end of line 6 in the book should be
    'column j' instead of 'column k'.)""";
      } simulate;
      // Symbol: vinecopulib::RVineStructure::str
      struct /* str */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:330
        const char* doc =
R"""(converts the structure to a string representation (most useful for
printing).)""";
      } str;
      // Symbol: vinecopulib::RVineStructure::struct_array
      struct /* struct_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:276
        const char* doc =
R"""(access elements of the structure array.

Parameter ``tree``:
    tree index.

Parameter ``edge``:
    edge index.

Parameter ``natural_order``:
    whether indices correspond to natural order.)""";
      } struct_array;
      // Symbol: vinecopulib::RVineStructure::struct_array_
      struct /* struct_array_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:141
        const char* doc = R"""()""";
      } struct_array_;
      // Symbol: vinecopulib::RVineStructure::to_json
      struct /* to_json */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:189
        const char* doc =
R"""(write the structure into a JSON file.

See to_ptree() for the structure of the file.

Parameter ``filename``:
    the name of the file to write.)""";
      } to_json;
      // Symbol: vinecopulib::RVineStructure::to_natural_order
      struct /* to_natural_order */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:504
        const char* doc =
R"""(converts ``struct_array_`` to natural order.)""";
      } to_natural_order;
      // Symbol: vinecopulib::RVineStructure::to_ptree
      struct /* to_ptree */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:172
        const char* doc =
R"""(converts the structure into a boost::property_tree::ptree object.

The ``ptree`` object contains two nodes: ``"array"`` for the structure
triangular array and ``"order"`` for the order vector.

Returns:
    the boost::property_tree::ptree object containing the structure.)""";
      } to_ptree;
      // Symbol: vinecopulib::RVineStructure::to_rvine_array
      struct /* to_rvine_array */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:488
        const char* doc =
R"""(extracts the structure array (entries above the diagonal in R-vine
array).

Parameter ``mat``:
    a array representing the R-vine array.)""";
      } to_rvine_array;
      // Symbol: vinecopulib::RVineStructure::trunc_lvl_
      struct /* trunc_lvl_ */ {
        // Source: vinecopulib/vinecop/rvine_structure.hpp:140
        const char* doc = R"""()""";
      } trunc_lvl_;
      // Symbol: vinecopulib::RVineStructure::truncate
      struct /* truncate */ {
        // Source: vinecopulib/vinecop/implementation/rvine_structure.ipp:316
        const char* doc =
R"""(truncates the R-vine structure.

Parameter ``trunc_lvl``:
    the truncation level.

If the structure is already truncated at a level less than
``trunc_lvl``, the function does nothing.)""";
      } truncate;
    } RVineStructure;
    // Symbol: vinecopulib::Vinecop
    struct /* Vinecop */ {
      // Source: vinecopulib/vinecop/class.hpp:25
      const char* doc =
R"""(A class for vine copula models

A vine copula model is characterized by the structure matrix (see
TriangularArray) and the pair-copulas.)""";
      // Symbol: vinecopulib::Vinecop::Vinecop
      struct /* ctor */ {
        // Source: vinecopulib/vinecop/class.hpp:29
        const char* doc_0args = R"""()""";
        // Source: vinecopulib/vinecop/class.hpp:31
        const char* doc_1args_d =
R"""(creates a D-vine on ``d`` variables with all pair-copulas set to
independence.

Parameter ``d``:
    the dimension (= number of variables) of the model.)""";
        // Source: vinecopulib/vinecop/class.hpp:55
        const char* doc_2args_filename_check =
R"""(creates from a JSON file.

Parameter ``filename``:
    the name of the JSON file to read (see to_ptree() for the
    structure of the file).

Parameter ``check``:
    whether to check if the ``"structure"`` node of the input
    represents a valid R-vine structure.)""";
        // Source: vinecopulib/vinecop/class.hpp:56
        const char* doc_2args_input_check =
R"""(creates from a boost::property_tree::ptree object

Parameter ``input``:
    the boost::property_tree::ptree object to convert from (see
    to_ptree() for the structure of the input).

Parameter ``check``:
    whether to check if the ``"structure"`` node represents a valid
    R-vine structure.)""";
        // Source: vinecopulib/vinecop/implementation/class.ipp:32
        const char* doc_3args_structure_pair_copulas_var_types =
R"""(creates an arbitrary vine copula model.

Parameter ``structure``:
    an RVineStructure object specifying the vine structure.

Parameter ``pair_copulas``:
    Bicop objects specifying the pair-copulas, see
    make_pair_copula_store().

Parameter ``var_types``:
    a vector specifying the types of the variables, e.g., ``{"c",
    "d"}`` means first variable continuous, second discrete. If empty,
    then all variables are set as continuous.)""";
        // Source: vinecopulib/vinecop/implementation/class.ipp:56
        const char* doc_3args_matrix_pair_copulas_var_types =
R"""(creates an arbitrary vine copula model.

Parameter ``matrix``:
    an R-vine matrix specifying the vine structure.

Parameter ``pair_copulas``:
    Bicop objects specifying the pair-copulas, see
    make_pair_copula_store().

Parameter ``var_types``:
    a vector specifying the types of the variables, e.g., ``{"c",
    "d"}`` means first variable continuous, second discrete. If empty,
    then all variables are set as continuous.)""";
        // Source: vinecopulib/vinecop/implementation/class.ipp:73
        const char* doc_4args_data_structure_var_types_controls =
R"""(constructs a vine copula model from data by creating a model and
calling select().

Parameter ``data``:
    an :math:`n \times d` matrix of observations.

Parameter ``structure``:
    an RVineStructure object specifying the vine structure. If empty,
    then it is selected as part of the fit.

Parameter ``var_types``:
    a vector specifying the types of the variables, e.g., ``{"c",
    "d"}`` means first variable continuous, second discrete. If empty,
    then all variables are set as continuous.

Parameter ``controls``:
    see FitControlsVinecop.)""";
        // Source: vinecopulib/vinecop/implementation/class.ipp:109
        const char* doc_4args_data_matrix_var_types_controls =
R"""(constructs a vine copula model from data by creating a model and
calling select().

Parameter ``data``:
    an :math:`n \times d` matrix of observations.

Parameter ``matrix``:
    either an empty matrix (default) or an R-vine structure matrix,
    see select(). If empty, then it is selected as part of the fit.

Parameter ``var_types``:
    a vector specifying the types of the variables, e.g., ``{"c",
    "d"}`` means first variable continuous, second discrete. If empty,
    then all variables are set as continuous.

Parameter ``controls``:
    see FitControlsVinecop.)""";
      } ctor;
      // Symbol: vinecopulib::Vinecop::aic
      struct /* aic */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:918
        const char* doc =
R"""(calculates the Akaike information criterion (AIC).

The AIC is defined as

.. math:: \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p,

where :math:`\mathrm{loglik}` is the log-liklihood and :math:`p` is
the (effective) number of parameters of the model, see loglik() and
get_npars(). The AIC is a consistent model selection criterion for
nonparametric models.

Parameter ``u``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``Vinecop::select()``).

Parameter ``num_threads``:
    the number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } aic;
      // Symbol: vinecopulib::Vinecop::bic
      struct /* bic */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:939
        const char* doc =
R"""(calculates the Bayesian information criterion (BIC).

The BIC is defined as

.. math:: \mathrm{BIC} = -2\, \mathrm{loglik} + \ln(n) p,

where :math:`\mathrm{loglik}` is the log-liklihood and :math:`p` is
the (effective) number of parameters of the model, see loglik() and
get_npars(). The BIC is a consistent model selection criterion for
nonparametric models.

Parameter ``u``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``Vinecop::select()``).

Parameter ``num_threads``:
    the number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } bic;
      // Symbol: vinecopulib::Vinecop::calculate_mbicv_penalty
      struct /* calculate_mbicv_penalty */ {
        // Source: vinecopulib/vinecop/class.hpp:179
        const char* doc = R"""(computes the penalty term for mBICV)""";
      } calculate_mbicv_penalty;
      // Symbol: vinecopulib::Vinecop::cdf
      struct /* cdf */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:823
        const char* doc =
R"""(calculates the cumulative distribution of the vine copula model.

Parameter ``u``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``Vinecop::select()``).

Parameter ``N``:
    integer for the number of quasi-random numbers to draw to evaluate
    the distribution (default: 1e4).

Parameter ``num_threads``:
    the number of threads to use for computations; if greater than 1,
    the function will generate ``n`` samples concurrently in
    ``num_threads`` batches.

Parameter ``seeds``:
    seeds to scramble the quasi-random numbers; if empty (default),
    the random number quasi-generator is seeded randomly.)""";
      } cdf;
      // Symbol: vinecopulib::Vinecop::check_data
      struct /* check_data */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:1200
        const char* doc =
R"""(checks if dimension d of the data matches the dimension of the vine.)""";
      } check_data;
      // Symbol: vinecopulib::Vinecop::check_data_dim
      struct /* check_data_dim */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:1177
        const char* doc =
R"""(checks if dimension d of the data matches the dimension of the vine.)""";
      } check_data_dim;
      // Symbol: vinecopulib::Vinecop::check_enough_data
      struct /* check_enough_data */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:1253
        const char* doc = R"""(checks if data size is large enough)""";
      } check_enough_data;
      // Symbol: vinecopulib::Vinecop::check_fitted
      struct /* check_fitted */ {
        // Source: vinecopulib/vinecop/class.hpp:184
        const char* doc = R"""()""";
      } check_fitted;
      // Symbol: vinecopulib::Vinecop::check_pair_copulas_rvine_structure
      struct /* check_pair_copulas_rvine_structure */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:1208
        const char* doc =
R"""(checks if pair copulas are compatible with the R-vine structure.)""";
      } check_pair_copulas_rvine_structure;
      // Symbol: vinecopulib::Vinecop::check_var_types
      struct /* check_var_types */ {
        // Source: vinecopulib/vinecop/class.hpp:185
        const char* doc = R"""()""";
      } check_var_types;
      // Symbol: vinecopulib::Vinecop::check_weights_size
      struct /* check_weights_size */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:1243
        const char* doc =
R"""(checks if weights are compatible with the data.)""";
      } check_weights_size;
      // Symbol: vinecopulib::Vinecop::collapse_data
      struct /* collapse_data */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:1305
        const char* doc =
R"""(removes superfluous columns for continuous data.)""";
      } collapse_data;
      // Symbol: vinecopulib::Vinecop::d_
      struct /* d_ */ {
        // Source: vinecopulib/vinecop/class.hpp:167
        const char* doc = R"""()""";
      } d_;
      // Symbol: vinecopulib::Vinecop::finalize_fit
      struct /* finalize_fit */ {
        // Source: vinecopulib/vinecop/class.hpp:180
        const char* doc = R"""()""";
      } finalize_fit;
      // Symbol: vinecopulib::Vinecop::get_aic
      struct /* get_aic */ {
        // Source: vinecopulib/vinecop/class.hpp:113
        const char* doc =
R"""(extracts the AIC.

The function throws an error if model has not been fitted to data.)""";
      } get_aic;
      // Symbol: vinecopulib::Vinecop::get_all_families
      struct /* get_all_families */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:387
        const char* doc =
R"""(extracts the families of all pair copulas.

Returns:
    a nested std::vector with entry ``[t][e]`` corresponding to edge
    ``e`` in tree ``t``.)""";
      } get_all_families;
      // Symbol: vinecopulib::Vinecop::get_all_pair_copulas
      struct /* get_all_pair_copulas */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:367
        const char* doc =
R"""(extracts all pair copulas.

Returns:
    a nested std::vector with entry ``[t][e]`` corresponding to edge
    ``e`` in tree ``t``.)""";
      } get_all_pair_copulas;
      // Symbol: vinecopulib::Vinecop::get_all_rotations
      struct /* get_all_rotations */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:415
        const char* doc =
R"""(extracts the rotations of all pair copulas.

Returns:
    a nested std::vector with entry ``[t][e]`` corresponding to edge
    ``e`` in tree ``t``.)""";
      } get_all_rotations;
      // Symbol: vinecopulib::Vinecop::get_all_taus
      struct /* get_all_taus */ {
        // Source: vinecopulib/vinecop/class.hpp:96
        const char* doc = R"""()""";
      } get_all_taus;
      // Symbol: vinecopulib::Vinecop::get_bic
      struct /* get_bic */ {
        // Source: vinecopulib/vinecop/class.hpp:114
        const char* doc =
R"""(extracts the BIC.

The function throws an error if model has not been fitted to data.)""";
      } get_bic;
      // Symbol: vinecopulib::Vinecop::get_dim
      struct /* get_dim */ {
        // Source: vinecopulib/vinecop/class.hpp:99
        const char* doc = R"""()""";
      } get_dim;
      // Symbol: vinecopulib::Vinecop::get_family
      struct /* get_family */ {
        // Source: vinecopulib/vinecop/class.hpp:77
        const char* doc =
R"""(extracts the family of a pair copula.

Parameter ``tree``:
    tree index (starting with 0).

Parameter ``edge``:
    edge index (starting with 0).)""";
      } get_family;
      // Symbol: vinecopulib::Vinecop::get_loglik
      struct /* get_loglik */ {
        // Source: vinecopulib/vinecop/class.hpp:111
        const char* doc =
R"""(extracts the log-likelihood (throws an error if model has not been
fitted to data).)""";
      } get_loglik;
      // Symbol: vinecopulib::Vinecop::get_matrix
      struct /* get_matrix */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:513
        const char* doc =
R"""(extracts the structure matrix of the vine copula model.)""";
      } get_matrix;
      // Symbol: vinecopulib::Vinecop::get_mbicv
      struct /* get_mbicv */ {
        // Source: vinecopulib/vinecop/class.hpp:115
        const char* doc =
R"""(extracts the log-likelihood.

The function throws an error if model has not been fitted to data.)""";
      } get_mbicv;
      // Symbol: vinecopulib::Vinecop::get_n_discrete
      struct /* get_n_discrete */ {
        // Source: vinecopulib/vinecop/class.hpp:188
        const char* doc =
R"""(returns the number of discrete variables.)""";
      } get_n_discrete;
      // Symbol: vinecopulib::Vinecop::get_nobs
      struct /* get_nobs */ {
        // Source: vinecopulib/vinecop/class.hpp:112
        const char* doc =
R"""(extracts the number of observations used for the fit.

The function throws an error if model has not been fitted to data.)""";
      } get_nobs;
      // Symbol: vinecopulib::Vinecop::get_npars
      struct /* get_npars */ {
        // Source: vinecopulib/vinecop/class.hpp:143
        const char* doc =
R"""(returns sum of the number of parameters for all pair copulas (see
Bicop::get_npars()).)""";
      } get_npars;
      // Symbol: vinecopulib::Vinecop::get_order
      struct /* get_order */ {
        // Source: vinecopulib/vinecop/class.hpp:101
        const char* doc = R"""()""";
      } get_order;
      // Symbol: vinecopulib::Vinecop::get_pair_copula
      struct /* get_pair_copula */ {
        // Source: vinecopulib/vinecop/class.hpp:75
        const char* doc =
R"""(extracts a pair copula.

Parameter ``tree``:
    tree index (starting with 0).

Parameter ``edge``:
    edge index (starting with 0).)""";
      } get_pair_copula;
      // Symbol: vinecopulib::Vinecop::get_parameters
      struct /* get_parameters */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:433
        const char* doc =
R"""(extracts the parameters of a pair copula.

Parameter ``tree``:
    tree index (starting with 0).

Parameter ``edge``:
    edge index (starting with 0).)""";
      } get_parameters;
      // Symbol: vinecopulib::Vinecop::get_rotation
      struct /* get_rotation */ {
        // Source: vinecopulib/vinecop/class.hpp:79
        const char* doc =
R"""(extracts the rotation of a pair copula.

Parameter ``tree``:
    tree index (starting with 0).

Parameter ``edge``:
    edge index (starting with 0).)""";
      } get_rotation;
      // Symbol: vinecopulib::Vinecop::get_rvine_structure
      struct /* get_rvine_structure */ {
        // Source: vinecopulib/vinecop/class.hpp:103
        const char* doc = R"""()""";
      } get_rvine_structure;
      // Symbol: vinecopulib::Vinecop::get_struct_array
      struct /* get_struct_array */ {
        // Source: vinecopulib/vinecop/class.hpp:107
        const char* doc =
R"""(extracts the above diagonal coefficients of the vine copula model.

Parameter ``natural_order``:
    whether indices correspond to natural order.)""";
      } get_struct_array;
      // Symbol: vinecopulib::Vinecop::get_tau
      struct /* get_tau */ {
        // Source: vinecopulib/vinecop/class.hpp:83
        const char* doc =
R"""(extracts the Kendall's :math:`tau` of a pair copula.

Parameter ``tree``:
    tree index (starting with 0).

Parameter ``edge``:
    edge index (starting with 0).)""";
      } get_tau;
      // Symbol: vinecopulib::Vinecop::get_threshold
      struct /* get_threshold */ {
        // Source: vinecopulib/vinecop/class.hpp:110
        const char* doc =
R"""(extracts the threshold (usually zero except ``select_threshold ==
TRUE`` in ``FitControlsVinecop()``).)""";
      } get_threshold;
      // Symbol: vinecopulib::Vinecop::get_trunc_lvl
      struct /* get_trunc_lvl */ {
        // Source: vinecopulib/vinecop/class.hpp:85
        const char* doc = R"""()""";
      } get_trunc_lvl;
      // Symbol: vinecopulib::Vinecop::get_var_types
      struct /* get_var_types */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:699
        const char* doc = R"""(extracts the variable types.)""";
      } get_var_types;
      // Symbol: vinecopulib::Vinecop::inverse_rosenblatt
      struct /* inverse_rosenblatt */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:1083
        const char* doc =
R"""(calculates the inverse Rosenblatt transform for a vine copula model.

The inverse Rosenblatt transform can be used for simulation: the
function applied to independent uniform variates resembles simulated
data from the vine copula model.

If the problem is too large, it is split recursively into halves
(w.r.t. n, the number of observations). "Too large" means that the
required memory will exceed 1 GB. An examplary configuration requiring
less than 1 GB is :math:`n = 1000`, :math:`d = 200`.

Only works for continous models.

Parameter ``u``:
    :math:`n \times d` matrix of evaluation points.

Parameter ``num_threads``:
    the number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } inverse_rosenblatt;
      // Symbol: vinecopulib::Vinecop::loglik
      struct /* loglik */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:893
        const char* doc =
R"""(calculates the log-likelihood.

The log-likelihood is defined as

.. math:: \mathrm{loglik} = \sum_{i = 1}^n \ln c(U_{1, i}, ..., U_{d,
i}),

where :math:`c` is the copula density pdf().

Parameter ``u``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``Vinecop::select()``).

Parameter ``num_threads``:
    the number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } loglik;
      // Symbol: vinecopulib::Vinecop::loglik_
      struct /* loglik_ */ {
        // Source: vinecopulib/vinecop/class.hpp:171
        const char* doc = R"""()""";
      } loglik_;
      // Symbol: vinecopulib::Vinecop::make_pair_copula_store
      struct /* make_pair_copula_store */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:222
        const char* doc =
R"""(initializes object for storing pair copulas.

Parameter ``d``:
    dimension of the vine copula.

Parameter ``trunc_lvl``:
    a truncation level (optional).

Returns:
    A nested vector such that ``pc_store[t][e]`` contains a Bicop.
    object for the pair copula corresponding to tree ``t`` and edge
    ``e``.)""";
      } make_pair_copula_store;
      // Symbol: vinecopulib::Vinecop::mbicv
      struct /* mbicv */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:966
        const char* doc =
R"""(calculates the modified Bayesian information criterion for vines
(mBICV).

The mBICV is defined as

.. math:: \mathrm{mBICV} = -2\, \mathrm{loglik} + \ln(n) \nu, - 2 *
\sum_{t=1}^(d - 1) \{q_t log(\psi_0^t) - (d - t - q_t) log(1
-\psi_0^t)\}

where :math:`\mathrm{loglik}` is the log-liklihood, :math:`\nu` is the
(effective) number of parameters of the model, :math:`t` is the tree
level :math:`\psi_0` is the prior probability of having a
non-independence copula in the first tree, and :math:`q_t` is the
number of non-independence copulas in tree :math:`t`; The vBIC is a
consistent model selection criterion for parametric sparse vine copula
models when :math:`d = o(\sqrt{n \ln n})`.

Parameter ``u``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``Vinecop::select()``).

Parameter ``psi0``:
    baseline prior probability of a non-independence copula.

Parameter ``num_threads``:
    the number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } mbicv;
      // Symbol: vinecopulib::Vinecop::nobs_
      struct /* nobs_ */ {
        // Source: vinecopulib/vinecop/class.hpp:172
        const char* doc = R"""()""";
      } nobs_;
      // Symbol: vinecopulib::Vinecop::pair_copulas_
      struct /* pair_copulas_ */ {
        // Source: vinecopulib/vinecop/class.hpp:169
        const char* doc = R"""()""";
      } pair_copulas_;
      // Symbol: vinecopulib::Vinecop::pdf
      struct /* pdf */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:715
        const char* doc =
R"""(calculates the density function of the vine copula model.

Parameter ``u``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    evaluation points, where :math:`k` is the number of discrete
    variables (see ``Vinecop::select()``).

Parameter ``num_threads``:
    the number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } pdf;
      // Symbol: vinecopulib::Vinecop::rosenblatt
      struct /* rosenblatt */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:1001
        const char* doc =
R"""(calculates the Rosenblatt transform for a vine copula model.

The Rosenblatt transform converts data from this model into
independent uniform variates. Only works for continuous data.

Parameter ``u``:
    :math:`n \times d` or :math:`n \times 2d` matrix of evaluation
    points.

Parameter ``num_threads``:
    the number of threads to use for computations; if greater than 1,
    the function will be applied concurrently to ``num_threads``
    batches of ``u``.)""";
      } rosenblatt;
      // Symbol: vinecopulib::Vinecop::rvine_structure_
      struct /* rvine_structure_ */ {
        // Source: vinecopulib/vinecop/class.hpp:168
        const char* doc = R"""()""";
      } rvine_structure_;
      // Symbol: vinecopulib::Vinecop::select
      struct /* select */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:262
        const char* doc =
R"""(automatically fits and selects a vine copula model.

``select()`` behaves differently depending on its current truncation
level and the truncation level specified in the controls, respectively
called ``trunc_lvl`` and ``controls.trunc_lvl`` in what follows.
Essentially, ``controls.trunc_lvl`` defines the object's truncation
level after calling ``select()``: - If ``controls.trunc_lvl <=
trunc_lvl``, the families and parameters for all pairs in trees
smaller or equal to ``controls.trunc_lvl`` are selected, using the
current structure. - If ``controls.trunc_lvl > trunc_lvl``, `select()`
behaves as above for all trees that are smaller or equal to
``trunc_lvl``, and then it selects the structure for higher trees
along with the families and parameters. This includes the case where
``trunc_lvl = 0``, namely where the structure is fully unspecified.

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

Parameter ``data``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    observations, where :math:`k` is the number of discrete variables.

Parameter ``controls``:
    the controls to the algorithm (see FitControlsVinecop).)""";
      } select;
      // Symbol: vinecopulib::Vinecop::select_all
      struct /* select_all */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:298
        const char* doc =
R"""(automatically fits and selects a vine copula model.

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

Parameter ``data``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    observations, where :math:`k` is the number of discrete variables.

Parameter ``controls``:
    the controls to the algorithm (see FitControlsVinecop).)""";
      } select_all;
      // Symbol: vinecopulib::Vinecop::select_families
      struct /* select_families */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:322
        const char* doc =
R"""(automatically selects all pair-copula families and fits all
parameters.

When at least one variable is discrete, two types of "observations"
are required: the first :math:`n \times d` block contains realizations
of :math:`F_Y(Y), F_X(X)`; the second :math:`n \times d` block
contains realizations of :math:`F_Y(Y^-), F_X(X^-), ...`. The minus
indicates a left-sided limit of the cdf. For continuous variables the
left limit and the cdf itself coincide. For, e.g., an integer-valued
variable, it holds :math:`F_Y(Y^-) = F_Y(Y - 1)`. Continuous variables
in the second block can be omitted.

Parameter ``data``:
    :math:`n \times (d + k)` or :math:`n \times 2d` matrix of
    observations, where :math:`k` is the number of discrete variables.

Parameter ``controls``:
    the controls to the algorithm (see FitControlsVinecop).)""";
      } select_families;
      // Symbol: vinecopulib::Vinecop::set_all_pair_copulas
      struct /* set_all_pair_copulas */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:633
        const char* doc =
R"""(sets all pair-copulas.

Parameter ``pair_copulas``:
    a vector of pair-copulas that has to be consistent with the
    current structure (see Vinecop()).)""";
      } set_all_pair_copulas;
      // Symbol: vinecopulib::Vinecop::set_continuous_var_types
      struct /* set_continuous_var_types */ {
        // Source: vinecopulib/vinecop/class.hpp:186
        const char* doc =
R"""(set all variable types to continuous. the function can be const,
because var_types_ is mutable.)""";
      } set_continuous_var_types;
      // Symbol: vinecopulib::Vinecop::set_var_types
      struct /* set_var_types */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:623
        const char* doc =
R"""(sets variable types.

Parameter ``var_types``:
    a vector specifying the types of the variables, e.g., ``{"c",
    "d"}`` means first varible continuous, second discrete.)""";
      } set_var_types;
      // Symbol: vinecopulib::Vinecop::set_var_types_internal
      struct /* set_var_types_internal */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:663
        const char* doc =
R"""(sets variable types.

Parameter ``var_types``:
    a vector specifying the types of the variables, e.g., ``{"c",
    "d"}`` means first varible continuous, second discrete.)""";
      } set_var_types_internal;
      // Symbol: vinecopulib::Vinecop::simulate
      struct /* simulate */ {
        // Source: vinecopulib/vinecop/implementation/class.ipp:866
        const char* doc =
R"""(simulates from a vine copula model, see inverse_rosenblatt().

Simulated data is always a continous :math:`n \times d` matrix.

Parameter ``n``:
    number of observations.

Parameter ``qrng``:
    set to true for quasi-random numbers.

Parameter ``num_threads``:
    the number of threads to use for computations; if greater than 1,
    the function will generate ``n`` samples concurrently in
    ``num_threads`` batches.

Parameter ``seeds``:
    seeds of the random number generator; if empty (default), the
    random number generator is seeded randomly.

Returns:
    An :math:`n \times d` matrix of samples from the copula model.)""";
      } simulate;
      // Symbol: vinecopulib::Vinecop::str
      struct /* str */ {
        // Source: vinecopulib/vinecop/class.hpp:164
        const char* doc =
R"""(summarizes the model into a string (can be used for printing).)""";
      } str;
      // Symbol: vinecopulib::Vinecop::threshold_
      struct /* threshold_ */ {
        // Source: vinecopulib/vinecop/class.hpp:170
        const char* doc = R"""()""";
      } threshold_;
      // Symbol: vinecopulib::Vinecop::to_json
      struct /* to_json */ {
        // Source: vinecopulib/vinecop/class.hpp:60
        const char* doc =
R"""(write the copula object into a JSON file.

See to_ptree() for the structure of the file.

Parameter ``filename``:
    the name of the file to write.)""";
      } to_json;
      // Symbol: vinecopulib::Vinecop::to_ptree
      struct /* to_ptree */ {
        // Source: vinecopulib/vinecop/class.hpp:59
        const char* doc =
R"""(converts the copula into a boost::property_tree::ptree object.

The ``ptree`` object contains two nodes : ``"structure"`` for the vine
structure, which itself contains nodes ``"array"`` for the structure
triangular array and ``"order"`` for the order vector, and ``"pair
copulas"``. The former two encode the R-Vine structure and the latter
is a list of child nodes for the trees (``"tree1"``, `"tree2"`, etc),
each containing a list of child nodes for the edges (``"pc1"``,
`"pc2"`, etc). See Bicop::to_ptree() for the encoding of pair-copulas.

Returns:
    the boost::property_tree::ptree object containing the copula.)""";
      } to_ptree;
      // Symbol: vinecopulib::Vinecop::truncate
      struct /* truncate */ {
        // Source: vinecopulib/vinecop/class.hpp:162
        const char* doc =
R"""(truncate the vine copula model.

Parameter ``trunc_lvl``:
    the truncation level. If the model is already truncated at a level
    less than ``trunc_lvl``, the function does nothing.)""";
      } truncate;
      // Symbol: vinecopulib::Vinecop::var_types_
      struct /* var_types_ */ {
        // Source: vinecopulib/vinecop/class.hpp:173
        const char* doc = R"""()""";
      } var_types_;
    } Vinecop;
    // Symbol: vinecopulib::tools_select
    struct /* tools_select */ {
    } tools_select;
  } vinecopulib;
} mkdoc_doc;

#if defined(__GNUG__)
#pragma GCC diagnostic pop
#endif
