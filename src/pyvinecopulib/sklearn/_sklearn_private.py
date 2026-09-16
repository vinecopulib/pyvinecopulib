"""The scikit-learn internals the estimators depend on, in one place.

`_parameter_constraints` is the validation mechanism scikit-learn's
third-party estimator guide prescribes, and `VineBase` follows it -- but
the constraint vocabulary that mechanism reads ships **only** from
`sklearn.utils._param_validation`. There is no public spelling: neither
`sklearn.utils` nor `sklearn.exceptions` re-exports `Interval` or
`Options`, checked against 1.9.

So the dependency is real, and isolating it here is what makes it one
maintenance obligation rather than one per estimator module. Supported
range: scikit-learn 1.4 (the floor in `pyproject.toml`) through 1.9.
Nothing here is imported by tests -- the constraints are exercised through
the estimators, which is where a moved name surfaces.
"""

from sklearn.utils._param_validation import (  # noqa: PLC2701
  Interval,
  Options,
)

__all__ = ["Interval", "Options"]
