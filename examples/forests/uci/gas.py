import builtins
import importlib
import os
import pickle
import sys

import numpy as np
import pandas as pd
import pandas.core.indexes as idx_pkg
import pandas.core.indexes.base as idx_base
import pandas.core.indexes.category as idx_cat
import pandas.core.indexes.datetimes as idx_dt
import pandas.core.indexes.period as idx_period
import pandas.core.indexes.range as idx_range
import pandas.core.indexes.timedeltas as idx_td

MODULE_ALIASES = {
  "pandas.index": idx_base,
  "pandas.indexes": idx_pkg,
  "pandas.indexes.base": idx_base,
  "pandas.indexes.numeric": idx_base,
  "pandas.indexes.datetimes": idx_dt,
  "pandas.indexes.timedeltas": idx_td,
  "pandas.indexes.period": idx_period,
  "pandas.indexes.range": idx_range,
  "pandas.indexes.category": idx_cat,
  "pandas.core.index": idx_base,
  "pandas.core.indexes.base": idx_base,
  "pandas.core.indexes.numeric": idx_base,
  "pandas.core.indexes.datetimes": idx_dt,
  "pandas.core.indexes.timedeltas": idx_td,
  "pandas.core.indexes.period": idx_period,
  "pandas.core.indexes.range": idx_range,
  "pandas.core.indexes.category": idx_cat,
}

CLASS_ALIASES = {
  ("pandas.core.index", "Index"): pd.Index,
  ("pandas.indexes.base", "Index"): pd.Index,
  ("pandas.core.indexes.base", "Int64Index"): pd.Index,
  ("pandas.core.indexes.numeric", "Int64Index"): pd.Index,
  ("pandas.indexes.numeric", "Int64Index"): pd.Index,
  ("pandas.core.indexes.base", "UInt64Index"): pd.Index,
  ("pandas.core.indexes.numeric", "UInt64Index"): pd.Index,
  ("pandas.indexes.numeric", "UInt64Index"): pd.Index,
  ("pandas.core.indexes.base", "Float64Index"): pd.Index,
  ("pandas.core.indexes.numeric", "Float64Index"): pd.Index,
  ("pandas.indexes.numeric", "Float64Index"): pd.Index,
  ("pandas.core.indexes.base", "RangeIndex"): pd.RangeIndex,
  ("pandas.core.indexes.range", "RangeIndex"): pd.RangeIndex,
  ("pandas.indexes.range", "RangeIndex"): pd.RangeIndex,
  ("pandas.core.indexes.datetimes", "DatetimeIndex"): pd.DatetimeIndex,
  ("pandas.indexes.datetimes", "DatetimeIndex"): pd.DatetimeIndex,
  ("pandas.core.indexes.timedeltas", "TimedeltaIndex"): pd.TimedeltaIndex,
  ("pandas.indexes.timedeltas", "TimedeltaIndex"): pd.TimedeltaIndex,
  ("pandas.core.indexes.period", "PeriodIndex"): pd.PeriodIndex,
  ("pandas.indexes.period", "PeriodIndex"): pd.PeriodIndex,
  ("pandas.core.indexes.category", "CategoricalIndex"): pd.CategoricalIndex,
  ("pandas.indexes.category", "CategoricalIndex"): pd.CategoricalIndex,
}

NAME_FALLBACKS = {
  "Int64Index": pd.Index,
  "UInt64Index": pd.Index,
  "Float64Index": pd.Index,
  "RangeIndex": pd.RangeIndex,
  "DatetimeIndex": pd.DatetimeIndex,
  "TimedeltaIndex": pd.TimedeltaIndex,
  "PeriodIndex": pd.PeriodIndex,
  "CategoricalIndex": pd.CategoricalIndex,
  "Index": pd.Index,
}


class PandasCompatUnpickler(pickle.Unpickler):
  def __init__(self, file, **kwargs):
    # force latin1 decoding for Python2 pickles
    kwargs.setdefault("encoding", "latin1")
    super().__init__(file, **kwargs)

  def find_class(self, module, name):
    key = (module, name)
    if key in CLASS_ALIASES:
      return CLASS_ALIASES[key]
    mod_obj = MODULE_ALIASES.get(module)
    if mod_obj is not None:
      if hasattr(mod_obj, name):
        return getattr(mod_obj, name)
      if name in NAME_FALLBACKS:
        return NAME_FALLBACKS[name]
    try:
      mod = importlib.import_module(module)
      return getattr(mod, name)
    except (ModuleNotFoundError, AttributeError):
      if name in NAME_FALLBACKS:
        return NAME_FALLBACKS[name]
      raise


def load_legacy_pickle(path):
  # patch sys.modules only inside this function
  had_builtin = "__builtin__" in sys.modules
  backup_builtin = sys.modules.get("__builtin__")
  sys.modules["__builtin__"] = builtins
  try:
    with open(path, "rb") as f:
      return PandasCompatUnpickler(f).load()
  finally:
    if had_builtin:
      sys.modules["__builtin__"] = backup_builtin
    else:
      del sys.modules["__builtin__"]


class GAS:
  class Data:
    def __init__(self, data):
      self.x = data.astype(np.float32)
      self.N = self.x.shape[0]

  def __init__(self, file):
    trn, val, tst = load_data_and_clean_and_split(file)

    self.trn = self.Data(trn)
    self.val = self.Data(val)
    self.tst = self.Data(tst)

    self.n_dims = self.trn.x.shape[1]


def load_data(file):
  # data = pd.read_pickle(file)

  # Check if file exists and use .pickle instead of .csv if not
  if not os.path.exists(file):
    file = file.replace(".csv", ".pickle")
    data = load_legacy_pickle(file)
    file = file.replace(".pickle", ".csv")
    data.to_csv(file, index=False)

  data = pd.read_csv(file, index_col=False)
  # data = pd.read_pickle(file).sample(frac=0.25)
  # data.to_pickle(file)
  data.drop("Meth", axis=1, inplace=True)
  data.drop("Eth", axis=1, inplace=True)
  data.drop("Time", axis=1, inplace=True)
  return data


def get_correlation_numbers(data):
  C = data.corr()
  A = C > 0.98
  B = A.values.sum(axis=1)
  return B


def load_data_and_clean(file):
  data = load_data(file)
  B = get_correlation_numbers(data)

  while np.any(B > 1):
    col_to_remove = np.where(B > 1)[0][0]
    col_name = data.columns[col_to_remove]
    data.drop(col_name, axis=1, inplace=True)
    B = get_correlation_numbers(data)
  # print(data.corr())
  data = (data - data.mean()) / data.std()

  return data


def load_data_and_clean_and_split(file):
  data = load_data_and_clean(file).values
  N_test = int(0.1 * data.shape[0])
  data_test = data[-N_test:]
  data_train = data[0:-N_test]
  N_validate = int(0.1 * data_train.shape[0])
  data_validate = data_train[-N_validate:]
  data_train = data_train[0:-N_validate]

  return data_train, data_validate, data_test
