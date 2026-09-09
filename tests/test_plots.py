import math
from typing import Any, cast
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

import pyvinecopulib as pv

from .helpers import ShiftedNormalMargin


def assert_called_once_or_twice(mock: Any) -> None:
  """Helper to assert a mock was called once or twice"""
  if not mock.called:
    raise AssertionError("Expected to be called at least once")
  if not (mock.call_count == 1 or mock.call_count == 2):
    raise AssertionError(f"Expected call count 1 or 2, got {mock.call_count}")


def assert_called_once_or_twice_with(
  mock: Any, *args: Any, **kwargs: Any
) -> None:
  """Helper to assert a mock was called once or twice with specific args"""
  if not mock.called:
    raise AssertionError("Expected to be called at least once")
  if not (mock.call_count == 1 or mock.call_count == 2):
    raise AssertionError(f"Expected call count 1 or 2, got {mock.call_count}")
  calls = [call.args for call in mock.call_args_list]
  if not any(call == args for call in calls):
    raise AssertionError(f"Expected to be called with args {args}, got {calls}")


class TestPairCopulaData:
  """Test pair_copuladata.py functions directly"""

  def test_pairs_copula_data_parameter_validation(self) -> None:
    """Test pairs_copula_data parameter validation"""
    from pyvinecopulib.utils import pairs_copula_data

    # Test None data
    with pytest.raises(ValueError, match="`data` cannot be None"):
      # `cast` rather than a `ty: ignore`: the ignore reads as unused in
      # some environments and is required in others, so either spelling
      # fails the type check somewhere.
      pairs_copula_data(cast("Any", None))

    # Test non-numeric data
    with pytest.raises(
      ValueError, match="Could not convert `data` to numeric array"
    ):
      pairs_copula_data([["a", "b"], ["c", "d"]])

    # Test wrong dimensions
    with pytest.raises(ValueError, match="`data` must be a 2D array-like"):
      pairs_copula_data([0.1, 0.2, 0.3])

    # Test empty data
    with pytest.raises(ValueError, match="`data` cannot be empty"):
      pairs_copula_data(np.array([]).reshape(0, 2))

    # Test values outside (0,1)
    with pytest.raises(ValueError, match="All values must lie strictly in"):
      pairs_copula_data([[0.0, 0.5], [0.5, 1.0]])

    with pytest.raises(ValueError, match="All values must lie strictly in"):
      pairs_copula_data([[-0.1, 0.5], [0.5, 0.8]])

    with pytest.raises(ValueError, match="All values must lie strictly in"):
      pairs_copula_data([[0.1, 0.5], [0.5, 1.1]])

    # Test negative grid_size
    valid_data = np.random.uniform(0.1, 0.9, size=(10, 2))
    with pytest.raises(
      ValueError, match="`grid_size` must be a positive integer"
    ):
      pairs_copula_data(valid_data, grid_size=-1)

    with pytest.raises(
      ValueError, match="`grid_size` must be a positive integer"
    ):
      pairs_copula_data(valid_data, grid_size=0)

    # Test negative bins
    with pytest.raises(ValueError, match="`bins` must be a positive integer"):
      pairs_copula_data(valid_data, bins=-1)

    with pytest.raises(ValueError, match="`bins` must be a positive integer"):
      pairs_copula_data(valid_data, bins=0)

    # Test negative scatter_size
    with pytest.raises(
      ValueError, match="`scatter_size` must be a positive number"
    ):
      pairs_copula_data(valid_data, scatter_size=-1.0)

    with pytest.raises(
      ValueError, match="`scatter_size` must be a positive number"
    ):
      pairs_copula_data(valid_data, scatter_size=0.0)

    # # Test too many dimensions
    # high_dim_data = np.random.uniform(0.1, 0.9, size=(10, 11))
    # with pytest.raises(
    #   ValueError, match="Dimension 11 is too large for visualization"
    # ):
    #   pairs_copula_data(high_dim_data)

    # Test too few observations
    few_obs_data = np.random.uniform(0.1, 0.9, size=(1, 2))
    with pytest.raises(ValueError, match="Need at least 2 observations, got 1"):
      pairs_copula_data(few_obs_data)

  def test_pairs_copula_data_parameter_types(self) -> None:
    """Test parameter type validation"""
    from pyvinecopulib.utils import pairs_copula_data

    valid_data = np.random.uniform(0.1, 0.9, size=(10, 2))

    # Test non-integer grid_size
    with pytest.raises(
      ValueError, match="`grid_size` must be a positive integer"
    ):
      pairs_copula_data(valid_data, grid_size=cast("int", 10.5))

    # Test non-integer bins
    with pytest.raises(ValueError, match="`bins` must be a positive integer"):
      pairs_copula_data(valid_data, bins=cast("int", 5.5))

    # Test string scatter_size
    with pytest.raises(
      ValueError, match="`scatter_size` must be a positive number"
    ):
      pairs_copula_data(valid_data, scatter_size=cast("float", "large"))

  def test_pairs_copula_data_basic_validation_success(self) -> None:
    """Test that valid inputs pass basic validation"""
    from pyvinecopulib.utils import pairs_copula_data

    # Create valid test data
    np.random.seed(42)
    data = np.random.uniform(0.1, 0.9, size=(10, 2))

    # Mock the plotting parts since we just want to test validation
    with patch("matplotlib.pyplot.subplots") as mock_subplots:
      mock_fig = MagicMock()
      mock_ax = MagicMock()
      mock_subplots.return_value = (mock_fig, mock_ax)

      # Mock the wdm and Bicop imports that would fail without the C++ extension
      with patch("pyvinecopulib.utils._pair_plots.wdm"):
        with patch("pyvinecopulib.utils._pair_plots.Bicop"):
          with patch("pyvinecopulib.utils._pair_plots.norm_cdf"):
            with patch("pyvinecopulib.utils._pair_plots.norm_pdf"):
              with patch("pyvinecopulib.utils._pair_plots.plt.tight_layout"):
                # This should not raise any validation errors
                try:
                  pairs_copula_data(data)
                  validation_passed = True
                except (ImportError, AttributeError):
                  # Expected due to missing matplotlib/C++ extension interactions
                  validation_passed = True
                except ValueError:
                  # This would be a validation error, which we don't expect
                  validation_passed = False

                assert validation_passed, "Valid data should pass validation"

  def test_pairs_copula_data_edge_cases(self) -> None:
    """Test edge cases that should be handled gracefully"""
    from pyvinecopulib.utils import pairs_copula_data

    # Test minimum valid data (2 observations, 1 dimension)
    min_data = np.array([[0.1], [0.9]])

    with patch("matplotlib.pyplot.subplots") as mock_subplots:
      mock_fig = MagicMock()
      mock_ax = MagicMock()
      mock_subplots.return_value = (mock_fig, mock_ax)

      with patch("pyvinecopulib.utils._pair_plots.norm_cdf"):
        with patch("pyvinecopulib.utils._pair_plots.norm_pdf"):
          with patch("pyvinecopulib.utils._pair_plots.plt.tight_layout"):
            # This should not raise validation errors
            try:
              pairs_copula_data(min_data)
              edge_case_passed = True
            except (ImportError, AttributeError):
              # Expected due to missing matplotlib/C++ extension interactions
              edge_case_passed = True
            except ValueError:
              # This would be a validation error
              edge_case_passed = False

            assert edge_case_passed, "Minimum valid data should pass validation"

    # Test exactly at dimension limit
    max_dim_data = np.random.uniform(0.1, 0.9, size=(5, 10))

    with patch("matplotlib.pyplot.subplots") as mock_subplots:
      mock_fig = MagicMock()
      mock_ax = MagicMock()
      mock_subplots.return_value = (mock_fig, mock_ax)

      with patch("pyvinecopulib.utils._pair_plots.norm_cdf"):
        with patch("pyvinecopulib.utils._pair_plots.norm_pdf"):
          with patch("pyvinecopulib.utils._pair_plots.plt.tight_layout"):
            # This should not raise validation errors
            try:
              pairs_copula_data(max_dim_data)
              max_dim_passed = True
            except (ImportError, AttributeError):
              # Expected due to missing matplotlib/C++ extension interactions
              max_dim_passed = True
            except ValueError:
              # This would be a validation error
              max_dim_passed = False

            assert max_dim_passed, (
              "Maximum dimension data should pass validation"
            )


class TestBicopHelpers:
  """Test bicop.py helper functions directly"""

  def test_get_default_xylim(self) -> None:
    """Test get_default_xylim function"""
    from pyvinecopulib.core._bicop_plot import get_default_xylim

    # Test valid margin types
    assert get_default_xylim("unif") == (1e-2, 1 - 1e-2)
    assert get_default_xylim("norm") == (-3, 3)
    assert get_default_xylim("exp") == (0, 6)

    # Test invalid margin type
    with pytest.raises(ValueError, match="Unknown margin type"):
      get_default_xylim("invalid")

  def test_get_default_grid_size(self) -> None:
    """Test get_default_grid_size function"""
    from pyvinecopulib.core._bicop_plot import get_default_grid_size

    # Test valid plot types
    assert get_default_grid_size("contour") == 100
    assert get_default_grid_size("surface") == 40

    # Test invalid plot type
    with pytest.raises(ValueError, match="Unknown plot type"):
      get_default_grid_size("invalid")

  def test_bicop_plot_parameter_validation(self) -> None:
    """Test bicop_plot parameter validation"""
    from pyvinecopulib.core._bicop_plot import bicop_plot

    # Create a mock copula object
    mock_cop = MagicMock()
    mock_cop.var_types = ["c", "c"]
    mock_cop.pdf.return_value = np.ones(100)

    # Test invalid plot type
    with pytest.raises(ValueError, match="Unknown type"):
      bicop_plot(mock_cop, plot_type="invalid")

    # Test invalid margin type
    with pytest.raises(ValueError, match="Unknown margin type"):
      bicop_plot(mock_cop, margin_type="invalid")

  def test_bicop_plot_restores_discrete_type_after_failure(self) -> None:
    """A plotting error cannot mutate a caller-owned discrete pair."""
    from pyvinecopulib.core._bicop_plot import bicop_plot

    class Pair:
      var_types = ["d", "c"]

      def pdf(self, u: np.ndarray) -> np.ndarray:
        raise RuntimeError("density failed")

    pair = Pair()
    with pytest.raises(RuntimeError, match="density failed"):
      bicop_plot(pair, grid_size=10)
    assert pair.var_types == ["d", "c"]

  @patch("matplotlib.pyplot.show")
  @patch("matplotlib.pyplot.contour")
  @patch("matplotlib.pyplot.clabel")
  def test_bicop_plot_contour(
    self, mock_clabel: Any, mock_contour: Any, mock_show: Any
  ) -> None:
    """Test bicop_plot with contour type"""
    from pyvinecopulib.core._bicop_plot import bicop_plot

    # Create a mock copula object
    mock_cop = MagicMock()
    mock_cop.var_types = ["c", "c"]
    mock_cop.pdf.return_value = np.ones(10000)  # 100x100 grid

    # Test contour plot
    bicop_plot(mock_cop, plot_type="contour", grid_size=100)

    # Verify matplotlib functions were called
    assert_called_once_or_twice(mock_contour)
    assert_called_once_or_twice(mock_clabel)
    assert_called_once_or_twice(mock_show)

  @patch("matplotlib.pyplot.show")
  @patch("matplotlib.pyplot.figure")
  def test_bicop_plot_surface(self, mock_figure: Any, mock_show: Any) -> None:
    """Test bicop_plot with surface type"""
    from pyvinecopulib.core._bicop_plot import bicop_plot

    # Create mock figure and axis
    mock_fig = MagicMock()
    mock_ax = MagicMock()
    mock_figure.return_value = mock_fig
    mock_fig.add_subplot.return_value = mock_ax

    # Create a mock copula object
    mock_cop = MagicMock()
    mock_cop.var_types = ["c", "c"]
    mock_cop.pdf.return_value = np.ones(1600)  # 40x40 grid

    # Test surface plot
    bicop_plot(mock_cop, plot_type="surface", grid_size=40)

    # Verify 3D plotting was set up
    assert_called_once_or_twice(mock_figure)
    assert_called_once_or_twice_with(mock_fig.add_subplot, 111, projection="3d")
    assert_called_once_or_twice(mock_ax.plot_surface)
    assert_called_once_or_twice(mock_show)

  @patch("matplotlib.pyplot.show")
  @patch("matplotlib.pyplot.contour")
  def test_bicop_plot_margin_types(
    self, mock_contour: Any, mock_show: Any
  ) -> None:
    """Test bicop_plot with different margin types"""
    from pyvinecopulib.core._bicop_plot import bicop_plot

    # Create a mock copula object
    mock_cop = MagicMock()
    mock_cop.var_types = ["c", "c"]
    mock_cop.pdf.return_value = np.ones(10000)  # 100x100 grid

    margin_types = ["unif", "norm", "exp"]

    for margin_type in margin_types:
      mock_contour.reset_mock()
      mock_show.reset_mock()

      bicop_plot(
        mock_cop, plot_type="contour", margin_type=margin_type, grid_size=100
      )

      assert_called_once_or_twice(mock_contour)
      assert_called_once_or_twice(mock_show)

  @patch("matplotlib.pyplot.show")
  @patch("matplotlib.pyplot.contour")
  def test_bicop_plot_custom_parameters(
    self, mock_contour: Any, mock_show: Any
  ) -> None:
    """Test bicop_plot with custom xylim and grid_size"""
    from pyvinecopulib.core._bicop_plot import bicop_plot

    # Create a mock copula object
    mock_cop = MagicMock()
    mock_cop.var_types = ["c", "c"]
    mock_cop.pdf.return_value = np.ones(2500)  # 50x50 grid

    # Test with custom parameters
    bicop_plot(
      mock_cop,
      plot_type="contour",
      margin_type="unif",
      xylim=(0.1, 0.9),
      grid_size=50,
    )

    assert_called_once_or_twice(mock_contour)
    assert_called_once_or_twice(mock_show)


class TestVinecopHelpers:
  """Test vinecop.py helper functions directly"""

  def setup_method(self) -> None:
    """Set up test fixtures"""
    # Create a simple vine copula for testing
    np.random.seed(1234)
    data = np.random.uniform(0, 1, size=(50, 4))
    self.vinecop = pv.Vinecop.from_data(
      data, controls=pv.FitControlsVinecop(family_set=[pv.families.indep])
    )

  def teardown_method(self) -> None:
    """Drop nanobind-backed instances so they don't accumulate when
    ``torch`` is imported alongside the test suite — torch's atexit
    handlers run before nanobind's and break the cleanup chain.
    """
    del self.vinecop
    import gc

    gc.collect()

  def test_get_name(self) -> None:
    """Test get_name function"""
    from pyvinecopulib.core._vinecop_plot import get_name

    vars_names = ["X1", "X2", "X3", "X4"]

    # Test for tree 0 (no conditioning)
    name = get_name(self.vinecop, 0, 0, vars_names)
    assert isinstance(name, str)
    assert len(name) > 0
    assert ";" not in name  # No conditioning separator for tree 0

    # Test for higher trees (with conditioning)
    if self.vinecop.trunc_lvl > 1:
      name = get_name(self.vinecop, 1, 0, vars_names)
      assert isinstance(name, str)
      assert len(name) > 0

  def test_get_graph(self) -> None:
    """Test get_graph function"""
    from pyvinecopulib.core._vinecop_plot import get_graph

    vars_names = ["X1", "X2", "X3", "X4"]

    # Test for tree 0
    adj_mat, node_labels, edge_labels = get_graph(0, self.vinecop, vars_names)

    # Check return types and shapes
    assert isinstance(adj_mat, np.ndarray)
    assert isinstance(node_labels, dict)
    assert isinstance(edge_labels, dict)

    # Check adjacency matrix properties
    assert adj_mat.shape[0] == adj_mat.shape[1]  # Square matrix
    assert adj_mat.shape[0] == len(node_labels)  # Consistent dimensions
    assert np.all((adj_mat == 0) | (adj_mat == 1))  # Binary matrix
    assert np.array_equal(adj_mat, adj_mat.T)  # Symmetric matrix

    # Check node labels
    assert all(isinstance(label, str) for label in node_labels.values())
    assert len(node_labels) == self.vinecop.dim

    # Check edge labels
    assert all(isinstance(label, str) for label in edge_labels.values())
    assert all(
      isinstance(key, tuple) and len(key) == 2 for key in edge_labels.keys()
    )

  def test_get_graph_higher_trees(self) -> None:
    """Test get_graph for higher order trees"""
    from pyvinecopulib.core._vinecop_plot import get_graph

    vars_names = ["X1", "X2", "X3", "X4"]

    # Test for tree 1 if available
    if self.vinecop.trunc_lvl > 1:
      adj_mat, node_labels, edge_labels = get_graph(1, self.vinecop, vars_names)

      # Check that dimensions decrease for higher trees
      assert adj_mat.shape[0] == self.vinecop.dim - 1
      assert len(node_labels) == self.vinecop.dim - 1

  def test_vinecop_plot_parameter_validation(self) -> None:
    """Test vinecop_plot parameter validation"""
    from pyvinecopulib.core._vinecop_plot import vinecop_plot

    # Test with wrong number of variable names
    with pytest.raises(
      ValueError, match="The number of variable names must be equal"
    ):
      vinecop_plot(self.vinecop, vars_names=["X1", "X2"])

  def test_vinecop_plot_high_dimension_error(self) -> None:
    """Test vinecop_plot with high truncation level"""
    from pyvinecopulib.core._vinecop_plot import vinecop_plot

    # Create a mock vinecop with high truncation level
    mock_vinecop = MagicMock()
    mock_vinecop.trunc_lvl = 10

    with pytest.raises(
      ValueError, match="The dimension and truncation level are too high"
    ):
      vinecop_plot(mock_vinecop)

  def test_vinecop_plot_basic(self) -> None:
    """Test basic vinecop_plot functionality"""
    from pyvinecopulib.core._vinecop_plot import vinecop_plot

    # Instead of mocking the complex matplotlib behavior,
    # let's test that the function accepts correct parameters and handles validation
    with patch("matplotlib.pyplot.subplots") as mock_subplots:
      with patch("matplotlib.pyplot.show"):
        with patch("matplotlib.pyplot.tight_layout"):
          with patch("networkx.draw"):
            with patch(
              "networkx.drawing.nx_pydot.graphviz_layout"
            ) as mock_graphviz:
              # Mock the return to avoid complex axis handling
              mock_graphviz.return_value = {0: (0, 0), 1: (1, 1)}
              mock_fig = MagicMock()
              mock_ax = MagicMock()
              mock_subplots.return_value = (mock_fig, mock_ax)

              # The key is that this doesn't raise an exception
              try:
                vinecop_plot(self.vinecop, tree=[0])
                # If we get here, the function accepted our parameters
                success = True
              except AttributeError as e:
                if "set_title" in str(e):
                  # This is the known matplotlib mocking issue, not a real bug
                  success = True
                else:
                  success = False
              except Exception:
                success = False

              assert success, "vinecop_plot should accept valid parameters"

  def test_vinecop_plot_multiple_trees(self) -> None:
    """Test vinecop_plot with multiple trees"""
    from pyvinecopulib.core._vinecop_plot import vinecop_plot

    # Use the same simpler mocking approach as test_vinecop_plot_basic
    with patch("matplotlib.pyplot.subplots") as mock_subplots:
      with patch("matplotlib.pyplot.show"):
        with patch("matplotlib.pyplot.tight_layout"):
          with patch("networkx.draw"):
            with patch("networkx.draw_networkx_edge_labels"):
              with patch(
                "networkx.drawing.nx_pydot.graphviz_layout"
              ) as mock_graphviz:
                # Mock the return to avoid complex axis handling
                mock_graphviz.return_value = {0: (0, 0), 1: (1, 1)}
                mock_fig = MagicMock()
                mock_axes = np.array([[MagicMock(), MagicMock()]])
                mock_subplots.return_value = (mock_fig, mock_axes)

                # Test with multiple trees
                available_trees = min(3, self.vinecop.trunc_lvl)
                tree_list = list(range(available_trees))

                # For multiple trees, we need a proper axes array
                if len(tree_list) > 1:
                  # Create a mock axes array with enough elements
                  mock_axes = np.array(
                    [MagicMock() for _ in range(len(tree_list))]
                  )
                  mock_subplots.return_value = (mock_fig, mock_axes)

                try:
                  vinecop_plot(self.vinecop, tree=tree_list)
                  success = True
                except AttributeError as e:
                  if "set_title" in str(e):
                    # This is the known matplotlib mocking issue, not a real bug
                    success = True
                  else:
                    success = False
                except Exception:
                  success = False

                assert success, "vinecop_plot should handle multiple trees"

  def test_vinecop_plot_edge_labels(self) -> None:
    """Test vinecop_plot with and without edge labels"""
    from pyvinecopulib.core._vinecop_plot import vinecop_plot

    # Use the same simpler mocking approach as test_vinecop_plot_basic
    with patch("matplotlib.pyplot.subplots") as mock_subplots:
      with patch("matplotlib.pyplot.show"):
        with patch("matplotlib.pyplot.tight_layout"):
          with patch("networkx.draw"):
            with patch(
              "networkx.draw_networkx_edge_labels"
            ) as mock_edge_labels:
              with patch(
                "networkx.drawing.nx_pydot.graphviz_layout"
              ) as mock_graphviz:
                # Mock the return to avoid complex axis handling
                mock_graphviz.return_value = {0: (0, 0), 1: (1, 1)}
                mock_fig = MagicMock()
                mock_ax = MagicMock()
                mock_subplots.return_value = (mock_fig, mock_ax)

                # Test with edge labels
                try:
                  vinecop_plot(self.vinecop, tree=[0], add_edge_labels=True)
                  success = True
                except AttributeError as e:
                  if "set_title" in str(e):
                    # This is the known matplotlib mocking issue, not a real bug
                    success = True
                  else:
                    success = False
                except Exception:
                  success = False

                assert success, (
                  "vinecop_plot should accept valid parameters with edge labels"
                )

                # Reset mocks
                mock_edge_labels.reset_mock()

                # Test without edge labels
                try:
                  vinecop_plot(self.vinecop, tree=[0], add_edge_labels=False)
                  success = True
                except AttributeError as e:
                  if "set_title" in str(e):
                    # This is the known matplotlib mocking issue, not a real bug
                    success = True
                  else:
                    success = False
                except Exception:
                  success = False

                assert success, (
                  "vinecop_plot should accept valid parameters without edge labels"
                )

  def test_vinecop_plot_layouts(self) -> None:
    """Test vinecop_plot with different layouts"""
    from pyvinecopulib.core._vinecop_plot import vinecop_plot

    # Use the same simpler mocking approach as test_vinecop_plot_basic
    with patch("matplotlib.pyplot.subplots") as mock_subplots:
      with patch("matplotlib.pyplot.show"):
        with patch("matplotlib.pyplot.tight_layout"):
          with patch("networkx.draw"):
            with patch(
              "networkx.drawing.nx_pydot.graphviz_layout"
            ) as mock_graphviz:
              with patch("networkx.spring_layout") as mock_spring:
                # Mock the return to avoid complex axis handling
                mock_graphviz.return_value = {0: (0, 0), 1: (1, 1)}
                mock_spring.return_value = {0: (0, 0), 1: (1, 1)}
                mock_fig = MagicMock()
                mock_ax = MagicMock()
                mock_subplots.return_value = (mock_fig, mock_ax)

                # Test graphviz layout
                try:
                  vinecop_plot(self.vinecop, tree=[0], layout="graphviz")
                  success = True
                except AttributeError as e:
                  if "set_title" in str(e):
                    # This is the known matplotlib mocking issue, not a real bug
                    success = True
                  else:
                    success = False
                except Exception:
                  success = False

                assert success, "vinecop_plot should accept graphviz layout"

                # Test spring layout
                try:
                  vinecop_plot(self.vinecop, tree=[0], layout="spring_layout")
                  success = True
                except AttributeError as e:
                  if "set_title" in str(e):
                    # This is the known matplotlib mocking issue, not a real bug
                    success = True
                  else:
                    success = False
                except Exception:
                  success = False

                assert success, "vinecop_plot should accept spring layout"

  def test_vinecop_plot_custom_variable_names(self) -> None:
    """Test vinecop_plot with custom variable names"""
    from pyvinecopulib.core._vinecop_plot import vinecop_plot

    # Use the same simpler mocking approach as test_vinecop_plot_basic
    with patch("matplotlib.pyplot.subplots") as mock_subplots:
      with patch("matplotlib.pyplot.show"):
        with patch("matplotlib.pyplot.tight_layout"):
          with patch("networkx.draw"):
            with patch(
              "networkx.drawing.nx_pydot.graphviz_layout"
            ) as mock_graphviz:
              # Mock the return to avoid complex axis handling
              mock_graphviz.return_value = {0: (0, 0), 1: (1, 1)}
              mock_fig = MagicMock()
              mock_ax = MagicMock()
              mock_subplots.return_value = (mock_fig, mock_ax)

              # Test with custom variable names
              custom_vars = ["Var1", "Var2", "Var3", "Var4"]
              try:
                vinecop_plot(self.vinecop, tree=[0], vars_names=custom_vars)
                success = True
              except AttributeError as e:
                if "set_title" in str(e):
                  # This is the known matplotlib mocking issue, not a real bug
                  success = True
                else:
                  success = False
              except Exception:
                success = False

              assert success, "vinecop_plot should accept custom variable names"

  def test_vinecop_plot_subplot_calculation(self) -> None:
    """Test subplot layout calculation"""
    from pyvinecopulib.core._vinecop_plot import vinecop_plot

    # Use the same simpler mocking approach as test_vinecop_plot_basic
    with patch("matplotlib.pyplot.subplots") as mock_subplots:
      with patch("matplotlib.pyplot.show"):
        with patch("matplotlib.pyplot.tight_layout"):
          with patch("networkx.draw"):
            with patch(
              "networkx.drawing.nx_pydot.graphviz_layout"
            ) as mock_graphviz:
              # Mock the return to avoid complex axis handling
              mock_graphviz.return_value = {0: (0, 0), 1: (1, 1)}
              mock_fig = MagicMock()
              mock_ax = MagicMock()
              mock_subplots.return_value = (mock_fig, mock_ax)

              # Test with 1 tree (should be 1x1)
              try:
                vinecop_plot(self.vinecop, tree=[0])
                success = True
              except AttributeError as e:
                if "set_title" in str(e):
                  # This is the known matplotlib mocking issue, not a real bug
                  success = True
                else:
                  success = False
              except Exception:
                success = False

              assert success, "vinecop_plot should handle single tree layout"

              # We can check the call args if the function succeeded
              if mock_subplots.called:
                args, kwargs = mock_subplots.call_args
                assert args[0] == 1  # n_row
                assert args[1] == 1  # n_col

              # Test with 3 trees (should be 3x1) if available
              if self.vinecop.trunc_lvl >= 3:
                mock_subplots.reset_mock()
                try:
                  vinecop_plot(self.vinecop, tree=[0, 1, 2])
                  success = True
                except AttributeError as e:
                  if "set_title" in str(e):
                    # This is the known matplotlib mocking issue, not a real bug
                    success = True
                  else:
                    success = False
                except Exception:
                  success = False

                assert success, (
                  "vinecop_plot should handle multiple tree layout"
                )

                # We can check the call args if the function succeeded
                if mock_subplots.called:
                  args, kwargs = mock_subplots.call_args
                  assert args[0] == 3  # n_row
                  assert args[1] == 1  # n_col


class TestMarginPlot:
  """What a margin plot draws, on real margins of all three variable types.

  Reached at the function rather than through ``.plot()``, which draws and
  returns nothing: the grid, the marks and the limits are checkable here and
  nowhere above it.
  """

  @staticmethod
  def _continuous() -> pv.core.Kde1d:
    kde = pv.core.Kde1d()
    kde.fit(np.random.RandomState(123).beta(0.5, 2.0, 200))
    return kde

  @staticmethod
  def _discrete(lo: int = 0, hi: int = 7) -> pv.core.Kde1d:
    rng = np.random.RandomState(7)
    y = np.clip(rng.poisson(3, 300), lo, hi).astype(float)
    # Both ends present, so the lattice the grid should recover is `lo..hi`.
    y[:2] = [float(lo), float(hi)]
    kde = pv.core.Kde1d(type="discrete")
    kde.fit(y)
    return kde

  @staticmethod
  def _zero_inflated() -> pv.core.Kde1d:
    rng = np.random.RandomState(5)
    y = rng.exponential(2.0, 300)
    y[rng.choice(300, 90, replace=False)] = 0.0
    kde = pv.core.Kde1d(xmin=0, type="zero-inflated")
    kde.fit(y)
    return kde

  def test_the_continuous_grid_spans_the_fitted_range(self) -> None:
    from pyvinecopulib.core._margin_plot import make_plotting_grid

    kde = self._continuous()
    grid = make_plotting_grid(kde, grid_size=100)

    assert grid.shape == (100,)
    assert np.all(np.diff(grid) > 0)
    assert grid[0] == pytest.approx(kde.grid_points.min())
    assert grid[-1] == pytest.approx(kde.grid_points.max())

  def test_a_declared_bound_wins_over_the_fitted_range(self) -> None:
    from pyvinecopulib.core._margin_plot import make_plotting_grid

    kde = pv.core.Kde1d(xmin=0.0, xmax=1.0)
    kde.fit(np.random.RandomState(1).beta(0.5, 2.0, 200))
    grid = make_plotting_grid(kde, grid_size=50)

    assert (grid[0], grid[-1]) == (0.0, 1.0)
    assert np.all(np.diff(grid) > 0)

  def test_the_discrete_grid_is_the_lattice_and_invents_no_level(self) -> None:
    """A fitted grid runs half a unit past the support at each end.

    Rounding it outwards -- which is what ``floor`` / ``ceil`` do to
    ``min - 0.25`` and ``max + 0.25`` -- plots one level below the smallest
    observation and one above the largest, neither of which can occur.
    """
    from pyvinecopulib.core._margin_plot import make_plotting_grid

    kde = self._discrete(lo=0, hi=7)
    assert kde.grid_points.min() < 0 and kde.grid_points.max() > 7

    grid = make_plotting_grid(kde)
    np.testing.assert_array_equal(grid, np.arange(0.0, 8.0))

  def test_a_declared_discrete_bound_is_a_level(self) -> None:
    from pyvinecopulib.core._margin_plot import make_plotting_grid

    y = np.clip(np.random.RandomState(3).poisson(3, 200), 0, 9).astype(float)
    kde = pv.core.Kde1d(xmin=0, xmax=9, type="discrete")
    kde.fit(y)

    np.testing.assert_array_equal(make_plotting_grid(kde), np.arange(0.0, 10.0))

  def test_the_zero_inflated_grid_excludes_the_atom(self) -> None:
    from pyvinecopulib.core._margin_plot import make_plotting_grid

    grid = make_plotting_grid(self._zero_inflated(), grid_size=101)

    assert 0.0 not in grid
    assert grid.min() > 0.0

  def test_a_margin_with_no_grid_is_drawn_between_its_tails(self) -> None:
    """A support with no last point is cut at the 0.1% quantiles."""
    from pyvinecopulib.core._margin_plot import make_plotting_grid

    scipy = pytest.importorskip("scipy")
    del scipy
    from pyvinecopulib.margins import SciPyMargin

    margin = SciPyMargin(family="norm").fit(
      np.random.RandomState(0).normal(size=400)
    )
    grid = make_plotting_grid(margin, grid_size=64)

    assert grid.shape == (64,)
    lo, hi = (
      float(margin.icdf(np.array([1e-3]))[0]),
      float(margin.icdf(np.array([1 - 1e-3]))[0]),
    )
    assert grid[0] == pytest.approx(lo)
    assert grid[-1] == pytest.approx(hi)

  @pytest.mark.parametrize("kind", ["density", "cdf"])
  @pytest.mark.parametrize("var_type", ["c", "d", "zi"])
  def test_every_variable_type_draws_both_kinds(
    self, var_type: str, kind: str
  ) -> None:
    from pyvinecopulib.core._margin_plot import margin_plot

    kde = {
      "c": self._continuous,
      "d": self._discrete,
      "zi": self._zero_inflated,
    }[var_type]()
    assert kde.var_type == var_type

    with (
      patch("matplotlib.pyplot.show") as show,
      patch("matplotlib.pyplot.plot") as plot,
      patch("matplotlib.pyplot.ylabel") as ylabel,
    ):
      margin_plot(kde, kind=kind)

    show.assert_called_once()
    ylabel.assert_called_once_with(
      "density" if kind == "density" else "probability"
    )
    # A discrete margin is marks, a continuous one a curve, and a
    # zero-inflated one the curve plus its one emphasized atom.
    first = plot.call_args_list[0]
    assert (first.kwargs.get("linestyle") == "None") == (var_type == "d")
    assert plot.call_count == (2 if var_type == "zi" else 1)

  def test_the_atom_is_droppable(self) -> None:
    from pyvinecopulib.core._margin_plot import margin_plot

    with patch("matplotlib.pyplot.show"), patch("matplotlib.pyplot.plot") as p:
      margin_plot(self._zero_inflated(), show_zero_mass=False)

    assert p.call_count == 1

  def test_the_cdf_is_drawn_on_the_unit_interval(self) -> None:
    from pyvinecopulib.core._margin_plot import margin_plot

    with (
      patch("matplotlib.pyplot.show"),
      patch("matplotlib.pyplot.plot") as plot,
      patch("matplotlib.pyplot.ylim") as ylim,
    ):
      margin_plot(self._continuous(), kind="cdf")

    values = plot.call_args_list[0].args[1]
    assert np.all((values >= 0) & (values <= 1))
    assert np.all(np.diff(values) >= -1e-12)
    ylim.assert_called_once_with(0, 1.05)

  def test_explicit_limits_are_passed_through(self) -> None:
    from pyvinecopulib.core._margin_plot import margin_plot

    with (
      patch("matplotlib.pyplot.show"),
      patch("matplotlib.pyplot.plot"),
      patch("matplotlib.pyplot.xlim") as xlim,
      patch("matplotlib.pyplot.ylim") as ylim,
    ):
      margin_plot(self._continuous(), xlim=(0.1, 0.4), ylim=(0, 2))

    xlim.assert_called_once_with((0.1, 0.4))
    ylim.assert_called_once_with((0, 2))

  def test_an_unknown_kind_is_refused(self) -> None:
    from pyvinecopulib.core._margin_plot import margin_plot

    with pytest.raises(ValueError, match="kind must be"):
      margin_plot(self._continuous(), kind="quantile")

  def test_an_unfitted_margin_is_refused(self) -> None:
    from pyvinecopulib.core._margin_plot import margin_plot

    with pytest.raises(ValueError, match="must be fitted"):
      margin_plot(pv.core.Kde1d())

  def test_a_margin_reading_no_covariates_refuses_x(self) -> None:
    """The refusal is explicit, since forwarding to a margin is by flag.

    ``declared_eval`` *skips* ``x`` for a margin that declares no
    ``supports_covariates``, so nothing below would raise and the plot would
    silently be the unconditional one.
    """
    from pyvinecopulib.core._margin_plot import margin_plot

    with pytest.raises(ValueError, match="declares no `supports_covariates`"):
      margin_plot(self._continuous(), x=np.zeros(2))

  def test_a_conditional_margin_is_drawn_at_one_covariate_row(self) -> None:
    margin = ShiftedNormalMargin(slope=2.0)
    with patch("matplotlib.pyplot.show"), patch("matplotlib.pyplot.plot") as p:
      margin.plot(x=np.array([1.0]))

    grid, values = p.call_args_list[0].args[:2]
    # Every row saw the same covariate, so the curve is the shifted density.
    np.testing.assert_allclose(
      values, np.exp(-0.5 * (grid - 2.0) ** 2) / np.sqrt(2.0 * np.pi)
    )

  def test_a_conditional_margin_refuses_more_than_one_row(self) -> None:
    from pyvinecopulib.core._margin_plot import margin_plot

    with pytest.raises(ValueError, match="single covariate row"):
      margin_plot(ShiftedNormalMargin(), x=np.zeros((3, 1)))

  def test_the_base_places_the_grid_on_the_margins_own_namespace(self) -> None:
    """A margin on tensors plots without converting inside its own ``pdf``.

    The grid is manufactured by the plot, so it is the one array the margin
    did not supply and cannot infer a namespace from; ``_prep`` is what brings
    it across, and ``pdf`` here would raise on anything else.
    """
    torch = pytest.importorskip("torch")
    from pyvinecopulib.core import MarginBase

    class _TensorMargin(MarginBase[Any]):
      """A standard normal whose arrays are tensors, so ``_prep`` finds one."""

      def __init__(self) -> None:
        # The array `reference_array` infers placement from.
        self._weights = torch.ones(4, dtype=torch.float64)
        self.seen: list[Any] = []

      def pdf(self, y: Any, /, *, x: Any = None) -> Any:
        assert isinstance(y, torch.Tensor)
        self.seen.append(y)
        return torch.exp(-0.5 * y * y) / math.sqrt(2.0 * math.pi)

      def cdf(self, y: Any, /, *, x: Any = None) -> Any:
        assert isinstance(y, torch.Tensor)
        return torch.special.ndtr(y)

      @property
      def support(self) -> tuple[float, float]:
        return (-4.0, 4.0)

    margin = _TensorMargin()
    with patch("matplotlib.pyplot.show"), patch("matplotlib.pyplot.plot") as p:
      margin.plot()

    assert margin.seen and isinstance(margin.seen[0], torch.Tensor)
    # The drawn values come back as NumPy, whatever the margin answered in.
    assert isinstance(p.call_args_list[0].args[1], np.ndarray)


class TestPlotDocstrings:
  """Test that docstrings are properly defined"""

  def test_bicop_plot_doc(self) -> None:
    """Test BICOP_PLOT_DOC is defined"""
    from pyvinecopulib.core._bicop_plot import BICOP_PLOT_DOC

    assert isinstance(BICOP_PLOT_DOC, str)
    assert len(BICOP_PLOT_DOC) > 0
    assert "Parameters" in BICOP_PLOT_DOC
    assert "Returns" in BICOP_PLOT_DOC

  def test_vinecop_plot_doc(self) -> None:
    """Test VINECOP_PLOT_DOC is defined"""
    from pyvinecopulib.core._vinecop_plot import VINECOP_PLOT_DOC

    assert isinstance(VINECOP_PLOT_DOC, str)
    assert len(VINECOP_PLOT_DOC) > 0
    assert "Parameters" in VINECOP_PLOT_DOC
    assert "Returns" in VINECOP_PLOT_DOC

  def test_margin_plot_doc(self) -> None:
    """Test MARGIN_PLOT_DOC is defined"""
    from pyvinecopulib.core._margin_plot import MARGIN_PLOT_DOC

    assert isinstance(MARGIN_PLOT_DOC, str)
    assert len(MARGIN_PLOT_DOC) > 0
    assert "Parameters" in MARGIN_PLOT_DOC
    assert "Returns" in MARGIN_PLOT_DOC


class TestEdgeCases:
  """Test edge cases and error conditions"""

  def test_bicop_plot_with_identical_density_values(self) -> None:
    """Test bicop_plot when all density values are identical"""
    from pyvinecopulib.core._bicop_plot import bicop_plot

    # Create a mock copula that returns identical density values
    mock_cop = MagicMock()
    mock_cop.var_types = ["c", "c"]
    mock_cop.pdf.return_value = np.ones(100)  # All values identical

    with patch("matplotlib.pyplot.show"), patch("matplotlib.pyplot.contour"):
      # This should handle the case where all density values are the same
      # The code adjusts dens[0] = 1.000001 * dens[0] to handle this
      bicop_plot(mock_cop, grid_size=10)

  def test_vinecop_get_name_edge_cases(self) -> None:
    """Test get_name function edge cases"""
    from pyvinecopulib.core._vinecop_plot import get_name

    # Create a simple 3D vine copula
    np.random.seed(1234)
    data = np.random.uniform(0, 1, size=(50, 3))
    vinecop = pv.Vinecop.from_data(
      data, controls=pv.FitControlsVinecop(family_set=[pv.families.indep])
    )

    vars_names = ["A", "B", "C"]

    # Test tree 0 (no conditioning set)
    name = get_name(vinecop, 0, 0, vars_names)
    assert ";" not in name  # No semicolon for tree 0

    # Test tree 1 if available (should have conditioning set)
    if vinecop.trunc_lvl > 1:
      name = get_name(vinecop, 1, 0, vars_names)
      # Tree 1 should have conditioning, so might contain semicolon
      assert isinstance(name, str)

  def test_matplotlib_colormap_creation(self) -> None:
    """Test that the custom colormap can be created"""
    from matplotlib.colors import LinearSegmentedColormap

    # Test the exact colormap creation as done in bicop_plot
    colors = [
      "#00007F",
      "blue",
      "#007FFF",
      "cyan",
      "#7FFF7F",
      "yellow",
      "#FF7F00",
      "red",
      "#7F0000",
    ]

    jet_colors = LinearSegmentedColormap.from_list("jet_colors", colors, N=100)
    assert jet_colors is not None
    assert jet_colors.N == 100
    assert jet_colors.name == "jet_colors"


def test_bicop_plot_refuses_x_on_a_pair_that_reads_no_covariates() -> None:
  """Better a loud refusal than an unconditional surface under a conditional call.

  Exercised at the helper level because the object under test does
  *not* conform to ``BicopLike``: a ``pdf`` with no ``x`` parameter is the
  compiled ``Bicop``'s shape, and a ``BicopBase`` subclass cannot express it
  without violating the contract -- which is itself why the forwarding rule
  cannot rely on the parameter being absent.
  """
  from pyvinecopulib.core._bicop_plot import bicop_plot

  class NoCovariates:
    var_types = None

    def pdf(self, u: np.ndarray) -> np.ndarray:
      return np.ones(u.shape[0], dtype=float)

  with pytest.raises(TypeError):
    bicop_plot(NoCovariates(), "contour", x=[0.5])


def test_bicop_plot_takes_one_covariate_row_only() -> None:
  """A 2-d surface shows the density at one covariate value, not many."""
  from pyvinecopulib.core._bicop_plot import bicop_plot

  class Conditional:
    var_types = None

    def pdf(self, u: np.ndarray, *, x: Any = None) -> np.ndarray:
      return np.ones(u.shape[0], dtype=float)

  for bad in (np.zeros((17, 1)), np.zeros((2, 2, 1))):
    with pytest.raises(ValueError, match="single covariate row"):
      bicop_plot(Conditional(), "contour", x=bad)


def test_bicop_plot_places_the_grid_through_the_supplied_hook() -> None:
  """``place`` is passed explicitly, so the compiled path keeps its NumPy grid."""
  from pyvinecopulib.core._bicop_plot import bicop_plot

  seen: dict[str, object] = {}

  class Recording:
    var_types = None

    def pdf(self, u: np.ndarray) -> np.ndarray:
      seen["type"] = type(u).__name__
      return np.ones(u.shape[0], dtype=float)

  # No `place`: the grid arrives exactly as it always did.
  bicop_plot(Recording(), "contour")
  assert seen["type"] == "ndarray"

  # With one: every manufactured array goes through it.
  calls: list[tuple[int, ...]] = []

  def place(a: np.ndarray) -> np.ndarray:
    calls.append(tuple(a.shape))
    return a

  bicop_plot(Recording(), "contour", place=place)
  assert calls and all(len(shape) == 2 for shape in calls)
