## Why

<!-- The problem or need. Link the issue or the upstream pull request. -->

## What changed

<!-- One short paragraph per change. Say what a reviewer should look at. -->

## Testing

<!-- What you ran, and what a reviewer should run. Numbers that moved, and
     why they are the right numbers now. -->

---

- [ ] Tests added or extended
- [ ] `CHANGELOG.md` updated under the `(unreleased)` heading
- [ ] Public-API changes reach the module docstring and the matching notebook
- [ ] `make check && make test && make docs` pass locally

<!-- Submodule bumps additionally run the numerics gate:
     tests/test_torch_bicop.py, tests/test_torch_vinecop.py and
     tests/test_structure_selection.py. Regenerate expected values rather than
     widening tolerances. See AGENTS.md. -->
