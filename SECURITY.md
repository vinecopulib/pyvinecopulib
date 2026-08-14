# Security policy

## Supported versions

Fixes go into the next release from `main`. Older tags are not patched.

## Reporting a vulnerability

Please report privately through GitHub's
[security advisories](https://github.com/vinecopulib/pyvinecopulib/security/advisories/new)
rather than a public issue.

`pyvinecopulib` is a numerical library. It wraps a C++ extension that parses
model files (JSON and CBOR) and otherwise operates on in-memory numeric arrays,
so the most plausible reports are crashes or memory errors from a malformed
model file, or from an array whose shape violates a documented precondition —
both worth reporting, since a native extension turns those into process
crashes rather than exceptions.

Vulnerabilities in the C++ core belong upstream in
[vinecopulib](https://github.com/vinecopulib/vinecopulib/security/advisories/new);
report them there and mention the `pyvinecopulib` version you saw them through.

Please include `pyvinecopulib.__version__`, your Python version and platform,
whether you installed a wheel or built from source, and the smallest input that
reproduces the problem.
