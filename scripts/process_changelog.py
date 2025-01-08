import re


def process_changelog(text, repo_url):
  # Match individual PR numbers (#123)
  pr_pattern = r"#(\d+)"

  # Match commit hashes (abcdef0) or similar
  commit_pattern = r"\b[0-9a-f]{7,40}\b"

  # Replace all PR references with Markdown links
  def replace_prs(match):
    prs = match.group(0)
    return (
      "("
      + ", ".join(
        f"[#{pr}](https://github.com/vinecopulib/vinecopulib/pull/{pr})"
        for pr in re.findall(r"\d+", prs)
      )
      + ")"
    )

  # Replace individual commit hashes with Markdown links
  def replace_commit(match):
    commit_hash = match.group(0)
    return f"([{commit_hash}]({repo_url}/commit/{commit_hash}))"

  # Replace PR patterns in parentheses
  text = re.sub(r"\(#(\d+(, ?#\d+)*)\)", replace_prs, text)

  # Replace individual commit hashes
  text = re.sub(commit_pattern, replace_commit, text)

  return text


# Example changelog
changelog_text = """
#### New features

* Use analytical derivatives in discrete pdf/hfuncs (#572)
* Allow for alternative for `"prim"` vs `"kruskal"` in MST-based model selection (#577)
* Improve the dependencies install script to use it in other projects (#576)
* Add tawn copula (#579)
* Improve doc (#580, #585, #607)
* Allow for the discrete Rosenblatt transform (#581)
* Add `Vinecop::fit()` (#584)
* Improve `Bicop::str()` (#588) and `Vinecop::str()` (#589)
* Properly handle discrete variables for the TLL family (#597)
* Weighted pseudo-observations (#602)
* Cross-platform random numbers and add seeds options to `to_pseudo_obs` (#603)
* Improve performance by
    * aligning with the `R` defaults (e.g., `BOOST_NO_AUTO_PTR`, `BOOST_ALLOW_DEPRECATED_HEADERS`, `BOOST_MATH_PROMOTE_DOUBLE_POLICY=false`, `std::string nonparametric_method = "constant"` for the TLL instead of `"quadratic"`, `-O3 -march=native` compiler flags) and add benchmarking example (#592, #611, #613),
    * using `Eigen` element-wise operations instead of `boost` whenever possible (#598, #612),
    * using binary search in the TLL for `get_indices` (#613).

#### Bug fixes

* Improve stability in BB7 PDF (#573)
* Revamped CI/CD pipeline, tests discoverable by CTest, boost version on windows (66cf8b0)
* Fix ASAN issues (#583)
* Fix interface includes and other CMake issue (#586, #599, #601, #608), by @jschueller
"""

# Repository URL
repository_url = "https://github.com/vinecopulib/vinecopulib"

# Process the changelog
processed_changelog = process_changelog(changelog_text, repository_url)
print(processed_changelog)
