import re


def process_changelog(text, repo_url):
  # # Match individual PR numbers (#123)
  # pr_pattern = r"#(\d+)"

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
### BREAKING API CHANGES

* The `mst_algorithm` option to `FitControlsVinecop` has been renamed to `tree_algorithm` to
  allow for alternative spanning tree algorithms (#637).
* `tree_algorithm`'s default value is now `"mst_prim"` instead of `"prim"`, and `"mst_kruskal"`
  replaces `"kruskal"` (#637).
* The CMake option `VINECOPULIB_BUILD_SHARED_LIBS` has been changed to `VINECOPULIB_PRECOMPILED`
  to better reflect its purpose (#641).

### NEW FEATURES

* Allow for random spanning trees as alternatives to the MST-based structure selection using
  `tree_algorithm` in `FitControlsVinecop` with `"random_weighted"` or `"random_unweighted"`
  (#637).

### BUG FIXES

* Decouple edge insertion from criterion computation in `VinecopSelector` to fix randomness
  issues in structure selection when using multiple threads (#640)
"""

# Repository URL
repository_url = "https://github.com/vinecopulib/vinecopulib"

# Process the changelog
processed_changelog = process_changelog(changelog_text, repository_url)
print(processed_changelog)
