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
### NEW FEATURES

* add `allow_rotation` option to `FitControlsBicop` and `FitControlsVinecop`
  to allow for the rotation of the pair copulas (#628).

* add a `FitControlsConfig` struct to create flexible and yet safe constructors
  for `FitControlsBicop` and `FitControlsVinecop` (#629).

### BUG FIXES

* restrict parameter range for fitting Tawn copulas; fix handling of their 
  shape/argument order (#620).

* compute and save loglik/nobs in `Vinecop::fit()` (#623)

* disable unwanted compiler output related to BOOST_CONCEPT checks (#624)
"""

# Repository URL
repository_url = "https://github.com/vinecopulib/vinecopulib"

# Process the changelog
processed_changelog = process_changelog(changelog_text, repository_url)
print(processed_changelog)
