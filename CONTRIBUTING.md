# Contributing to paftacular

Bug reports, feature proposals, documentation improvements, and code contributions are welcome.

## Issues and support

Use the [GitHub issue tracker](https://github.com/tacular-omics/paftacular/issues) for reproducible bug reports, feature requests, and usage questions. Include the paftacular and Python versions, a minimal example, the expected behavior, and the observed behavior. For suspected security vulnerabilities, follow [SECURITY.md](SECURITY.md) instead of opening a public issue.

## Development setup

Install [uv](https://docs.astral.sh/uv/) and [just](https://just.systems/), then run:

```bash
git clone https://github.com/tacular-omics/paftacular.git
cd paftacular
just install
just pre-release
```

`just install` is `uv sync --all-extras`, so it installs the three optional extras
(`smiles`, `peptacular` and `mcp`) alongside the dev group. The tests for those parts
call `pytest.importorskip` in two places, and the two behave differently. At module
scope a whole test file is never collected, which lowers the collected total and
produces no skip. Inside a single test the test is collected and then skipped, which
leaves the total unchanged and moves the skip count. A missing extra therefore shows
up in one number or the other, so read the collected total together with the passed
and skipped counts. Neither on its own tells you whether the parts you changed ran.
Some tests are guarded by neither and fail with paftacular's own `ImportError` naming
the extra, so a partial install can also show failures that are not about your change.

Install the commit hooks once per clone with `uvx pre-commit install`. On every
commit they run `ruff check`, `ruff format --check` and a few file checks
(`uvx pre-commit run --all-files` runs them on everything). Type checks and
tests are not in the hooks. Run `just check` for those.

Changes should include tests for new behavior and user-facing documentation when applicable. Keep pull requests focused, explain the motivation, and describe how the change was validated.

By participating, you agree to follow the [Code of Conduct](CODE_OF_CONDUCT.md).
