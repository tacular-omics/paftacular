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
(`smiles`, `peptacular` and `mcp`, with `all` as their union) alongside the dev group.
The tests for those parts call `pytest.importorskip` in two places, and the two behave
differently. At module scope the file's tests are never collected and the run reports a
single skip standing for the whole file, so the collected total falls by however many
tests that file holds and the skip count rises by one. Inside a single test the test is
collected and then skipped, so the collected total is unchanged and the skip count rises
by one. Neither count on its own tells the two apart. A subtraction does: add the
passed, failed, skipped, xfailed, xpassed and error counts, leave out any subtest count,
and subtract the collected total. The difference is the number of files skipped at
collection, and it is zero when every extra is installed. pytest prints that number on
the collection line too, as `collected N items / M skipped`, but `-q` drops the line and
the final summary never carries it. Some tests are guarded by neither and fail with
paftacular's own `ImportError` naming the extra, so a partial install can also show
failures that are not about your change.

Install the commit hooks once per clone with `uvx pre-commit install`. On every
commit they run `ruff check`, `ruff format --check` and a few file checks
(`uvx pre-commit run --all-files` runs them on everything). Type checks and
tests are not in the hooks. Run `just check` for those.

Changes should include tests for new behavior and user-facing documentation when applicable. Keep pull requests focused, explain the motivation, and describe how the change was validated. AI-assisted contributions are welcome; you are responsible for reviewing and testing what you submit.

By participating, you agree to follow the [Code of Conduct](CODE_OF_CONDUCT.md).
