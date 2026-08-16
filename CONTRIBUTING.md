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

Changes should include tests for new behavior and user-facing documentation when applicable. Keep pull requests focused, explain the motivation, and describe how the change was validated.

By participating, you agree to follow the [Code of Conduct](CODE_OF_CONDUCT.md).
