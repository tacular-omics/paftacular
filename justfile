default: lint format check test

# Install all dependencies (dev + all extras)
install:
    uv sync --all-extras

install-all: install
    

# Install minimal dependencies (no dev, no extras)
install-prod:
    uv sync --no-dev --frozen

# Run linting checks
lint:
    uv run ruff check src

# Format code
format:
	uv run ruff check --select I --fix src tests
	uv run ruff format src tests

# Run ty type checker
ty:
    uv run ty check src


# Run type checking
check:
    just lint
    just ty
    just test

# Run tests
test:
    uv run pytest tests

test-cov:
    uv run pytest tests --cov=src --cov-branch --cov-report=term-missing --cov-report=html --cov-report=xml

codecov-tests:
    uv run pytest tests --cov=src --junitxml=junit.xml -o junit_family=legacy

# Build documentation
docs:
    cd docs && uv run sphinx-build -b html . _build/html

docs-test:
    cd docs && uv run sphinx-build -b doctest . _build/doctest


# Clean documentation build
docs-clean:
    rm -rf docs/_build

# Build and open documentation
docs-open:
    just docs
    python -c "import webbrowser; webbrowser.open('file://{{justfile_directory()}}/docs/_build/html/index.html')"
