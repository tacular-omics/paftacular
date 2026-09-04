default: check

install:
    uv sync --all-extras

install-all: install

install-prod:
    uv sync --no-dev --frozen

lint:
    uv run ruff check src tests benchmarks

format:
    uv run ruff check --select I --fix src tests benchmarks
    uv run ruff format src tests benchmarks

format-check:
    uv run ruff format --check src tests benchmarks

ty:
    uv run ty check src

check: lint format-check ty test

test:
    uv run pytest tests

test-cov:
    uv run pytest tests --cov=src --cov-branch --cov-report=term-missing --cov-report=html --cov-report=xml

codecov-tests:
    uv run pytest tests --cov=src --cov-branch --junitxml=junit.xml

docs:
    uv run sphinx-build -W -b html docs docs/_build/html

docs-test:
    uv run sphinx-build -W -b doctest docs docs/_build/doctest

docs-clean:
    rm -rf docs/_build

docs-open: docs
    uv run python -m webbrowser "file://{{justfile_directory()}}/docs/_build/html/index.html"

benchmark:
    uv run python benchmarks/parse.py

build:
    uv build

pre-release: format check docs docs-test build
