default: lint format check test

# Install dependencies
install:
    uv sync

# Run linting checks
lint:
    uv run ruff check src

# Format code
format:
	uv run ruff check --select I --fix src tests
	uv run ruff format src tests

# Run type checking
check:
    uv run ty check src

# Run tests
test:
    uv run pytest tests

test-cov:
    uv run pytest tests --cov=src --cov-branch --cov-report=term-missing --cov-report=html --cov-report=xml

codecov-tests:
    uv run pytest tests --cov=src --junitxml=junit.xml -o junit_family=legacy