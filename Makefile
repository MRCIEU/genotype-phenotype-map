.PHONY: format lint lint-check lint-fix lint-summary test help

help:
	@echo "Available commands:"
	@echo "  make format        - Format all R files using formatR"
	@echo "  make lint          - Run lintr on all R files and show results"
	@echo "  make lint-check    - Run lintr and exit with error code if issues found"
	@echo "  make lint-fix      - Auto-fix linting issues (formatting, line length, etc.)"
	@echo "  make lint-summary  - Show summary of linting issues by type"
	@echo "  make test          - Run testthat tests"

format:
	@echo "Formatting R files..."
	@Rscript scripts/lint_fix.R
	@echo "Formatting complete!"

lint:
	@echo "Linting R files..."
	@Rscript scripts/lint.R

lint-check:
	@echo "Running lint check..."
	@Rscript scripts/lint_check.R

lint-fix:
	@echo "Auto-fixing linting issues..."
	@Rscript scripts/lint_fix.R
	@echo "Auto-fix complete!"

lint-summary:
	@echo "Analyzing linting issues..."
	@Rscript scripts/lint_summary.R

test:
	@echo "Running tests..."
	@Rscript tests/testthat/run_tests.R
