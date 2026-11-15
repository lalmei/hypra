UV ?= uv
PYTHON_RUN = $(UV) run --locked
BUMP_VERSION = $(UV) tool run --from bump-my-version==1.5.1 bump-my-version
CURRENT_VERSION = $(BUMP_VERSION) show current_version
PYTEST_ARGS ?=
DOCS_ADDR ?= 127.0.0.1:8000
DIST_DIR ?= dist

.DEFAULT_GOAL := help
.PHONY: help sync test lint check docs docs-serve build version bump-preview
.PHONY: bump-major bump-minor bump-patch bump-pre bump-stage bump-release require-prerelease

help: ## List available targets
	@awk -F ':.*## ' '/^[a-z][a-z-]*:.*## / { printf "%-15s %s\n", $$1, $$2 }' $(MAKEFILE_LIST)

sync: ## Install dependencies from uv.lock
	$(UV) sync --locked

test: ## Run the CI test suite (optional PYTEST_ARGS)
	$(PYTHON_RUN) pytest -q -m 'not real_backend' $(PYTEST_ARGS)

lint: ## Lint the Python package with Ruff
	$(PYTHON_RUN) ruff check src

check: lint test ## Run lint and tests

docs: ## Build and verify the documentation
	$(PYTHON_RUN) properdocs build --strict
	$(PYTHON_RUN) python scripts/check_docs_build.py

docs-serve: ## Serve the documentation locally (optional DOCS_ADDR)
	$(PYTHON_RUN) properdocs serve --dev-addr $(DOCS_ADDR)

build: ## Build the wheel and source distribution (optional DIST_DIR)
	$(UV) build --out-dir "$(DIST_DIR)"

version: ## Show the current package version
	@$(CURRENT_VERSION)

bump-preview: ## Show possible version bumps
	@current="$$($(CURRENT_VERSION))" && \
	if [ "$$current" = "$${current%%[a-z]*}" ]; then \
		for part in major minor patch; do \
			next="$$($(BUMP_VERSION) show --increment "$$part" new_version)" || exit; \
			printf '%s -> %s (%s)\n' "$$current" "$$next" "$$part"; \
		done; \
	else \
		$(BUMP_VERSION) show-bump; \
	fi

bump-major: ## Increment major and start at alpha 0
bump-minor: ## Increment minor and start at alpha 0
bump-patch: ## Increment patch and start at alpha 0
bump-pre: require-prerelease ## Increment the current prerelease number
bump-stage: require-prerelease ## Advance from alpha to beta to rc to stable
bump-major bump-minor bump-patch bump-pre bump-stage:
	$(BUMP_VERSION) bump $(patsubst bump-%,%,$@)
	$(UV) lock

bump-release: ## Promote the current prerelease directly to stable
	@current="$$($(CURRENT_VERSION))" && \
	stable="$${current%%[a-z]*}" && \
	if [ "$$current" = "$$stable" ]; then \
		printf '%s\n' "$$current is already stable"; \
	else \
		$(BUMP_VERSION) bump --new-version "$$stable" && $(UV) lock; \
	fi

require-prerelease:
	@current="$$($(CURRENT_VERSION))" && \
	if [ "$$current" = "$${current%%[a-z]*}" ]; then \
		printf '%s\n' "$$current is stable. Use bump-major, bump-minor, or bump-patch to start a prerelease." >&2; \
		exit 1; \
	fi
