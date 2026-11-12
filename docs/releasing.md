# Preparing a PyPI release

The release workflow validates distributions first and keeps production uploads
behind a separate manual action and protected environment review. Creating the
workflow does not upload a package.

## Current readiness audit

Read-only checks on 12 November 2026 found:

- The public PyPI JSON endpoint for `hypra` returned HTTP 404. This suggests no
  current published project, but does not reserve the name or prove ownership.
- GitHub reported no repository secret names and no deployment environments
  for `lalmei/hypra`. Only names/settings were queried, never secret values.
- No authenticated PyPI account or configured Trusted Publisher was verified.
  A repository owner and a PyPI project owner are not automatically the same.

These are setup blockers to resolve before a production release. The package's project licence/distribution policy must also
be settled before release; no restricted HypExp/HPL source or tables may be
added to the wheel or source distribution.

## What the workflow does

`.github/workflows/release.yml` pins third-party Actions to full commit hashes.
A version-tag push or a manual run with `publish=false` only:

1. Checks that a tag, when selected, exactly equals `v` plus `project.version`.
2. Installs locked Python dependencies and runs tests excluding the explicitly
   external `real_backend` marker.
3. Builds the wheel and source distribution with `uv build`.
4. Runs `twine check --strict` and uploads the distributions as a review artifact.

It never publishes automatically on a tag push. A production upload requires
**all** of the following:

- A separate manual dispatch with `publish=true`, selecting a version tag.
- Successful tests, build and metadata checks for that exact run.
- An already existing GitHub environment named `pypi`, with required reviewers.
  A preflight query fails before the publishing job if this protection is absent;
  the workflow does not create an unprotected environment as a fallback.
- Approval of the environment's deployment review.
- An explicitly configured secret named `PYPI_API_TOKEN` in the protected
  environment or repository, scoped to the authorized PyPI project where possible.
- Confirmed PyPI ownership/name availability and a maintainer-approved release.

The workflow has no repository-write or OIDC permission. Because no existing
Trusted Publisher was verified, this draft uses an API-token path and disables
OIDC-based attestations. Do not paste tokens into source, workflow inputs,
issues or logs. If Trusted Publishing is configured later, review a separate
change to use its exact repository/workflow/environment identity, remove the
password path and grant only the publishing job `id-token: write`.

## Release procedure after approval and setup

First confirm that every file intended for the distribution is publishable,
choose the version in `pyproject.toml`, update the lock if metadata changes,
and run local tests. For the current prerelease version, the corresponding tag
would be `v0.1.0a0`; this page does not create it.

After an authorized tag push, inspect the validation artifact and its metadata.
When the maintainer has approved the release, select that
same tag in the workflow's manual Run form, set `publish=true`, and review the
protected-environment request. A normal branch selection cannot publish.

PyPI versions cannot be overwritten. The workflow does not skip an existing
version or silently select a new version. A collision, missing credential,
missing protection, unsuccessful check or failed upload must be resolved
explicitly before retrying. `twine check` validates metadata, not account
ownership or legal permission.

Local metadata check, without uploading:

```sh
uv build --out-dir dist
uvx --from twine==7.0.0 twine check --strict dist/*
```

This workflow does not enable GitHub Pages, deploy the documentation, modify
repository settings, create tokens or set up Trusted Publishing.
