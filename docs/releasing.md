# Preparing a tokenless PyPI release

The release workflow validates distributions first and keeps production uploads
behind a separate manual action and protected environment review. Creating the
workflow or a Pending Publisher does not upload a package.

## Configure the Pending Publisher

For a new project, open [PyPI's account publishing settings](https://pypi.org/manage/account/publishing/)
and add a **GitHub** Pending Publisher with exactly these fields:

| Field | Value |
| --- | --- |
| PyPI project name | `hypra` |
| GitHub owner | `lalmei` |
| GitHub repository | `hypra` |
| Workflow filename | `release.yml` |
| Environment name | `pypi` |

The workflow filename is **not** its display title or full `.github/workflows/`
path. The environment must match the workflow's `environment.name` exactly.
A Pending Publisher does not reserve the project name. Its first successful
upload creates the project under the authorizing PyPI account. If the project
already exists under that account, configure the same fields under the
project's **Publishing** settings instead. An unrelated existing project
cannot be claimed this way.

Public name checks initially returned HTTP 404 for `hypra`; recheck before
release. An authenticated PyPI owner and successfully saved publisher must be
verified separately. No API token or GitHub token secret is required. Do not
paste account credentials into source, workflow inputs, issues or logs.

## GitHub environment

Create an environment named `pypi` in `lalmei/hypra` with at least one required
reviewer. The workflow checks this protection before scheduling the publishing
job, so it cannot silently create an unprotected environment. Environment
approval is a separate step from manually requesting publication. If the sole
maintainer also dispatches the run, do not enable a prevention-of-self-review
rule unless another reviewer will approve deployments.

The environment's existence and required reviewers must be verified before a
release; creating the workflow alone does not configure repository settings.

## What runs automatically

`.github/workflows/release.yml` pins third-party Actions to full commit hashes.
A version-tag push or a manual run with `publish=false` only:

1. Checks that a tag, when selected, exactly equals `v` plus `project.version`.
2. Installs locked Python dependencies and runs tests excluding the explicitly
   external `real_backend` marker.
3. Builds the wheel and source distribution with `uv build`.
4. Runs `twine check --strict` and uploads distributions as a review artifact.

Production upload requires a separate manual dispatch with `publish=true`,
selecting a version tag, passing validation and approving the `pypi` environment.
Only the publishing job gets `id-token: write`. The official PyPI publish action
exchanges GitHub's short-lived OIDC identity for a project-scoped credential and
creates distribution attestations. There is no stored password/token path.
PyPI will reject the upload unless the saved Trusted Publisher matches this
repository, workflow and environment.

## Release readiness

Before publishing, settle the project license/distribution policy, verify the
PyPI account and publisher, and inspect every distribution file. Restricted
HypExp/HPL source and lookup tables must not enter either archive. Hypra's
capability map states the remaining upstream-parity gaps; choose prerelease
metadata honestly rather than implying full compatibility.

Choose the version in `pyproject.toml`, update the lock if metadata changes,
and run tests and local metadata validation. The current prerelease would use
tag `v0.1.0a0`; this page does not create it. After an authorized tag push,
inspect its validation artifact before requesting a manual production upload.
A successful build is not evidence of publication: verify the resulting PyPI
project/version and downloadable artifacts afterward. PyPI versions cannot be
overwritten, so choose the release number deliberately.
