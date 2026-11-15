# Preparing a tokenless PyPI release

Publishing a GitHub release starts validation and requests approval to upload
the distributions to PyPI. A reviewer must approve the `pypi` environment
before the upload runs.

## Track a release

[Open a release tracker](https://github.com/lalmei/hypra/issues/new?template=release-tracker.md)
for each version. Use it to record the exact tag, validation run, release draft,
remaining blockers, and links to the published distributions.

After merging the version bump, use [Draft release notes](release-drafter.md)
to create the GitHub draft and generate its categorized PR list.

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
job. If the sole maintainer publishes the GitHub release or dispatches the run,
do not enable a prevention-of-self-review rule unless another reviewer can
approve deployments.

The environment's existence and required reviewers must be verified before a
release; creating the workflow alone does not configure repository settings.

## What runs automatically

`.github/workflows/release.yml` pins third-party Actions to full commit hashes.
Version-tag pushes, published GitHub releases and manual runs all validate the
distributions:

1. Checks that a tag, when selected, exactly equals `v` plus `project.version`.
2. Installs locked Python dependencies and runs tests excluding the explicitly
   external `real_backend` marker.
3. Builds the wheel and source distribution with `uv build`.
4. Runs `twine check --strict` and uploads distributions as a review artifact.

Publishing a stable release or prerelease on GitHub starts the PyPI publication
path. Tag pushes and manual runs with `publish=false` stop after validation.
To request publication manually, dispatch the workflow on the version tag with
`publish=true`. Both publication paths require successful validation and `pypi`
environment approval.

GitHub runs the workflow from the release tag. The frozen `v0.1.0a0` tag has the
earlier workflow, so use manual dispatch with `publish=true` for that version.
For subsequent releases, tag a commit that includes the updated workflow.

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
and run tests and local metadata validation. After pushing the version tag,
inspect its validation artifact before publishing the GitHub release or
requesting a manual production upload.
A successful build is not evidence of publication: verify the resulting PyPI
project/version and downloadable artifacts afterward. PyPI versions cannot be
overwritten, so choose the release number deliberately.
