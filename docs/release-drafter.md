# Maintaining draft release notes

## Create or refresh a release draft

Run **Draft release notes** from GitHub Actions after merging the package version
bump. Enter the previous version tag. The workflow uses the selected `main`
commit by default. To refresh a draft pinned to an earlier commit, enter that
commit's full SHA in `target`.

The workflow reads the version from the target commit and lists PRs merged into
`main` between the previous tag and that commit. Each entry has the original PR
title, linked number and `@author`. Labels take priority over conventional PR
title prefixes such as `feat:` or `fix:`. Unclassified entries go under
Maintenance. Empty categories are omitted.

Write a short feature summary above the generated section. Later runs preserve
text outside the markers and refresh the PR list and comparison link. The
workflow saves a GitHub draft. Publication and PyPI approval follow the
[release process](releasing.md).

To print the notes without saving them, run this locally with GitHub CLI auth:

```sh
uv run --no-sync python scripts/draft_release.py --repo lalmei/hypra \
  --previous-tag v0.1.0a0 --target a10e1b02a46b30b09dded44a4cd2adc748d9238b
```

## Preview the frozen alpha notes

The separate **Preview draft release notes** workflow is for `v0.1.0a0`. The default
`apply=false` generates a downloadable Markdown proposal with a read-only token.
Selecting `apply=true` additionally updates **only** the reserved generated-notes
region in the existing unpublished `v0.1.0a0` draft. There is no push/merge trigger,
automatic version bump, release publication, asset upload/deletion or tag write.

## Protecting the frozen alpha

The `v0.1.0a0` tag must resolve to
`e1c13fc82a9077720a6b27350a5729993e772a83`. Release Drafter's comparison target
is pinned to that commit, not the moving tip of `main`. The optional apply job
also checks each proposed PR independently: its merged commit must be an
ancestor of, or identical to, the immutable tag commit. Entries merged after
the tag, or on a divergent history, are excluded. Malformed entries, unmerged
PRs and failed API comparisons stop the update rather than guessing.

The workflow snapshots the draft ID, timestamp and body hash before generating
notes, then rechecks them, its draft status and the immutable tag before applying.
If someone edited or published the draft during generation, it fails and asks
for a new run. Do not manually edit the draft while an apply run is active:
GitHub's release update API does not provide an atomic compare-and-swap, so a
final check cannot eliminate a simultaneous edit in the last request interval.

Exactly one ordered marker pair must already exist:

```text title="Reserved generated-notes region"
<!-- release-drafter:start -->
...generated pull-request notes...
<!-- release-drafter:end -->
```

The wrapper preserves all bytes outside that region, including the curated
overview, checksums, asset links and limitations. Its only release update field
is `body`; it never changes draft/publication state, the release target or assets.
Only the optional apply job has `contents: write`. Generation always runs with
read-only permissions and Release Drafter's own `dry-run: true`.

## Why a preservation wrapper is necessary

A normal Release Drafter run replaces the release body. Its `$BODY` variable
means a **pull request body**, not the draft's curated notes. We therefore use
its generated output as data and perform the narrower, ancestry-checked edit
ourselves. Do not replace this setup with the action's normal write mode.

An empty PR list does not mean the release has no changes: direct commits do
not necessarily have associated PRs. The generated region says so and leaves
the curated release overview intact. PR titles are escaped as Markdown data;
they are never evaluated as shell or JavaScript source.

## PEP 440 alpha version

The package version is `0.1.0a0` and actual tag is `v0.1.0a0`. Release Drafter's
internal SemVer input is `0.1.0-alpha.0`; explicit `tag` and `name` overrides
preserve the PEP 440 spelling. Future versions use **Draft release notes**
and read the package version from their own target commits.

## Verification and boundaries

The official pinned [v7.10.0 action](https://github.com/release-drafter/release-drafter/blob/978bb0b30c5ee8cdd332ced253758157150d43f0/action.yml)
declares `dry-run` and the generated `body` output. See its
[configuration documentation](https://github.com/release-drafter/release-drafter/blob/v7.10.0/README.md).
Run `node scripts/test_release_drafter.cjs` for offline safety regressions.
Local checks validate YAML and exercise the actual inline wrapper with mocked GitHub data:
included/post-tag PRs, unchanged curated notes, changed drafts, published
releases and malformed entries. They do not claim a successful live GitHub run.

No personal access token is required. This workflow does not interact with
PyPI or bypass a pending publication approval. Installing or preparing it does
not itself dispatch it or change the current draft.
