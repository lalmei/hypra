---
name: Release tracker
about: Track release readiness, approval, and publication for one Hypra version.
title: "Release tracker for v"
---

### Release target

- Version and tag:
- Scope and known limitations:
- GitHub release draft:
- Tag validation run:

### Before publication

- [ ] Confirm the release scope, package version, and prerelease status. If the tag already exists, verify its commit and keep the tag fixed.
- [ ] Verify that validation passed for the exact release tag.
- [ ] Confirm the project's license and distribution policy for this release. Inspect the wheel and source distribution for restricted HypExp/HPL sources or tables.
- [ ] Review the release notes, distribution files, and checksums in the GitHub draft.
- [ ] Verify that the PyPI Trusted Publisher matches this repository, `release.yml`, and the `pypi` environment.
- [ ] Verify that the GitHub `pypi` environment has a required reviewer.

### Publication

- [ ] Dispatch the release workflow on the version tag with `publish=true`. Record the run link and obtain the `pypi` environment approval.
- [ ] Verify the version and both distributions on PyPI. Record the PyPI link.
- [ ] Publish the GitHub release with the intended prerelease status. Record the public release link.

### Blockers and evidence

<!-- Link decisions, validation results, and publication approvals here. -->
