# Releasing Opal

Opal uses **tag-driven releases**: the git tag is the single source of truth
for a release, and cutting one is a single repeatable action. This eliminates
the version drift that previously existed across `CMakeLists.txt`,
`include/opal/version.hpp` and `python/setup.py` (issue #9).

## Versioning policy

Opal follows [Semantic Versioning](https://semver.org/) — `MAJOR.MINOR.PATCH`:

- **MAJOR** — incompatible public API changes (C++ headers, CLI flags, or the
  Python API).
- **MINOR** — new functionality added in a backward-compatible way (a new
  model, engine, instrument or binding).
- **PATCH** — backward-compatible bug fixes and numerical-accuracy corrections.

Pre-`1.0.0` the API is still settling, so minor versions may carry small
breaking changes; these are called out in the release notes.

## Single source of truth

The version lives in exactly one place:

```c
// include/opal/version.hpp
#define OPAL_VERSION_STRING "0.2.0"
```

Everything else derives from it, so the files can never disagree:

| Consumer | How it gets the version |
|---|---|
| C++ library / CLI | `OPAL_VERSION_STRING` directly (header-only) |
| CMake (`project(opal VERSION …)`) | parses the header at configure time |
| Python wheel / sdist (`opal-pricing`) | `setup.py` parses the bundled header |

The pushed git tag is authoritative for a *release*: the CI release gate fails
if the tag and the header disagree (see below), so a release can never be
published against a stale version.

## Cutting a release

### Option A — the helper script (recommended)

```sh
scripts/release.sh 0.3.0
```

This bumps `include/opal/version.hpp`, commits the bump (`release: v0.3.0`), and
creates an annotated tag `v0.3.0`. It does **not** push. Review, then:

```sh
git push origin main && git push origin v0.3.0
```

### Option B — by hand

1. Edit `include/opal/version.hpp` so `OPAL_VERSION_STRING` (and the
   `MAJOR`/`MINOR`/`PATCH` macros) read the new version.
2. Commit: `git commit -am "release: v0.3.0"`.
3. Tag: `git tag -a v0.3.0 -m "opal v0.3.0"`.
4. Push the branch and the tag: `git push origin main && git push origin v0.3.0`.

## What happens on a tag push

Pushing a `v*` tag triggers `.github/workflows/release.yml`:

1. **version-check** — fails the whole release if `v<tag>` does not match
   `OPAL_VERSION_STRING` in `include/opal/version.hpp`.
2. **build-cli** — builds the `opal` CLI on Linux, macOS and Windows.
3. **wheels** — reuses the `wheels.yml` reusable workflow to build the
   `cibuildwheel` matrix of wheels plus an sdist (not rebuilt elsewhere).
4. **publish** — creates a **GitHub Release** for the tag with auto-generated
   notes (from the merged PRs since the previous tag) and attaches the CLI
   binaries, wheels and sdist.

Publishing to PyPI is intentionally **out of scope** — releases are
GitHub-only. The wheels are attached to the release for direct download.

## If the release gate fails

A failed `version-check` means the tag and the header disagree. Either you
tagged without bumping the header, or you bumped to a different version. Fix it
by deleting the bad tag and re-cutting:

```sh
git tag -d v0.3.0 && git push origin :refs/tags/v0.3.0   # remove remote tag
scripts/release.sh 0.3.0                                  # bump header + retag
git push origin main && git push origin v0.3.0
```
