#!/usr/bin/env bash
# Cut an opal release (issue #19).
#
#   scripts/release.sh X.Y.Z
#
# Updates the single source of truth (include/opal/version.hpp), commits the
# bump, and creates an annotated `vX.Y.Z` tag. CMake and the Python package
# derive their versions from that header, so nothing else needs editing. The
# git tag is authoritative for releases; pushing it triggers the Release
# workflow, which re-verifies tag == header before publishing.
#
# This script does not push. Review the commit and tag, then:
#
#   git push origin main && git push origin vX.Y.Z
#
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "usage: $0 X.Y.Z" >&2
  exit 2
fi

version="$1"
if ! printf '%s' "$version" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$'; then
  echo "error: '$version' is not a X.Y.Z semantic version" >&2
  exit 2
fi

root="$(git rev-parse --show-toplevel)"
header="$root/include/opal/version.hpp"

if [ -n "$(git -C "$root" status --porcelain)" ]; then
  echo "error: working tree is not clean; commit or stash first" >&2
  exit 1
fi

tag="v$version"
if git -C "$root" rev-parse "$tag" >/dev/null 2>&1; then
  echo "error: tag $tag already exists" >&2
  exit 1
fi

major="${version%%.*}"
rest="${version#*.}"
minor="${rest%%.*}"
patch="${rest#*.}"

cat > "$header" <<EOF
#pragma once

#define OPAL_VERSION_MAJOR $major
#define OPAL_VERSION_MINOR $minor
#define OPAL_VERSION_PATCH $patch
#define OPAL_VERSION_STRING "$version"
EOF

git -C "$root" add "$header"
git -C "$root" commit -m "release: $tag"
git -C "$root" tag -a "$tag" -m "opal $tag"

echo "Bumped include/opal/version.hpp to $version, committed, and tagged $tag."
echo "Push to publish: git push origin main && git push origin $tag"
