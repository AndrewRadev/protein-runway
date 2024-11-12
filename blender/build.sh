#! /bin/sh

set -e

blender_dir="$(dirname "${BASH_SOURCE[0]}")"

pushd "$blender_dir/extension"
zip -r ../extension.zip wheels/ __init__.py blender_manifest.toml
popd
