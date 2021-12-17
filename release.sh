#!/bin/bash
set -euo pipefail
V=$1
echo $V

# broken for some reason
# target=aarch64-unknown-linux-musl

mkdir -p dists

for target in x86_64-apple-darwin aarch64-apple-darwin x86_64-unknown-linux-musl; do
    echo $target
    cross build --release --target $target
    tar -czvf ./dists/rustybam_v${V}-${target}.tar.gz \
        -C ./target/$target/release/ \
        rustybam rb
done

gh release create \
    -t "Release v${V}" \
    -n "v${V}" \
    "v${V}" \
    ./dists/rustybam_v${V}-*tar.gz
