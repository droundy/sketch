#!/bin/bash

set -ev

mkdir -p release

rustup target add x86_64-unknown-linux-musl
cargo bundle --release --format deb --target x86_64-unknown-linux-musl
cp target/x86_64-unknown-linux-musl/release/bundle/deb/color-splotch_*_amd64.deb release/

rustup target add x86_64-pc-windows-gnu
cargo build --release --target x86_64-pc-windows-gnu
# apt install gcc-mingw-w64
rm -f release/color-splotch.exe
cp target/x86_64-pc-windows-gnu/release/color-splotch.exe release/

rustup target add x86_64-apple-darwin aarch64-apple-darwin
cargo bundle --release --format osx --target x86_64-apple-darwin
cargo bundle --release --format osx --target aarch64-apple-darwin
