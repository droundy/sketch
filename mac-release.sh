#!/bin/bash

set -ev

mkdir -p release

rustup target add x86_64-apple-darwin aarch64-apple-darwin
cargo bundle --release --format osx --target x86_64-apple-darwin
cargo bundle --release --format osx --target aarch64-apple-darwin

APP_NAME='Color Splotch'

rm -rf "$APP_NAME"
mkdir "$APP_NAME"
cp -r target/x86_64-apple-darwin/release/bundle/osx/ColorSplotch.app "$APP_NAME/"
ln -s /Applications "$APP_NAME/Applications"
rm -rf "$APP_NAME/.Trashes"

hdiutil create "release/$APP_NAME-x86-64.dmg" -srcfolder "$APP_NAME" -ov

rm -rf "$APP_NAME"
mkdir "$APP_NAME"
cp -r target/aarch64-apple-darwin/release/bundle/osx/ColorSplotch.app "$APP_NAME/"
ln -s /Applications "$APP_NAME/Applications"
rm -rf $APP_NAME/.Trashes

hdiutil create "release/$APP_NAME-aarch64.dmg" -srcfolder "$APP_NAME" -ov
