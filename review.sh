#!/bin/bash

set -ev

echo running regressions
sleep 2
for j in regression/*.json; do
    echo $j
    cargo run --release $j
done

echo running todo
sleep 2
for j in todo/*.json; do
    echo $j
    cargo run --release $j
done
