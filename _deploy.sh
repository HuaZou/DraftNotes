#!/bin/sh

set -e

git add -A
git commit -m "Update the book" || true
#git push -q origin gh-pages
git push origin main

