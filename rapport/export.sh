#!/bin/bash
pandoc rapport.md -t html -s --katex=https://cdn.jsdelivr.net/npm/katex@0.16.25/dist/ -o o.html
