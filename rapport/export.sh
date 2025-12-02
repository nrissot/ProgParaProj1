#!/bin/bash
pandoc rapport.md -t html -s --katex=https://cdn.jsdelivr.net/npm/katex@0.16.25/dist/ --metadata title="ProgPara - MPI Projet 1" -o o.html
