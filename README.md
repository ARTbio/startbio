[![Release Documentation](https://img.shields.io/badge/release-documentation-blue.svg)](https://artbio.github.io/startbio/)
[![Docs](https://github.com/ARTbio/startbio/actions/workflows/docs.yml/badge.svg)](https://github.com/ARTbio/startbio/actions/workflows/docs.yml)
![GitHub issues](https://img.shields.io/github/issues/ARTbio/startbio?color=purple)
![GitHub pull requests](https://img.shields.io/github/issues-pr-raw/ARTbio/startbio?color=blue)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/artbio/startbio?color=darkgreen)

# startbio
#### Training Material of the STARTbio program   [![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/) ![GitHub watchers](https://img.shields.io/github/watchers/artbio/startbio?style=social) ![Twitter URL](https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Ftwitter.com%2FARTbio_IBPS)


The documentation is built with MkDocs and Material for MkDocs, and published at [https://artbio.github.io/startbio/](https://artbio.github.io/startbio/)

## Build the site locally

The build environment is pinned: Python 3.12 and the exact package versions in `requirements.txt`
(top-level dependencies are listed in `requirements.in`).

```bash
micromamba create -f environment.yml   # once
micromamba activate startbio
mkdocs serve --strict                  # live preview at http://127.0.0.1:8000/startbio/
mkdocs build --strict                  # the same check that CI runs on every pull request
```

Without conda/micromamba: `python3.12 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt`.

## Publication

The workflow `.github/workflows/docs.yml` runs `mkdocs build --strict` on every pull request to `master`.
The build fails on any warning: broken internal link or anchor, or page missing from the `nav`.
When a pull request is merged (push to `master`), the same build runs again and `mkdocs gh-deploy`
publishes the site to the `gh-pages` branch, which GitHub Pages serves.
