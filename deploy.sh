#!/usr/bin/env bash
pip install mkdocs mkdocs-material
mkdir mkdocs_build
cd mkdocs_build
# Initialize gh-pages checkout
DATE=`date`
git clone https://github.com/ARTbio/startbio.git
cd startbio
git config credential.helper "store --file=.git/credentials"
echo "https://${GH_TOKEN}:@github.com" > .git/credentials
mkdocs gh-deploy --clean -m "gh-deployed by travis $DATE"
