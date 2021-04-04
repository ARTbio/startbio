#!/usr/bin/env bash
pip install mkdocs mkdocs-material pygments pymdown-extensions
mkdir mkdocs_build
cd mkdocs_build
# Initialize gh-pages checkout
DATE=`date`
git clone https://github.com/ARTbio/startbio.git
cd startbio
shuf docs/Covid-19/youtube.md -o shuffle.txt
mv shuffle.txt docs/Covid-19/youtube.md
git config credential.helper "store --file=.git/credentials"
echo "https://${GH_TOKEN}:@github.com" > .git/credentials
mkdocs gh-deploy --clean -m "gh-deployed by travis $DATE"
