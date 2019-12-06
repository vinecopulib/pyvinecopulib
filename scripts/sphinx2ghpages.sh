#!/bin/bash

doc_dictory=docs/_build
branch_name=gh-pages

# get a clean dev branch
git checkout dev 
git pull origin dev 
git clean -df
git checkout -- .
git fetch --all

# install the library
pip install .

# build html docs from sphinx files
cd docs
python3 serve_sphinx.py --generates=True
cd ..

# create or use orphaned gh-pages branch
if [ $(git branch --list "$branch_name") ]
then
	git stash
	git checkout $branch_name
	git pull origin $branch_name
	#git stash apply
	git checkout stash -- . # force git stash to overwrite added files
else
	git checkout --orphan "$branch_name"
fi

if [ -d "$doc_directory" ]
then
	mv docs/_build/* . && rm -rf docs lib src examples tests .github MANIFEST.in requirements.txt setup.py
  mv _build/* . && rm -rf _build
  touch .nojekyll
	git add -A
	git commit -m "new pages version $(date)"
	git push origin gh-pages
	# github.com recognizes gh-pages branch and create pages
	# url scheme https//:[github-handle].github.io/[repository]
else
	echo "directory $doc_directory does not exists"
fi
