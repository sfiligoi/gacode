#!/bin/bash
set -e
if [ $# -ne 2 ] ;
then
  echo "Usage: $0 <new_version> <from_branch>"
  echo ""
  echo "Merge <from_branch> into stable to make <new_version>"
  echo "Last tag:" `git tag | tail -1`
  exit 1
fi
cd $GACODE_ROOT
git checkout $2
git pull
git checkout stable
git pull
git merge $2
python shared/bin/gacode_regression.py -clean
git tag -a $1
git push
git push --tags
git checkout master
git pull
#git merge stable
#git push
