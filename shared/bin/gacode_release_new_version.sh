#!/bin/bash
set -e
if [ $# -ne 1 ] ;
then
  echo "Usage: $0 <new_version>"
  exit 1
fi
cd $GACODE_ROOT
git checkout master
git pull
git checkout stable
git pull
git merge master
python shared/bin/gacode_regression.py -clean
git tag -a $1
git push
git push --tags
