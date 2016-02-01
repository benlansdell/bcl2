#!/bin/sh
echo "=======Deploy to espresso========"
git push espresso dev
echo "========Deploy to mocha========"
git push mocha dev
echo "========Deploy to galao========"
git push galao dev
echo "========Deploy to latte========"
git push latte dev
#echo "=======Deploy to americano========"
#git push americano dev
#git push hyak dev

