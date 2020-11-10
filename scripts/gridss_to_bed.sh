#!/bin/sh

# see https://stackoverflow.com/questions/59895/how-to-get-the-source-directory-of-a-bash-script-from-within-the-script-itself
DIR="$(dirname $0)"

perl ${DIR}/gridss2bed.pl $1 $2 $3

result=$?
# if [ $result -eq 0]; then # further processing
# print OK
# else
# print -u2 ERROR in xxx
exit $result
# fi