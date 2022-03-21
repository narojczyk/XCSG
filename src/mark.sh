#!/bin/bash
chkFile="checksum.h"
mkFile="Makefile"
echo "Updating source SHA1 hashes in $chkFile"

gitName=`git config user.name 2>/dev/null`
bldrName=${gitName:-`whoami`}

gitMail=`git config user.email 2>/dev/null`
bldrMail=${gitMail:-`whoami`"@localhost"}

bldHost=`uname -n`

gitTag=`git tag -l | tail -n 1`
gitVersion=`echo ${gitTag} | sed 's/^v//'`
gitTagCommit=`git rev-list -n 1 ${gitTag}`
gitTagCommitNum=`git rev-list --count ${gitTagCommit}`
gitHeadCommitNum=`git rev-list --count HEAD`
let commAfetrTag=gitHeadCommitNum-gitTagCommitNum
gitVersion=`echo "${gitVersion}.${commAfetrTag}" | sed 's/\.$//'`
buildVersion=${gitVersion:-"no version tag in git"}

gitComitDate=`git log -1 --format=%ci 2>/dev/null`
commitDate=${gitComitDate:-"\(unspecified\)"}

gitComitId=`git rev-parse HEAD 2>/dev/null`
commitId=${gitComitId:-"\(unspecified\)"}

gitBranchName=`git rev-parse --abbrev-ref HEAD 2>/dev/null`
# use `git branch –-show-current` in git version > 2.22
branch=${gitBranchName:-"\(unknown\)"}

ccVendor=`cat ${mkFile} | grep ^CC | sed 's/^.*\ =\ //'`
ccVersion=`${ccVendor} -v 2>&1 |\
  grep ' version ' | eval "sed 's/^.*${ccVendor}\ version\ //'"`
ccFlags=`cat ${mkFile} |\
  grep ^CFLAGS | sed 's/^.*=\ /\ /' | tr -d '\n' | sed 's/^\ //'`;

# List of files to hash
fileList=(`ls -1 *.[ch] | grep -v $chkFile`)
N=${#fileList[@]}

eval "sed -i 's/\(builder\).*$/\1\ =\ \"$bldrName\ <$bldrMail>\";/' $chkFile"
eval "sed -i 's/\(buildAt\).*$/\1\ =\ \"$bldHost\";/' $chkFile"
eval "sed -i 's/\(build\)\ =.*\"/\1\ =\ \"`date`\"/' $chkFile"
eval "sed -i 's/\(code_version\).*$/\1\ =\ \"$buildVersion\";/' $chkFile"
eval "sed -i 's/\(code_commit_id\).*$/\1\ =\ \"$commitId\";/' $chkFile"
eval "sed -i 's/\(code_branch_name\).*$/\1\ =\ \"${branch}\";/' $chkFile"
eval "sed -i 's/\(code_commit_date\).*$/\1\ =\ \"${commitDate}\";/' $chkFile"
eval "sed -i 's/\(cc_vendor\).*$/\1\ =\ \"${ccVendor}\";/' $chkFile"
eval "sed -i 's/\(cc_version\).*$/\1\ =\ \"${ccVersion}\";/' $chkFile"
eval "sed -i 's/\(cc_flags\).*$/\1\ =\ \"${ccFlags}\";/' $chkFile"

i=0
while [ $i -lt $N ]
do
  file_i=`echo ${fileList[$i]} | sed 's/\.\([ch]\)/_\1_SHA1/'`

  lineNumber=`cat -n $chkFile|grep $file_i |sed 's/\tco.*//'|sed 's/\ \ *//'`
#   printf " %-18s\t%d\n" $file_i $lineNumber; echo

  if [ $lineNumber ]
  then
    chk_string=`sha1sum ${fileList[$i]} | sed 's/\ .*//'`
#     echo " ${fileList[$i]} $chk_string"
    eval "sed -i '${lineNumber}s/\".*\"/\"${chk_string}\"/'" $chkFile
  fi

  (( i++ ))
done

exit 0

