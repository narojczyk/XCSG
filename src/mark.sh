#!/bin/bash
chkFile="checksum.h"
echo "Updating source SHA1 hashes in $chkFile"


builderName=`git config user.name 2>/dev/null`
if [ $? -ne 0 ]; then
  builderName=`whoami`
  echo $builderName
fi
builderMail=`git config user.email 2>/dev/null`
if [ $? -ne 0 ]; then
  builderMail="unknown"
fi

buildHost=`uname -n`

buildVersion=`git tag -l | tail -n 1 | sed 's/^v//'`
if [ $? -ne 0 ]; then
  buildVersion="untagged in git"
fi

comitDate=`git log -1 --format=%ci 2>/dev/null`
if [ $? -ne 0 ]; then
  comitDate="no git repository found"
fi

comitId=`git rev-parse HEAD 2>/dev/null`
if [ $? -ne 0 ]; then
  comitId="no git repository found"
fi

src_base=(`ls -1 *.[ch] | grep -v $chkFile`)
N=${#src_base[@]}

eval "sed -i 's/builder.*$/builder\ =\ \"$builderName\ <$builderMail>\";/' $chkFile"
eval "sed -i 's/buildAt.*$/buildAt\ =\ \"$buildHost\";/' $chkFile"
eval "sed -i 's/build\ =.*\"/build\ =\ \"`date`\"/' $chkFile"
eval "sed -i 's/code_version.*$/code_version\ =\ \"$buildVersion\";/' $chkFile"
eval "sed -i 's/code_comit_id.*$/code_comit_id\ =\ \"$comitId\";/' $chkFile"
eval "sed -i 's/code_comit_date.*$/code_comit_date\ =\ \"$comitDate\";/' $chkFile"

i=0
while [ $i -lt $N ]
do
  var_name=`echo ${src_base[$i]} | sed 's/\.\([ch]\)/_\1_SHA1/'`

  ln=`cat -n $chkFile|grep $var_name |sed 's/\tch.*//'|sed 's/\ \ *//'`
#   printf " %-18s\t%d\n" $var_name $ln; echo

  if [ $ln ]
  then
    chk_string=`sha1sum ${src_base[$i]} | sed 's/\ .*//'`
#     echo " ${src_base[$i]} $chk_string"

    eval "sed -i '${ln}s/\".*\"/\"${chk_string}\"/'" $chkFile
  fi

  (( i++ ))
done

exit 0

