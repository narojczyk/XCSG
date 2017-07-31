#!/bin/bash
echo "Computing source files SHA1 checksums"
eval "sed -i 's/build\ =.*\"/build\ =\ \"`date`\"/' checksum.h"

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

src_base=(`ls -1 *.[ch] | grep -v checksum.h`)
N=${#src_base[@]}

eval "sed -i 's/builder.*$/builder\ =\ \"$builderName\ <$builderMail>\";/' checksum.h"
eval "sed -i 's/buildAt.*$/buildAt\ =\ \"$buildHost\";/' checksum.h"

i=0
while [ $i -lt $N ]
do
  var_name=`echo ${src_base[$i]} | sed 's/\.\([ch]\)/_\1_SHA1/'`

  ln=`cat -n checksum.h|grep $var_name |sed 's/\tch.*//'|sed 's/\ \ *//'`
#   printf " %-18s\t%d\n" $var_name $ln; echo

  if [ $ln ]
  then
    chk_string=`sha1sum ${src_base[$i]} | sed 's/\ .*//'`
#     echo " ${src_base[$i]} $chk_string"

    eval "sed -i '${ln}s/\".*\"/\"${chk_string}\"/'" checksum.h
  fi

  (( i++ ))
done

exit 0

