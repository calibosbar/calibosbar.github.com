
#!/bin/bash
while read -r FILE
do
  FOUND="$(find . -name "$FILE" -print -quit)"
  if [ "x$FOUND" != "x" ]
  then
    echo "FOUND: $FILE"
  else
    echo "NOT FOUND: $FILE"
  fi
done <redshift2.txt
