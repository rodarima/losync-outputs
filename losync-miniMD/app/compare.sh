#!/bin/bash
NAME=$1

# Sort local output and compare with known good
echo "Sorting test version of $NAME..."
sort -o $NAME -g -k1,1 -k2,2 -k3,3 $NAME
echo "Done"
echo "Sorting known good version of $NAME..."
sort -o ../miniMD_ref_clean/${NAME} -g -k1,1 -k2,2 -k3,3 ../miniMD_ref_clean/${NAME}
echo "Done"
diff $NAME ../miniMD_ref_clean/${NAME}

if [ $? -eq 0 ]; then
  echo "No differences found"
fi
