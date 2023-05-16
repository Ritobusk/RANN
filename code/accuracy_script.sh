#!/bin/bash

for i in  20; do
for (( c=5; c<=20; c +=5 ));
do
   i64s="i64"
   i32s="i32"
   var_ktmp="-b -g $i$i32s -g $c$i64s -g [100000][150]f32 -g [500][150]f32"
   futhark dataset $var_ktmp > tmpdata$c
   echo "Data set with k = $c , t = $i created. Acc for configuration of RANN is "
    ./knn-comparer -t time_with_k$c < tmpdata$c
done;
var_ktmp="-b -g $i$i32s -g 100$i64s -g [100000][150]f32 -g [500][150]f32"
futhark dataset $var_ktmp > tmpdata100
echo "Data set with k = 100, t = $i created. Acc for configuration of RANN is "
./knn-comparer -t time_with_k100 < tmpdata100
done
