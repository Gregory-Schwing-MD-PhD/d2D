for i in top-weights.*; do
cp $i weights.tab
cd ../../test-refDB
./do_ref.sh
awk '{sum+=$3; c++} END {print 100-sum/c, "'$i'"}' refDB.dat
cd ../programs/other-db
done
