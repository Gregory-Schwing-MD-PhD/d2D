for i in `jot - 0.00000 0.80000 0.00100`; do
sed 's/_SEEEED_/'$i'/' top-weights.tab > weights.tab
cd ../../../ENSEMBLES/drkensembles/A/sparta+
d2D.l.x -file cs-drkA.dat -ppii -nopsipred >& /dev/null
grep -v "#" SS-results.dat | tail -n 59 > tmp
echo -n "$i "
paste tmp ../../ss-count-A.dat | awk '{sh+=($3-$9)^2; sp+=($6-$12)^2; c++} END {printf("%lf %lf ", sqrt(sh/c), sqrt(sp/c))}'
paste tmp ../../ss-count-A.dat | awk '{print $3, $9}' | awk -f ~/bin/cor.awk | awk '{printf("%lf ", $1)}'
paste tmp ../../ss-count-A.dat | awk '{print $6, $12}' | awk -f ~/bin/cor.awk | awk '{printf("%lf \n", $1)}'
cd ../../../../d2D/programs/other-db
#cd ../../test-refDB
#./do_ref.sh
#awk '{sum+=$3; c++} END {print 100-sum/c}' refDBb.dat
#cd ../programs/other-db
done
