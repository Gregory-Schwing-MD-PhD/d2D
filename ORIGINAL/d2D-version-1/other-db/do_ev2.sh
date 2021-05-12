for i in `jot - 0.00000 3.50000 0.010000`; do
sed 's/_SEEEED_/'$i'/' top-weights.tab > weights.tab
cd ../../../ENSEMBLES/Sintetico/AB40
d2D.l.x -file cs-sparta.dat -ppii -nopsipred >& /dev/null
grep -v "#" SS-results.dat | tail -n 40 > tmp
echo -n "$i "
paste tmp l-count.dat | awk '{sh+=($3-$9)^2; sp+=($6-$12)^2; se+=($4-$11)^2; c++} END {printf("%lf %lf %lf ", sqrt(sh/c), sqrt(sp/c), sqrt(se/c))}'
paste tmp l-count.dat | awk '{print $3, $9}' | awk -f ~/bin/cor.awk | awk '{printf("%lf ", $1)}'
paste tmp l-count.dat | awk '{print $6, $12}' | awk -f ~/bin/cor.awk | awk '{printf("%lf ", $1)}'
paste tmp l-count.dat | awk '{print $4, $11}' | awk -f ~/bin/cor.awk | awk '{printf("%lf \n", $1)}'
cd ../../../d2D/programs/other-db
#cd ../../test-refDB
#./do_ref.sh
#awk '{sum+=$3; c++} END {print 100-sum/c}' refDBb.dat
#cd ../programs/other-db
done
