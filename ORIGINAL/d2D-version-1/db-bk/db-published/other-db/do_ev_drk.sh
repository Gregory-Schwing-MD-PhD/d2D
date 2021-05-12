for i in `jot - 0.80000 3.00000 0.200000`; do
for j in `jot - 0.80000 3.00000 0.200000`; do
for k in `jot - 0.80000 3.00000 0.200000`; do
for l in `jot - 0.10000 0.90000 0.100000`; do
sed 's/_SEEEEA_/'$i'/' top-weights.tab > tmp 
sed 's/_SEEEEB_/'$j'/' tmp > tmp2
sed 's/_SEEEEC_/'$k'/' tmp2 > tmp3
sed 's/_SEEEED_/'$l'/' tmp3 > weights.tab
cd ../../../ENSEMBLES/drkensembles

d2D.x -file cs-sparta.dat  >& /dev/null
grep -v "#" SS-results.dat | tail -n 57 > tmp
echo -n "$i $j $k $l " > test
paste tmp ss-a.dat | awk '{sh+=($3-$9)^2; sp+=($6-$12)^2; se+=($5-$11)^2; sc+=($4-$10)^2; c++} END {printf("%lf %lf %lf %lf ", sqrt(sh/c), sqrt(sp/c), sqrt(se/c), sqrt(sc/c))}' >> test
paste tmp ss-a.dat | awk '{print $3, $9}' | awk -f ~/bin/cor.awk | awk '{printf("%lf ", $1)}' >> test
paste tmp ss-a.dat | awk '{print $6, $12}' | awk -f ~/bin/cor.awk | awk '{printf("%lf ", $1)}' >> test
paste tmp ss-a.dat | awk '{print $4, $10}' | awk -f ~/bin/cor.awk | awk '{printf("%lf ", $1)}' >> test
paste tmp ss-a.dat | awk '{print $5, $11}' | awk -f ~/bin/cor.awk | awk '{printf("%lf \n", $1)}' >> test

d2D.x -file cs-shiftx2.dat  >& /dev/null
grep -v "#" SS-results.dat | tail -n 57 > tmp
echo -n "$i $j $k $l " >> test
paste tmp ss-a.dat | awk '{sh+=($3-$9)^2; sp+=($6-$12)^2; se+=($5-$11)^2; sc+=($4-$10)^2; c++} END {printf("%lf %lf %lf %lf", sqrt(sh/c), sqrt(sp/c), sqrt(se/c), sqrt(sc/c))}' >> test
paste tmp ss-a.dat | awk '{print $3, $9}' | awk -f ~/bin/cor.awk | awk '{printf("%lf ", $1)}' >> test
paste tmp ss-a.dat | awk '{print $6, $12}' | awk -f ~/bin/cor.awk | awk '{printf("%lf ", $1)}' >> test
paste tmp ss-a.dat | awk '{print $4, $10}' | awk -f ~/bin/cor.awk | awk '{printf("%lf ", $1)}' >> test
paste tmp ss-a.dat | awk '{print $5, $11}' | awk -f ~/bin/cor.awk | awk '{printf("%lf \n", $1)}' >> test
awk '{a+=$1; b+=$2; c+=$3; d+=$4; e+=$5; f+=$6; g+=$7; h+=$8; i+=$9; l+=$10; m+=$11; n+=$12} END {print a/2,b/2,c/2,d/2,e/2,f/2,g/2,h/2,i/2,l/2,m/2,n/2}' test
cd ~/Projects/0-Michele/1-d2D-predictor/d2D/d2D-program/other-db
done 
done
done
done
