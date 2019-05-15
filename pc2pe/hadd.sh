# SK4
./hadd -f -t output/pc2pe_tst061889.root output/61889/*.root &
./hadd -f -t output/pc2pe_tst061892_to_5.root output/6189[2-5]/*.root &

# SK5
./hadd -f -t output/pc2pe_tst080871_to_5.root output/8087[1-5]/*.root &
./hadd -f -t output/pc2pe_tst080877.root output/80877/*.root &
./hadd -f -t output/pc2pe_tst080884_and_6.root output/8088[4,6]/*.root &
./hadd -f -t output/pc2pe_tst080885.root output/80885/*.root &

./hadd -f -t output/pc2pe_tst081028.root output/81028/*.root &
./hadd -f -t output/pc2pe_tst081030.root output/81030/*.root &

wait 

cp output/pc2pe_tst081028.root output/pc2pe_tst081028_kludge4avg.root
