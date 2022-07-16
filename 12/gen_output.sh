mkdir -p output
for i in {1..1000}
do
    ./monte_carlo > output/output_${i}.txt
done
