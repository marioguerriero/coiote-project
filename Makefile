build:
		g++ heuristic.cpp main.cpp -o Coiote_heuristic -fpermissive -std=c++11 -pthread -O3

run:
		./Coiote_heuristic -i input/Co_100_1_ST_0.txt -o output/summary.csv

feasibility-check:
		./Coiote_heuristic -i input/Co_30_1_ST_0.txt -o output/summary.csv  -s solutions/30.txt
		./Coiote_heuristic -i input/Co_30_1_ST_0.txt -o output/summary.csv  -s solutions/30.txt -test
