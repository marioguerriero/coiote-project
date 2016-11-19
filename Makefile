build:
		g++ heuristic.cpp main.cpp -o Coiote_heuristic -fpermissive -std=c++11 -lga

build-tune:
		g++ Tuner.cpp heuristic.cpp -o tuner -fpermissive -std=c++11 -lga

run:
		./Coiote_heuristic -i input/Co_300_20_NT_0.txt -o output/summary.csv

tune:
		./tuner

feasibility-check:
		./Coiote_heuristic -i input/Co_300_20_NT_0.txt -o output/summary.csv  -s solutions/30.txt
		./Coiote_heuristic -i input/Co_300_20_NT_0.txt -o output/summary.csv  -s solutions/30.txt -test
