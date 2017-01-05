INSTANCE = input/Co_100_1_T_7.txt
THREADS_NO = 4

EXE_NAME = Coiote_heuristic
OUTPUT_PATH = output/summary.csv
SOLUTIONS_PATH = solutions

CXXFLAGS += heuristic.cpp main.cpp -o $(EXE_NAME) -fpermissive -std=c++11 -pthread -O3

ifeq ($(THREADS_NO), 2)
	CXXFLAGS += -DTWO_THREADS
else ifeq ($(THREADS_NO), 3)
	CXXFLAGS += -DTHREE_THREADS
endif

build:
	$(CXX) $(CXXFLAGS) 
	
run:
	./$(EXE_NAME) -i $(INSTANCE) -o $(OUTPUT_PATH)

feasibility-check:
	mkdir -p $(SOLUTIONS_PATH)
	./$(EXE_NAME) -i $(INSTANCE) -o $(OUTPUT_PATH)  -s $(SOLUTIONS_PATH)/sol.txt
	./$(EXE_NAME) -i $(INSTANCE) -o $(OUTPUT_PATH)  -s $(SOLUTIONS_PATH)/sol.txt -test
