graph_reduce:
	g++ gr.cpp -O2 -o gr -I/opt/ibm/ILOG/CPLEX_Studio129/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio129/concert/include -DIL_STD -L/opt/ibm/ILOG/CPLEX_Studio129/concert/lib/x86-64_linux/static_pic -L/opt/ibm/ILOG/CPLEX_Studio129/cplex/lib/x86-64_linux/static_pic  -lilocplex -lconcert -lcplex -lm -lpthread -ldl
