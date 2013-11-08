CC=g++
#CFLAGS=-std=c++11 -pg
CFLAGS=-std=c++11
#CFLAGS=-std=c++11 -g -pg

all: modexp basic_impl exp_opt mult_opt

modexp: modexp.cpp
	$(CC) $(CFLAGS) modexp.cpp -o modexp

basic_impl: basic_impl.cpp
	$(CC) $(CFLAGS) basic_impl.cpp -o basic_impl

exp_opt: exp_opt.cpp
	$(CC) $(CFLAGS) exp_opt.cpp -o exp_opt

mult_opt: mult_opt.cpp
	$(CC) $(CFLAGS) mult_opt.cpp -o mult_opt

clean:
	rm -rf modexp basic_impl exp_opt mult_opt

run:
	./modexp
	gprof modexp gmon.out > g.txt
