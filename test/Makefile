FLAG = -g -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=leak -fsanitize=undefined 

all: main
	g++ ${FLAG} main.o -o unit_test

clean:
	rm *.o
	rm unit_test

main: main.cpp
	g++ -c ${FLAG} main.cpp
