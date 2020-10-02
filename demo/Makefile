FLAGS = -I ../include/ -O3 -march=native

all:
	g++ -L.. -Wl,-rpath,.. ${FLAGS} main.cpp -o demo -ltecio -lcartdg

clean:
	rm demo
	rm *.szplt

