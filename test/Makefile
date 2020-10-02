FLAGS = -I ../include/ -g -Wall -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=leak -fsanitize=undefined

SRCS := ${wildcard *.cpp}
OBJS := ${SRCS:%.cpp=%.o}

all: ${OBJS}
	g++ -L.. -Wl,-rpath,.. ${FLAGS} ${OBJS} -o test -ltecio -lcartdg

%.o: %.cpp
	g++ -c ${FLAGS} $< -o $@

clean:
	rm *.o
	rm test
	rm *.szplt

