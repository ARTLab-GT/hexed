#FLAGS = -I include/ -g -Wall -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=leak -fsanitize=undefined
FLAGS = -I include/ -O3 -march=native

SRCS := ${wildcard src/*.cpp}
OBJS := ${SRCS:%.cpp=%.o}

library: ${OBJS}
	g++ -shared ${FLAGS} ${OBJS} -o libcartdg.so -ltecio

%.o: %.cpp
	g++ -c -fpic ${FLAGS} $< -o $@

clean:
	rm src/*.o
	rm libcartdg.so

