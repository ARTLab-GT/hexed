FLAGS = -I ../include/ -g -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=leak -fsanitize=undefined

SRCS := ${wildcard ../src/*.cpp}
SRCS := ${wildcard *.cpp} ${SRCS}
OBJS := ${SRCS:%.cpp=%.o}

all: ${OBJS}
	g++ ${FLAGS} ${OBJS} -o unit_test -ltecio

%.o: %.cpp
	g++ -c ${FLAGS} $< -o $@

clean:
	rm *.o
	rm unit_test
	rm *.szplt

