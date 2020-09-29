FLAGS = -I ../include/ -O3 -march=native

SRCS := ${wildcard ../src/*.cpp}
SRCS := ${wildcard *.cpp} ${SRCS}
OBJS := ${SRCS:%.cpp=%.o}

all: ${OBJS}
	g++ ${FLAGS} ${OBJS} -o demo -ltecio

%.o: %.cpp
	g++ -c ${FLAGS} $< -o $@

clean:
	rm *.o
	rm ../src/*.o
	rm demo
	rm *.szplt

