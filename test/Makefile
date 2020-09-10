FLAGS = -I ../include/ -g -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=leak -fsanitize=undefined 

SRCS := Basis.cpp Equidistant.cpp
SRCS := ${SRCS:%=../src/%}
SRCS := ${wildcard *.cpp} ${SRCS}
OBJS := ${SRCS:%.cpp=%.o}

all: ${OBJS}
	g++ ${FLAGS} ${OBJS} -o unit_test

%.o: %.cpp
	g++ -c ${FLAGS} $< -o $@

clean:
	@echo ${SRCS}
	@echo ${OBJS}
	rm *.o
	rm unit_test

