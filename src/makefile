CC = g++

CCFLAGS = -O3
OBJS = mdl.o

all	: mdl

clean   :
	rm *.o mdl

mdl  :	${OBJS}
	${CC} -o $@ ${CCFLAGS} ${OBJS}

mdl.o :	mdl.h mdl.C
	${CC} -c ${CCFLAGS} mdl.C

