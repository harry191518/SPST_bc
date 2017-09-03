CC = gcc
CFLAGS = -O4
INC = -I include
INC_PATH = include

a.out: main.o header.o structure.o function.o result.o
	${CC} main.o header.o structure.o function.o result.o ${CFLAGS} ${INC} -o a.out
main.o: main.c ${INC_PATH}/header.h ${INC_PATH}/structure.h ${INC_PATH}/function.h ${INC_PATH}/result.h
	${CC} main.c ${CFLAGS} ${INC} -c
header.o: header.c ${INC_PATH}/header.h 
	${CC} header.c ${CFLAGS} ${INC} -c
structure.o: structure.c ${INC_PATH}/header.h ${INC_PATH}/structure.h ${INC_PATH}/function.h
	${CC} structure.c ${CFLAGS} ${INC} -c 
function.o: function.c ${INC_PATH}/structure.h ${INC_PATH}/function.h 
	${CC} function.c ${CFLAGS} ${INC} -c
result.o: result.c ${INC_PATH}/header.h ${INC_PATH}/structure.h ${INC_PATH}/function.h ${INC_PATH}/result.h
	${CC} result.c ${CFLAGS} ${INC} -c

clean:
	@rm -rf *.o a.out
