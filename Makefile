all: projection

LDFLAGS=-lm

projection: projection.c

.PHONY=format
format: projection.c
	clang-format -i projection.c
