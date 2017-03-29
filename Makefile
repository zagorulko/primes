SOURCES=main.c find_prime.c
OUTPUT=prime
LDFLAGS=-lgmp -lm -lpthread

default: build

build:
	gcc -Wall -g $(SOURCES) $(LDFLAGS) -o $(OUTPUT) -O3
