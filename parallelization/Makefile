# Builds the per-character parallelization drivers against the existing,
# unmodified object files of the main build.  Run `make PARI=<prefix>` in
# the repository root first so that ../build/*.o exist.

PARI ?= /usr/local
CC = clang
CPPFLAGS = -I.. -I../headers -I$(PARI)/include
CFLAGS = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -pipe -march=native -pthread -g
LDFLAGS = -Wl,-O3 -Wl,-rpath,$(PARI)/lib -pthread
LIBS = -L$(PARI)/lib -lpari -lm

CORE_OBJECTS = $(filter-out ../build/massey.o, $(wildcard ../build/*.o))

all: character_driver finish_driver

character_driver: character_driver.c $(CORE_OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ character_driver.c $(CORE_OBJECTS) $(LDFLAGS) $(LIBS)

finish_driver: finish_driver.c $(CORE_OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ finish_driver.c $(CORE_OBJECTS) $(LDFLAGS) $(LIBS)

.PHONY: all clean
clean:
	$(RM) character_driver finish_driver
