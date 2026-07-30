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

# Multithreaded variant: a second PARI prefix built with --mt=pthread.
# The core sources are recompiled against its headers into a private
# object directory, so neither src/ nor the production build is touched.
PARI_MT ?= $(HOME)/.local-mt
MTOBJDIR = build-mt
CORE_SOURCES = $(filter-out ../src/massey.c, $(wildcard ../src/*.c))
MT_OBJECTS = $(patsubst ../src/%.c,$(MTOBJDIR)/%.o,$(CORE_SOURCES))
MT_CPPFLAGS = -I.. -I../headers -I$(PARI_MT)/include

all: character_driver finish_driver

$(MTOBJDIR):
	mkdir -p $(MTOBJDIR)

$(MTOBJDIR)/%.o: ../src/%.c | $(MTOBJDIR)
	$(CC) -c -o $@ $< $(CFLAGS) $(MT_CPPFLAGS)

character_driver_mt: character_driver.c $(MT_OBJECTS)
	$(CC) $(CFLAGS) $(MT_CPPFLAGS) -o $@ character_driver.c $(MT_OBJECTS) \
		-Wl,-O3 -Wl,-rpath,$(PARI_MT)/lib -L$(PARI_MT)/lib -lpari -lm -pthread

character_driver: character_driver.c $(CORE_OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ character_driver.c $(CORE_OBJECTS) $(LDFLAGS) $(LIBS)

finish_driver: finish_driver.c $(CORE_OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ finish_driver.c $(CORE_OBJECTS) $(LDFLAGS) $(LIBS)

.PHONY: all clean
clean:
	$(RM) character_driver finish_driver
