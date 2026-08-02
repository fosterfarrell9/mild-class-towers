PARI ?= $(HOME)/.local
CC ?= cc
CFLAGS ?= -O2 -g -std=c11 -Wall -Wextra -Wpedantic -fno-strict-aliasing
CPPFLAGS += -I$(PARI)/include
LDFLAGS += -L$(PARI)/lib -Wl,-rpath,$(PARI)/lib
LDLIBS += -lpari -lm

.PHONY: all clean

CORE_OBJECTS = $(filter-out ../build/massey.o, $(wildcard ../build/*.o))

all: cm_construct audit_cm_fields cm_character_driver

cm_construct: cm_construct.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

audit_cm_fields: audit_cm_fields.c
	$(CC) $(CPPFLAGS) -I.. $(CFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

cm_character_driver: ../parallelization/character_driver.c \
		bnrclassfield_cm_wrap.c $(CORE_OBJECTS)
	$(CC) $(CPPFLAGS) -I.. -I../headers $(CFLAGS) \
		../parallelization/character_driver.c bnrclassfield_cm_wrap.c \
		$(CORE_OBJECTS) $(LDFLAGS) -Wl,--wrap=bnrclassfield $(LDLIBS) -pthread -o $@

clean:
	$(RM) cm_construct audit_cm_fields cm_character_driver
