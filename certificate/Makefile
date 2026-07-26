PARI ?= /usr/local
CC ?= cc
CFLAGS ?= -O2 -Wall -Wextra -fno-strict-aliasing

verify_certificate: verify_certificate.c
	$(CC) $(CFLAGS) -I$(PARI)/include -o $@ $< \
		-Wl,-rpath,$(PARI)/lib -L$(PARI)/lib -lpari -lm

.PHONY: clean
clean:
	$(RM) verify_certificate
