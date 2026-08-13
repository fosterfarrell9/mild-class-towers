PARI ?= /usr/local
CC ?= cc
CFLAGS ?= -O2 -Wall -Wextra -fno-strict-aliasing

all: verify_certificate

verify_certificate: verify_certificate.c
	$(CC) $(CFLAGS) -I$(PARI)/include -o $@ $< \
		-Wl,-rpath,$(PARI)/lib -L$(PARI)/lib -lpari -lm

.PHONY: all check clean
check: verify_certificate
	./verify_certificate ../certificates/K-3640387-p3/certificate.gp
	python3 test_rejections.py

clean:
	$(RM) verify_certificate
