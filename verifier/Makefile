PARI ?= /usr/local
CC ?= cc
CFLAGS ?= -O2 -Wall -Wextra -fno-strict-aliasing

all: verify_certificate

verify_certificate: verify_certificate.c
	$(CC) $(CFLAGS) -I$(PARI)/include -o $@ $< \
		-Wl,-rpath,$(PARI)/lib -L$(PARI)/lib -lpari -lm

.PHONY: all check clean
check: verify_certificate
	./verify_certificate ../certificates/p3/000/K-3640387-p3/certificate.gp
	./verify_certificate ../certificates/p5/K-2800905-p5/certificate.gp
	./verify_certificate ../certificates/p7/K-501510767-p7/certificate.gp
	python3 test_rejections.py

clean:
	$(RM) verify_certificate
