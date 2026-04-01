CC = gcc
CFLAGS = -Wall -Wextra -g

all: mcf_section3 mcf_section4

mcf_section3: main_section3.c
	$(CC) $(CFLAGS) -o $@ $<

mcf_section4: main_section4.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f mcf_section3 mcf_section4

.PHONY: all clean
