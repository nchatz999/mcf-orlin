CC = gcc
CFLAGS = -Wall -Wextra -g

BUILD = build

all: mcf_section3 mcf_section4

mcf_section3: main_section3.c
	$(CC) $(CFLAGS) -o $@ $<

mcf_section4: main_section4.c
	$(CC) $(CFLAGS) -o $@ $<

thesis: thesis.tex
	mkdir -p $(BUILD)
	latexmk -xelatex -shell-escape -output-directory=$(BUILD) thesis.tex
	cp $(BUILD)/thesis.pdf .

clean:
	rm -f mcf_section3 mcf_section4

clean-tex:
	rm -rf $(BUILD) _minted thesis.pdf

.PHONY: all clean clean-tex thesis
