
PHONY:=all test clean
all: test

test: test/test-svgtoipe.ipe test/test-markers.ipe

test/%.ipe: test/%.svg
	python svgtoipe.py -v1 $< $@

clean: 
	rm -f test/test-svgtoipe.ipe
	rm -f test/test-markers.ipe
