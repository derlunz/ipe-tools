
PHONY:=all test clean
all: test

test: test/test-svgtoipe.ipe -v1 test/test-markers.ipe

test/%.ipe: test/%.svg
	python svgtoipe.py -v0	$< $@

clean: 
	rm -f test/test-svgtoipe.ipe
	rm -f test/test-markers.ipe
