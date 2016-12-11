
PHONY:=all test clean
all: test

test: test/test-svgtoipe.ipe test/test-svgtoipe-normalsvg.ipe test/test-markers.ipe test/svg-tiger.ipe test/svg-tiger-mml.ipe

test/%.ipe: test/%.svg
	python svgtoipe.py -v2 $< $@

clean:
	rm -f test/test-svgtoipe.ipe
	rm -f test/test-markers.ipe
