.FORCE:

.FORCE:

build: .FORCE
	@cmake -S . -B $@ > /dev/null

clean: .FORCE
	git clean -ifdx -e .vscode

build/examples/%.sqlite: build
	cd $^ && make $* --quiet > /dev/null && cd examples && ./$* $*_config.json --process --simulate -n 1000
