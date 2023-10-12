.FORCE:

.FORCE:

build: .FORCE
	@mkdir $@ && cd $@
	@cmake ..

install: build
	@cd $^ && make install

clean: .FORCE
	git clean -ifdx -e .vscode

build/examples/%.sqlite: build
	cd $^ && make $* --quiet > /dev/null && cd examples && ./$* $*_config.json --process --simulate -n 1000
