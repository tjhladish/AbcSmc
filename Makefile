.FORCE:

build: .FORCE
	mkdir -p $@ && cd $@ && cmake ..

build-shared: .FORCE
	mkdir -p $@ && cd $@ && cmake .. -DBUILD_SHARED_LIBS=YES -DCMAKE_BUILD_TYPE=Release && cmake --build .

build-static: .FORCE
	mkdir -p $@ && cd $@ && cmake .. -DBUILD_SHARED_LIBS=NO -DCMAKE_BUILD_TYPE=Release && cmake --build .

build-debug: .FORCE
	mkdir -p $@ && cd $@ && cmake .. -DBUILD_SHARED_LIBS=NO -DCMAKE_BUILD_TYPE=Debug && cmake --build .

install-exec: .FORCE
	cd build && make && sudo cmake --install .

clean: .FORCE
	git clean -ifdx -e .vscode

build/examples/%.sqlite: build
	cd $^ && make $* --quiet > /dev/null && cd examples && ./$* $*_config.json --process --simulate -n 1000
