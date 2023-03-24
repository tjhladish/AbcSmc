.FORCE:

build: .FORCE
	cmake -S . -B $@

clean: .FORCE
	git clean -ifdx -e .vscode