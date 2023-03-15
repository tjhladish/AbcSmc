
.FORCE:

build: .FORCE
	cmake -S . -B $@
