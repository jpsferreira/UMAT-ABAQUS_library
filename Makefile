# Get all directories with a number at the beginning. ignore .git
DIRS := $(shell find . -type d -name "[0-9]*" ! -path "./.git/*")

# Rules to build and test
.PHONY: all build test test_in_abaqus clean $(DIRS)

all: build test test_in_abaqus

build: $(addsuffix _build,$(DIRS))
test: $(addsuffix _test,$(DIRS))
test_in_abaqus: $(addsuffix _test_in_abaqus,$(DIRS))

%_build: 
	cd $(@:_build=) && sh build*.sh

%_test: 
	cd $(@:_test=) && ./*.o

%_test_in_abaqus: 
	cd $(@:_test_in_abaqus=)/test_in_abaqus ; INP_FILE=`ls cube*.inp | head -n 1` ; FOR_FILE=`ls umat*.f* | head -n 1` ; abaqus job=$${INP_FILE%.*} user=$${FOR_FILE%.*} interactive ask_delete=OFF

clean:
	for dir in $(DIRS); do \
		$(RM) $$dir/*.o ; \
		find $$dir/test_in_abaqus -type f ! \( -name "*.inc" -o -name "*.inp" -o -name "*.py" -o -name "*.for" \) -delete ; \
	done