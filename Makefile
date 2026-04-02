# UMAT-ABAQUS Library Makefile
#
# Primary workflow: generate material directories with generate.py
# This Makefile provides convenience targets for testing.

PYTHON = uv run python
EXAMPLES = neo_hooke mooney_rivlin ogden_3term humphrey_hgo humphrey_fiber \
           neo_hooke_damage humphrey_hgo_damage humphrey_fiber_damage \
           neo_hooke_visco mooney_rivlin_visco ogden_visco \
           humphrey_hgo_visco humphrey_fiber_visco \
           affine_network nonaffine_network

.PHONY: all examples test validate clean

all: test

# Generate all built-in examples
examples:
	@for ex in $(EXAMPLES); do \
		echo "=== Generating $$ex ===" ; \
		$(PYTHON) generate.py --example $$ex ; \
	done

# Generate and run standalone tests for all examples
test: examples
	@pass=0; fail=0; \
	for ex in $(EXAMPLES); do \
		echo "=== Testing $$ex ===" ; \
		if (cd $$ex && bash -o pipefail -c 'make run 2>&1 | tail -1'); then \
			pass=$$((pass+1)); \
		else \
			echo "FAIL: $$ex"; fail=$$((fail+1)); \
		fi; \
	done; \
	echo ""; echo "$$pass passed, $$fail failed"

# Run validation against legacy code
validate:
	$(PYTHON) validate.py

# Clean generated directories
clean:
	@for ex in $(EXAMPLES); do rm -rf $$ex; done
	rm -rf validation_tmp

# --- Legacy targets (soft_tissues/ and biofilaments/) ---
LEGACY_DIRS := $(shell find . -type d -name "[0-9]*" ! -path "./.git/*" 2>/dev/null)

legacy_build: $(addsuffix _build,$(LEGACY_DIRS))
legacy_test: $(addsuffix _test,$(LEGACY_DIRS))

%_build:
	cd $(@:_build=) && sh build*.sh

%_test:
	cd $(@:_test=) && ./*.o

legacy_clean:
	for dir in $(LEGACY_DIRS); do \
		$(RM) $$dir/*.o ; \
	done
