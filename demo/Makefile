BUILD_DIR ?= ../build
OUT_DIR ?= outputs
CMAKE ?= cmake
PYTHON ?= python3

.PHONY: all demo build run plot clean

all: demo

build:
	$(CMAKE) -S .. -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -DCPPQUANT_BUILD_DEMO=ON
	$(CMAKE) --build $(BUILD_DIR) --target surface_demo

run: build
	$(BUILD_DIR)/surface_demo $(OUT_DIR)

plot: run
	$(PYTHON) plot.py --data-dir $(OUT_DIR) --out-dir $(OUT_DIR)

demo: plot

clean:
	rm -f $(OUT_DIR)/*.csv $(OUT_DIR)/*.png