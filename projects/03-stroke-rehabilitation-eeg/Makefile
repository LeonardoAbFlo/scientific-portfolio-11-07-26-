.PHONY: install test pipeline figures

install:
	python -m pip install -e .

test:
	pytest

pipeline:
	python scripts/run_eeg_pipeline.py --data-dir data/raw --output-dir results

figures:
	python scripts/stage_05_create_figures.py --output-dir results
