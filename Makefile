test:
	python -m pytest tests/ -s

lint:
	python -m pylint lightecc/ --fail-under=10