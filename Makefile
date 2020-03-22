install:
	pip install -r requirements.txt --upgrade
	pip install -e .
	pip install pre-commit
	pre-commit install

lint:
	pyflakes litho
