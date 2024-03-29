
install:
	pip install -r requirements.txt --upgrade
	pip install -r requirements_dev.txt --upgrade
	pip install -e .
	pre-commit install

test:
	pytest

cov:
	pytest --cov= litho

mypy:
	mypy . --ignore-missing-imports

lint:
	flake8

pylint:
	pylint litho

lintd2:
	flake8 --select RST

lintd:
	pydocstyle litho

doc8:
	doc8 docs/

update:
	pur

update2:
	pre-commit autoupdate --bleeding-edge
