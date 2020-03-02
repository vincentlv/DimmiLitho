
install: 
	pip install -r requirements.txt --upgrade
	pip install -e .

hook-black:
	cp .hooks/pre-commit .git/hooks/pre-commit

hook-pytest:
	cp .hooks/pre-push .git/hooks/pre-push

unhook:
	rm .git/hooks/*

lint:
	pyflakes dimmilitho 

