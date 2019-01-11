pipenv run python -m cProfile -s tottime tests.py > report.txt
head -40 report.txt
