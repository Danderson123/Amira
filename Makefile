# Inspired by https://github.com/snakemake/snakefmt/blob/master/Makefile
PROJECT = amira_prototype
OS := $(shell uname -s)
VERSION := $(shell poetry version -s)
BOLD := $(shell tput bold)
NORMAL := $(shell tput sgr0)

.PHONY: all
all: install

.PHONY: install
install:
	poetry lock
	poetry install

.PHONY: install-ci
install-ci:
	poetry install --no-interaction
	poetry run amira --help

.PHONY: pre-commit
pre-commit:
	poetry run pre-commit autoupdate
	poetry run pre-commit run --all-files -v

.PHONY: test
test:
	pip install .
	pytest tests

build:
	poetry build