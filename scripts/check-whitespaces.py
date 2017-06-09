#!/usr/bin/env python
# coding=utf-8
import os
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
DIRECTORIES = ["src", "tests"]
IGNORED_EXTENSIONS = ["pyc"]

ERRORS = 0


def error(message):
    global ERRORS
    ERRORS += 1
    print(message)


def check_whitespace(path):
    path = os.path.relpath(os.path.realpath(path), ROOT)
    with open(path) as fd:
        lines = fd.readlines()

    line_number = 0
    for line in lines:
        line_number += 1
        if line.endswith(" \n") or line.endswith("\t\n"):
            error(
                "whitespace at the end of the line at {}:{}"
                .format(path, line_number)
            )
    if not lines[-1].endswith("\n"):
        error("missing new line at the end of the file in {}".format(path))


if __name__ == '__main__':
    for directory in DIRECTORIES:
        for (root, _, paths) in os.walk(os.path.join(ROOT, directory)):
            for path in paths:
                extension = path.split(".")[-1]
                if extension not in IGNORED_EXTENSIONS:
                    check_whitespace(os.path.join(root, path))

    if ERRORS != 0:
        print("------------------\n{} whitespace errors".format(ERRORS))
        sys.exit(1)
