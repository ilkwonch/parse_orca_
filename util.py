import re


class ParseError(Exception):
    pass


def goto(f, regex):
    for line in f:
        if re.search(regex, line):
            return line


def iter_until(f, regex):
    for line in f:
        if re.search(regex, line):
            break
        yield line

