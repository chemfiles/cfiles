# -*- codding: utf8 -*-
from chfltests import cfiles


def help_flag():
    stdout, stderr = cfiles("--help")

    assert(stderr == "")
    assert("help <subcommand>" in stdout)

    stdout, stderr = cfiles("-h")

    assert(stderr == "")
    assert("help <subcommand>" in stdout)


def help_command():
    stdout, stderr = cfiles("help")

    assert(stderr == "")
    assert("help <subcommand>" in stdout)


if __name__ == '__main__':
    help_flag()
    help_command()
