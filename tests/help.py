# -*- codding: utf8 -*-
from chfltests import cfiles


def help_command():
    stdout, stderr = cfiles("help")

    assert(stderr == "")
    assert("help <subcommand>" in stdout)


if __name__ == '__main__':
    help_command()
