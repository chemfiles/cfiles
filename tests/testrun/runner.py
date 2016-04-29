# -*- codding: utf8 -*-
from subprocess import Popen, PIPE


class CfilesError(Exception):
    pass


def cfiles(*args):
    command = ["./cfiles"]
    command.extend(args)
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise CfilesError("Process exited with non-zero return code")

    return stdout, stderr
