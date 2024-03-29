from subprocess import PIPE, Popen


class CfilesError(Exception):
    pass


def cfiles(*args):
    command = ["./cfiles"]
    command.extend(args)
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        command = " ".join(command)
        raise CfilesError(
            "Process '{}' exited with non-zero return code".format(command)
        )

    return stdout.decode("utf8"), stderr.decode("utf8")
