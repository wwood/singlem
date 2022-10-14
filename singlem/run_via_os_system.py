import tempfile
import os
import logging

class RunViaOsSystemException(Exception):
    def __init__(self, cmd, stderr_tempfile):
        self.cmd = cmd
        self.stderr_tempfile = stderr_tempfile

    def __str__(self):
        msg = "Command failed: '{}'.".format(self.cmd)
        if self.stderr_tempfile is not None:
            with open(self.stderr_tempfile.name) as f:
                msg = "{} STDERR was '{}'".format(msg, f.read())
        return msg

def run_via_os_system(cmd, add_stderr=True, collect_stdout=True):
    '''Running with os.system() rather than subprocess is rarely useful. It is
    useful in instances where a file descriptor presented to the parent python
    process needs to be used in a child process.

    e.g. 
    $ script.py <(echo yes)
    where inside script.py is:
    extern.run_via_os_system('cat {}'.format(sys.argv[1]))

    Code showing various other methods failing:
    # cmd = "cat {}".format(sys.argv[1])
    # logging.info(extern.run("cat {}".format(sys.argv[1])))
    # => failed
    # logging.info(subprocess.check_output("cat {}".format(sys.argv[1])))
    # => failed
    # os.system("cat {}".format(sys.argv[1]))
    # => works
    # process = os.popen(cmd)
    # logging.info(process.read())
    # => fails
    # logging.info(subprocess.check_output("cat {}".format(sys.argv[1]), shell=True))
    # => fails
    '''

    stderr_tf = None
    if add_stderr:
        stderr_tf = tempfile.NamedTemporaryFile(prefix='extern_run_via_os_system')
        cmd = "{} 2>{}".format(cmd, stderr_tf.name)
    if collect_stdout:
        stdout_tf = tempfile.NamedTemporaryFile(prefix='extern_run_via_os_system')
        cmd = "{} >{}".format(cmd, stdout_tf.name)
    logging.debug("Running os.system command via tempfile: {}".format(cmd))
    with tempfile.NamedTemporaryFile(prefix='os.system') as script_tf:
        script_tf.write(cmd.encode())
        script_tf.flush()

        wait_status = os.system('bash {}'.format(script_tf.name))
        if wait_status != 0:
            raise RunViaOsSystemException(cmd, stderr_tf)
        if collect_stdout:
            with open(stdout_tf.name) as f:
                return f.read()
        else:
            return