import csv
import os
import re
import shlex
import subprocess
import time
import logging
import sys
from tqdm import tqdm

ZSTD_EXTENSIONS = ('.zst', '.zstd')

# awk in the chunking pipe exits with this sentinel code when it deliberately
# stops early (after printing the last requested line) instead of reading to
# EOF. wrap_chunked_extraction_command() uses it to tell a deliberate early stop
# apart from a genuine upstream extractor failure.
CHUNKING_EARLY_EXIT_CODE = 99

class OrfMUtils:
    def un_orfm_name(self, name):
        return re.sub(r'_\d+_\d+_\d+$', '', name)
    
    def un_orfm_start_frame_number(self, name):
        match = re.search(r'^.*_(\d+)_(\d+)_(\d+)$', name)
        start, frame, number = match.groups()[:4]
        return int(start), int(frame), int(number)

class TaxonomyFile:
    def __init__(self, taxonomy_file_path):
        self.sequence_to_taxonomy = {}
        utils = OrfMUtils()
        with open(taxonomy_file_path) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                self.sequence_to_taxonomy[\
                      utils.un_orfm_name(row[0])] = row[1]

    def __getitem__(self, item):
        return self.sequence_to_taxonomy[item]

    def merge(self, another_taxonomy_file):
        for key, value in another_taxonomy_file.sequence_to_taxonomy.items():
            if key not in self.sequence_to_taxonomy:
                self.sequence_to_taxonomy[key] = value


class FastaNameToSampleName:
    @staticmethod
    def fasta_to_name(query_sequences_file):
        sample_name = os.path.basename(query_sequences_file)
        if sample_name.endswith('.fifo'):
            sample_name = sample_name[:-5]
        compressed_extensions = ('.gz',) + ZSTD_EXTENSIONS
        # Put compressed extensions first so they are stripped off before FASTA/Q suffixes.
        for extension in compressed_extensions + ('.fna', '.fq', '.fastq', '.fasta', '.fa', '.sra'):
            if sample_name.endswith(extension):
                sample_name = sample_name[0:(len(sample_name)-len(extension))]
                if extension not in compressed_extensions:
                    break
        return sample_name


def prepare_zstd_fifos(file_paths, temp_dir, sleep_after_mkfifo=None):
    """Return (new_paths, processes) where any .zst/.zstd files are streamed into FIFOs.

    The caller is responsible for waiting on the returned processes.
    """
    prepared_paths = []
    processes = []
    for path in file_paths or []:
        if path is None:
            prepared_paths.append(path)
            continue
        if not path.endswith(ZSTD_EXTENSIONS):
            prepared_paths.append(path)
            continue

        base = os.path.basename(path)
        fifo_path = os.path.join(temp_dir, f"{base}.fifo")
        attempt = 1
        while os.path.exists(fifo_path):
            attempt += 1
            fifo_path = os.path.join(temp_dir, f"{base}.{attempt}.fifo")
        os.mkfifo(fifo_path)
        if sleep_after_mkfifo:
            time.sleep(sleep_after_mkfifo)
        cmd = f"zstdcat {shlex.quote(path)} > {shlex.quote(fifo_path)}"
        proc = subprocess.Popen(['bash', '-c', cmd], stderr=subprocess.PIPE, text=True)
        processes.append((proc, cmd))
        prepared_paths.append(fifo_path)
    return prepared_paths, processes

def add_chunking_pipe(read_chunk_size, read_chunk_number, lines_per_read=2):
    # Pipe to extract a chunk of reads.
    # Use awk rather than tail|head: awk exits actively after printing the last
    # required line, so it never blocks waiting for upstream EOF when the
    # remaining data is shorter than the requested chunk (e.g. the last chunk
    # of a file whose size is not a multiple of the chunk size).
    #
    # awk exits CHUNKING_EARLY_EXIT_CODE when it stops early after printing the
    # last requested line, and 0 when it instead reaches EOF first (the final or
    # a short chunk). wrap_chunked_extraction_command() relies on this to
    # distinguish a deliberate early stop from a genuine upstream failure.
    start_line = (read_chunk_size * (read_chunk_number - 1)) * lines_per_read + 1
    end_line = start_line + read_chunk_size * lines_per_read - 1
    logging.debug(f"Chunking pipe: awk NR>={start_line} to NR<={end_line} (chunk_size={read_chunk_size}, chunk_num={read_chunk_number}, lines_per_read={lines_per_read})")
    return f" | awk 'NR>={start_line}{{print; if(NR>={end_line})exit {CHUNKING_EARLY_EXIT_CODE}}}'"

def wrap_chunked_extraction_command(read_command, chunking_pipe, output_path):
    """Build a bash command running `read_command <chunking_pipe> > output_path`
    under `set -o pipefail`, but whose final exit status distinguishes a
    deliberate early stop from a genuine failure.

    When awk stops early (CHUNKING_EARLY_EXIT_CODE) the upstream extractor (e.g.
    kingfisher, gzip, cat) is killed mid-stream and exits non-zero via
    SIGPIPE/BrokenPipeError. Without this wrapper, pipefail would report that as
    a failure for every non-final chunk. The wrapper treats the awk sentinel as
    success, and otherwise propagates the real exit code so genuine extractor
    failures still surface.
    """
    inner = f"{read_command}{chunking_pipe} > {shlex.quote(output_path)}"
    # PIPESTATUS after the pipeline: [0]=read_command, [1]=awk.
    return (
        "set -o pipefail; "
        f"{inner}; "
        's=("${PIPESTATUS[@]}"); '
        f'[ "${{s[1]}}" -eq {CHUNKING_EARLY_EXIT_CODE} ] && exit 0; '
        '[ "${s[1]}" -ne 0 ] && exit "${s[1]}"; '
        'exit "${s[0]}"'
    )

def prepare_chunking_fifos(file_paths, temp_dir, read_chunk_size, read_chunk_number, sleep_after_mkfifo=None):
    """Return (new_paths, processes) where any files are streamed into FIFOs with chunking.

    The caller is responsible for waiting on the returned processes.
    """
    prepared_paths = []
    processes = []
    for path in file_paths or []:
        if path is None:
            prepared_paths.append(path)
            continue

        base = os.path.basename(path)
        for suffix in ('.gz', '.fifo') + ZSTD_EXTENSIONS:
            if base.endswith(suffix):
                base = base[:-len(suffix)]
        if base.endswith(('.fq', '.fastq')):
            lines_per_read = 4
        elif base.endswith(('.fa', '.fasta', '.fna')):
            lines_per_read = 2
        else:
            raise Exception(f"Cannot determine format (FASTA or FASTQ) for chunking from file extension: {path}")
        logging.debug(f"Detected {'FASTQ' if lines_per_read == 4 else 'FASTA'} format for {path} ({lines_per_read} lines/read)")
        chunking_dir = os.path.join(temp_dir, "chunking")
        os.makedirs(chunking_dir, exist_ok=True)
        fifo_path = os.path.join(chunking_dir, base)
        attempt = 1
        while os.path.exists(fifo_path):
            attempt += 1
            name, ext = os.path.splitext(base)
            fifo_path = os.path.join(chunking_dir, f"{name}.{attempt}{ext}")
        os.mkfifo(fifo_path)
        if sleep_after_mkfifo:
            # On kubernetes this seems to be required, at least in some circumstances
            logging.debug("Sleeping for {} seconds after mkfifo".format(sleep_after_mkfifo))
            time.sleep(sleep_after_mkfifo)
        prepared_paths.append(fifo_path)

        if path.endswith('.gz'):
            read_cmd = f"gzip -dc {shlex.quote(path)}"
        else:
            read_cmd = f"cat {shlex.quote(path)}"
        cmd = wrap_chunked_extraction_command(
            read_cmd,
            add_chunking_pipe(read_chunk_size, read_chunk_number, lines_per_read),
            fifo_path)
        logging.debug("Running chunking command: {}".format(cmd))
        process = subprocess.Popen(
            ['bash','-c',cmd],
            stdout=None,
            stderr=subprocess.PIPE,
            text=True)
        processes.append((process, cmd))

    return prepared_paths, processes

def terminate_processes(processes, description):
    """Send SIGTERM to any still-running background processes (e.g. after a downstream failure)."""
    for proc, _cmd in processes:
        if proc.poll() is None:
            logging.debug(f"Terminating {description} process (pid {proc.pid}) due to downstream failure")
            proc.terminate()

def finish_processes(processes, description):
    """Wait for streaming processes to finish and raise on non-zero exit codes."""
    for proc, cmd in processes:
        _, stderr_output = proc.communicate()
        # Negative codes mean the process was killed by a signal (e.g. our own
        # terminate_processes after a downstream failure) - tolerate those.
        # Deliberate early stops in chunked pipelines are already normalised to
        # exit 0 by wrap_chunked_extraction_command().
        if proc.returncode > 0:
            raise Exception(f"{description} command failed (exit {proc.returncode}): {cmd}\nSTDERR: {stderr_output}")

class LoggingTqdm(tqdm):
    """Drop-in replacement for tqdm that logs progress periodically when not in a TTY."""

    def __init__(self, *args, logger=None, log_interval=5, **kwargs):
        """
        Parameters
        ----------
        *args
            Passed through to tqdm (typically the iterable).
        logger : logging.Logger, optional
            Logger to write progress messages to. Defaults to the module logger.
        log_interval : int
            Minimum seconds between log messages. Default 60.
        **kwargs
            Passed through to tqdm. 'disable' is set automatically based on
            whether stderr is a TTY, but can be overridden explicitly.
        """
        self._logger = logger or logging.getLogger(__name__)
        self._log_interval = log_interval
        self._last_log = 0
        self._counter_for_logging = 0
        kwargs.setdefault('disable', not sys.stderr.isatty())
        super().__init__(*args, **kwargs)

    def update(self, n=1):
        """
        Advance the progress counter by n steps.

        When output is a TTY, delegates entirely to tqdm for live rendering.
        When output is not a TTY (e.g. log file), emits a log line at most
        once per log_interval seconds.

        Parameters
        ----------
        n : int
            Number of steps to advance. Default 1.
        """
        super().update(n)
        if self.disable:
            # tqdm doesn't update .n when disabled, so we track it separately for logging purposes.
            self._counter_for_logging += n
            now = time.time()
            
            if now - self._last_log > self._log_interval:
                pct = f"{100*self._counter_for_logging//self.total}%" if self.total else "?"
                self._logger.info(f"Progress: {self._counter_for_logging}/{self.total} ({pct})")
                self._last_log = now