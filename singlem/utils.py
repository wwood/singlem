import csv
import os
import re
import shlex
import signal
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


GZIP_EXTENSIONS = ('.gz',)

def prepare_zstd_fifos(file_paths, temp_dir, sleep_after_mkfifo=None, decompress_gzip=False):
    """Return (new_paths, processes) where any .zst/.zstd files (and .gz files,
    if decompress_gzip) are decompressed by a separate process into FIFOs.

    gzip input is decompressed this way because DIAMOND's own gzip reading is
    slow: it took 140s to load the queries of a 20M read FASTQ.gz file, versus
    18s when the same file was given to it already decompressed through zcat.

    The caller is responsible for waiting on the returned processes.
    """
    prepared_paths = []
    processes = []
    for path in file_paths or []:
        if path is None:
            prepared_paths.append(path)
            continue
        if path.endswith(ZSTD_EXTENSIONS):
            decompressor = 'zstdcat'
        elif decompress_gzip and path.endswith(GZIP_EXTENSIONS):
            decompressor = 'gzip -dc'
        else:
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
        cmd = f"{decompressor} {shlex.quote(path)} > {shlex.quote(fifo_path)}"
        proc = subprocess.Popen(['bash', '-c', cmd], stderr=subprocess.PIPE, text=True,
                                start_new_session=True)
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
        # In its own process group, so that terminate_processes() can kill the
        # whole pipeline, not just bash.
        process = subprocess.Popen(
            ['bash','-c',cmd],
            stdout=None,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True)
        processes.append((process, cmd))

    return prepared_paths, processes

def terminate_processes(processes, description):
    """Send SIGTERM to any still-running background processes (e.g. after a downstream failure)."""
    for proc, _cmd in processes:
        if proc.poll() is None:
            logging.debug(f"Terminating {description} process (pid {proc.pid}) due to downstream failure")
        # The process is the leader of its own process group (see
        # start_new_session above), so kill the whole group. Killing only bash
        # would leave the commands in its pipeline running, and those can block
        # forever opening a FIFO that is never going to be read, holding open the
        # stderr pipe that finish_processes() then waits on. This is done even
        # when bash itself has already exited, for the same reason.
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            if proc.poll() is None:
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
def first_sequence_lengths(path, num_sequences=5):
    '''Return the lengths of the first num_sequences sequences in a FASTA or
    FASTQ file, which may be gzip or zstd compressed. Fewer lengths are returned
    when the file has fewer sequences. Returns None when the file cannot be read
    from the start (e.g. a named pipe, which reading would consume).'''
    import gzip
    import stat
    if not os.path.exists(path) or not stat.S_ISREG(os.stat(path).st_mode):
        return None

    zstd_process = None
    if path.endswith(ZSTD_EXTENSIONS):
        zstd_process = subprocess.Popen(
            ['zstdcat', path], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        handle = zstd_process.stdout
    else:
        with open(path, 'rb') as f:
            is_gzip = f.read(2) == b'\x1f\x8b'
        handle = gzip.open(path, 'rt') if is_gzip else open(path)

    lengths = []
    try:
        first_line = handle.readline()
        if first_line.startswith('@'):
            # FASTQ, assumed unwrapped
            while len(lengths) < num_sequences:
                sequence = handle.readline()
                if sequence == '':
                    break
                lengths.append(len(sequence.strip()))
                handle.readline()
                handle.readline()
                if handle.readline() == '':
                    break
        elif first_line.startswith('>'):
            current = 0
            for line in handle:
                if line.startswith('>'):
                    lengths.append(current)
                    if len(lengths) >= num_sequences:
                        break
                    current = 0
                else:
                    current += len(line.strip())
            else:
                lengths.append(current)
    finally:
        handle.close()
        if zstd_process is not None:
            zstd_process.kill()
            zstd_process.wait()
    return lengths

def first_sra_spot_read_lengths(path, num_spots=5):
    '''Return the lengths of the biological reads in the first num_spots spots
    of an SRA file, as a list (one entry per spot) of lists of read lengths, by
    running vdb-dump. A spot with 2 or more biological reads is a read pair.
    Returns None if vdb-dump fails or its output cannot be parsed.'''
    try:
        output = subprocess.run(
            ['vdb-dump', '-C', 'READ_LEN,READ_TYPE', '-R', '1-{}'.format(num_spots), path],
            capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError) as e:
        logging.debug("Could not determine read lengths of {} with vdb-dump: {}".format(path, e))
        return None

    spots = []
    read_lengths = None
    for line in output.splitlines():
        line = line.strip()
        if line.startswith('READ_LEN:'):
            read_lengths = [int(x) for x in line.split(':', 1)[1].split(',') if x.strip() != '']
        elif line.startswith('READ_TYPE:') and read_lengths is not None:
            read_types = [x.strip() for x in line.split(':', 1)[1].split(',')]
            if len(read_types) != len(read_lengths):
                return None
            spots.append([length for (length, read_type) in zip(read_lengths, read_types)
                          if 'BIOLOGICAL' in read_type and length > 0])
            read_lengths = None
    if len(spots) == 0:
        return None
    return spots
