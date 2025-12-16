import csv
import os
import re
import shlex
import subprocess
import time

ZSTD_EXTENSIONS = ('.zst', '.zstd')

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
        for extension in compressed_extensions + ('.fna','.fq','.fastq','.fasta','.fa'):
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


def finish_processes(processes, description):
    """Wait for streaming processes to finish and raise on non-zero exit codes."""
    for proc, cmd in processes:
        _, stderr_output = proc.communicate()
        if proc.returncode not in (0, 141):  # 141 = SIGPIPE when consumer exits early
            raise Exception(f"{description} command failed (exit {proc.returncode}): {cmd}\nSTDERR: {stderr_output}")
