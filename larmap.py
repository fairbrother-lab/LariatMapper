import subprocess
from dataclasses import dataclass, fields

@dataclass
class larmap:

    fastq_dir: str
    read_1_a: str
    read_2_a: str
    read_1_b: str
    read_2_b: str
    output_base_name_a: str
    output_base_name_b: str
    num_cpus: int

    ref_b2index: str
    ref_fasta: str
    ref_gtf: str
    ref_5p_fasta: str
    ref_3p_b2index: str
    ref_3p_lengths: str
    ref_introns: str
    ref_repeatmasker: str

    results_filename: str

    def __init__(self, fastq_dir: str, read_1_a: str, read_2_a: str, read_1_b: str, read_2_b: str,
                 output_base_name_a: str, output_base_name_b: str, num_cpus: int, ref_b2index: str,
                 ref_fasta: str, ref_gtf: str, ref_5p_fasta: str, ref_3p_b2index: str, ref_3p_lengths: str,
                 ref_introns: str, ref_repeatmasker: str, results_filename: str):
        self.fastq_dir = fastq_dir
        self.read_1_a = read_1_a
        self.read_2_a = read_2_a
        self.read_1_b = read_1_b
        self.read_2_b = read_2_b
        self.output_base_name_a = output_base_name_a
        self.output_base_name_b = output_base_name_b
        self.num_cpus = num_cpus
        self.ref_b2index = ref_b2index
        self.ref_fasta = ref_fasta
        self.ref_gtf = ref_gtf
        self.ref_5p_fasta = ref_5p_fasta
        self.ref_3p_b2index = ref_3p_b2index
        self.ref_3p_lengths = ref_3p_lengths
        self.ref_introns = ref_introns
        self.ref_repeatmasker = ref_repeatmasker
        self.results_filename = results_filename
    
    def get_args(self):
        attributes = [field.name for field in fields(self)]
        values = [getattr(self, attr) for attr in attributes]
        return values


def map_lariats(larmap_obj):
    top_script = 'larmap_top.sh'
    print(larmap_obj.get_args())
    command = ['sh', top_script] + larmap_obj.get_args()
    subprocess.run(' '.join(command), shell=True, check=True)