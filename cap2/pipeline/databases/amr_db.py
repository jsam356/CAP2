
import luigi
from os.path import join, basename, dirname, abspath
import subprocess

from ..config import PipelineConfig
from ..utils.conda import CondaPackage


class GrootDB(luigi.Task):
    config_filename = luigi.Parameter()
    cores = luigi.IntParameter(default=1)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pkg = CondaPackage(
            package="groot",
            executable="groot",
            channel="bioconda"
        )
        self.config = PipelineConfig(self.config_filename)
        self.db_dir = self.config.db_dir
        self.msas = None

    @property
    def groot_index(self):
        return self.output()['groot_index'].path

    def output(self):
        groot_index = luigi.LocalTarget(join(self.db_dir, 'index_groot-db.90'))
        groot_index.makedirs()
        return {
            'groot_index': groot_index,
        }

    def get_msas(self):
        if self.msas is not None:
            return self.msas
        dl_cmd = f'{self.pkg.bin} get -o {self.db_dir} -d groot-db'
        subprocess.check_call(dl_cmd, shell=True)
        self.msas = join(self.db_dir, 'groot-db.90')
        return self.msas

    def run(self):
        groot_index = self.groot_index  # call first so directories get made
        index_cmd = f'cd {dirname(self.groot_index)} && '
        index_cmd += f'{abspath(self.pkg.bin)} index -i {self.get_msas()} -o ./{basename(self.groot_index)} '
        index_cmd += f' --logFile /dev/null -l 250 --containment -j 0.5'
        subprocess.check_call(index_cmd, shell=True)


class MegaResDB(luigi.Task):

    config_filename = luigi.Parameter()
    cores = luigi.IntParameter(default=1)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pkg = CondaPackage(
            package="mash",
            executable="mash sketch",
            channel="bioconda"
        )
        self.config = PipelineConfig(self.config_filename)
        self.db_dir = self.config.db_dir
        self.fastqs = []

    @property
    def bowtie2_index(self):
        return 'hmp_mash_sketch.msh'

    @property
    def fasta(self):
        return self._fasta

    @property
    def annotations(self):
        return self._annotations

    def output(self):
        sketch = luigi.LocalTarget(join(self.db_dir, self.mash_sketch))
        sketch.makedirs()
        return {
            'hmp_sketch': sketch,
        }

    def run(self):
        pass


class CardDB(luigi.Task):

    config_filename = luigi.Parameter()
    cores = luigi.IntParameter(default=1)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pkg = CondaPackage(
            package="mash",
            executable="mash sketch",
            channel="bioconda"
        )
        self.config = PipelineConfig(self.config_filename)
        self.db_dir = self.config.db_dir
        self.fastqs = []

    @property
    def mash_sketch(self):
        return 'hmp_mash_sketch.msh'

    def output(self):
        sketch = luigi.LocalTarget(join(self.db_dir, self.mash_sketch))
        sketch.makedirs()
        return {
            'hmp_sketch': sketch,
        }

    def run(self):
        pass
