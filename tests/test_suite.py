import pandas as pd
from argparse import Namespace
from hamplicons.main import main


class TestClass:
	def test_version(self):
		"""
		Purpose it to remind myself to update version
		"""
		with open('pyproject.toml') as f:
			for line in f.read():
				if line.startswith("version"):
					assert line.endswith("4.2.0")
					break

	def test_defaults(self, tmp_path):
		"""
		test basic functionality

		:param tmp_path:
		:return:
		"""
		args = Namespace(
			t=8, log_level='INFO', targets_fasta='tests/data/targets.fa', data_directory='tests/data',
			data_directory_recursive=False, skip_merging=False, deactivate_relaxed_matching=False, hamming_distance=4,
			o=tmp_path, output_formats='xlsx', convert_to_wells=False, plate_layout=[8, 12]
		)
		main(args)
		files = set([x.name for x in tmp_path.glob("*")])
		assert {"hamp_out.log", "hamp_out.merged", "hamp_out.xlsx"} == files
		assert pd.read_excel(tmp_path / "hamp_out.xlsx").equals(pd.read_excel("tests/output/hamp_out.xlsx"))
