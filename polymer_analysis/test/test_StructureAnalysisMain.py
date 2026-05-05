import tempfile
from pathlib import Path
from unittest.mock import Mock

import pytest

from MakroLyzer.structure_modules import analyzer_registry
from MakroLyzer.structure_modules import structureAnalysisMain
from MakroLyzer.structure_modules.Anisotropy import AnisotropyAnalyzer
from MakroLyzer.structure_modules.ConvexHullVolume import ConvexHullVolumeAnalyzer


class TestAnalyzerRegistry:
	def test_create_anisotropy_none_and_instance(self, tmp_path):
		# No flag -> factory returns None
		args = {}
		assert analyzer_registry.create_anisotropy(args) is None

		# With flag -> returns AnisotropyAnalyzer instance
		temp_file = tmp_path / "anisotropy.csv"
		args = {"anisotropyFactor": True, "anisotropy_file": temp_file}
		instance = analyzer_registry.create_anisotropy(args)
		assert isinstance(instance, AnisotropyAnalyzer)
		# initialize_output should create the header file (streaming mode)
		instance.initialize_output()
		assert temp_file.exists()
		assert temp_file.read_text() == "Frame, Anisotropy Factor (κ²)\n"

	def test_create_convex_hull_volume_none_and_instance(self, tmp_path):
		args = {}
		assert analyzer_registry.create_convex_hull_volume(args) is None

		temp_file = tmp_path / "convhull.csv"
		args = {"convexHullVolume": str(temp_file)}
		instance = analyzer_registry.create_convex_hull_volume(args)

		assert isinstance(instance, ConvexHullVolumeAnalyzer)
		instance.initialize_output()
		assert temp_file.exists()
		assert temp_file.read_text() == (
			"Frame, Convex Hull - Volume / Å³, Mass / g/mol, Density / g/cm³\n"
		)


class TestStructureAnalysisMain:
	def test_main_uses_registry_and_runs_analyzers(self, monkeypatch, tmp_path):
		# Prepare args expected by structureAnalysisMain
		xyz_path = str(tmp_path / "traj.xyz")
		temp_out = tmp_path / "out.csv"

		args = {
			'xyzFile': xyz_path,
			'lmpFile': None,
			'nthStep': 1,
			'BoxSize': None,
			'vibFactor': 1.15,
			'orderParameter': ([1.0], None, None),
			# enable a test analyzer name used in the patched registry
			'test_flag': True,
		}

		# Create a mock analyzer with initialize_output, run and finalize_output
		analyzer_mock = Mock()
		analyzer_mock.initialize_output = Mock()
		analyzer_mock.run = Mock()
		analyzer_mock.finalize_output = Mock()

		# Create a factory that returns our mock analyzer
		def factory(a, **ctx):
			return analyzer_mock

		# Patch the registry used by structureAnalysisMain to only contain our factory
		monkeypatch.setattr(structureAnalysisMain, 'ANALYZERS_REGISTRATION', {'test_analyzer': factory})

		# Patch readXYZ and estimateFrames to yield two frames
		def fake_read(path):
			yield 'frame1'
			yield 'frame2'

		monkeypatch.setattr(structureAnalysisMain, 'read', None, raising=False)
		# Instead, patch the read function used by main through the module import names
		monkeypatch.setattr('MakroLyzer.structure_modules.structureAnalysisMain.readXYZ.readXYZ', lambda p: fake_read(p))
		monkeypatch.setattr('MakroLyzer.structure_modules.structureAnalysisMain.estimateFrames.EstimateFrames.estimateFramesXYZ', lambda p: 2)

		# Patch GraphManager to accept the frame and return a dummy graph
		monkeypatch.setattr(
			'MakroLyzer.structure_modules.structureAnalysisMain.graphs.GraphManager',
			lambda frame, boxSize=None, vib_factor=None: Mock(),
		)

		# Run main
		structureAnalysisMain.main(args)

		# validate factory produced analyzer and lifecycle methods were called
		analyzer_mock.initialize_output.assert_called()
		# run should be called once per yielded frame
		assert analyzer_mock.run.call_count == 2
		analyzer_mock.finalize_output.assert_called()

	def test_main_static_topology_and_selection(self, monkeypatch, tmp_path):
		xyz_path = str(tmp_path / "traj.xyz")
		args = {
			'xyzFile': xyz_path,
			'lmpFile': None,
			'nthStep': 1,
			'BoxSize': None,
			'orderParameter': None,
			'staticTopology': True,
			'subgraphSelection': ['C2H6'],
		}

		analyzer_mock = Mock()
		analyzer_mock.initialize_output = Mock()
		analyzer_mock.run = Mock()
		analyzer_mock.finalize_output = Mock()

		def factory(a, **ctx):
			return analyzer_mock

		monkeypatch.setattr(structureAnalysisMain, 'ANALYZERS_REGISTRATION', {'test_analyzer': factory})

		def fake_read(path):
			yield 'frame1'
			yield 'frame2'

		monkeypatch.setattr('MakroLyzer.structure_modules.structureAnalysisMain.readXYZ.readXYZ', lambda p: fake_read(p))
		monkeypatch.setattr('MakroLyzer.structure_modules.structureAnalysisMain.estimateFrames.EstimateFrames.estimateFramesXYZ', lambda p: 2)

		class FakeGraphManager:
			created = []

			def __init__(self, data=None, boxSize=None, vib_factor=None):
				self.data = data
				self.boxSize = boxSize
				self.vib_factor = vib_factor
				self.update_calls = 0
				self.select_calls = 0
				FakeGraphManager.created.append(self)

			def update_coordinates(self, frame, boxSize=None):
				self.update_calls += 1

			def select_subgraph_nodes(self, selections):
				self.select_calls += 1
				return {0}

			def subgraph(self, nodes):
				return {'nodes': nodes}

		monkeypatch.setattr('MakroLyzer.structure_modules.structureAnalysisMain.graphs.GraphManager', FakeGraphManager)

		structureAnalysisMain.main(args)

		base_graph = FakeGraphManager.created[0]
		assert base_graph.update_calls == 1
		assert base_graph.select_calls == 1
		assert analyzer_mock.run.call_count == 2

	def test_full_graph_analyzer_does_not_override_global_selection(self, monkeypatch, tmp_path):
		xyz_path = str(tmp_path / "traj.xyz")
		args = {
			'xyzFile': xyz_path,
			'lmpFile': None,
			'nthStep': 1,
			'BoxSize': None,
			'vibFactor': 1.15,
			'orderParameter': None,
			'staticTopology': False,
			'subgraphSelection': ['C2H6'],
		}

		class FakeAnalyzer:
			requires_full_graph = False

			def __init__(self):
				self.graphs = []

			def initialize_output(self):
				pass

			def run(self, graph, frame_idx):
				self.graphs.append(graph)

			def finalize_output(self):
				pass

		class FakeFullGraphAnalyzer(FakeAnalyzer):
			requires_full_graph = True

		selected_analyzer = FakeAnalyzer()
		full_graph_analyzer = FakeFullGraphAnalyzer()
		full_graph_analyzer.requires_full_graph = True

		def selected_factory(a, **ctx):
			return selected_analyzer

		def full_graph_factory(a, **ctx):
			return full_graph_analyzer

		monkeypatch.setattr(
			structureAnalysisMain,
			'ANALYZERS_REGISTRATION',
			{
				'selected_analyzer': selected_factory,
				'full_graph_analyzer': full_graph_factory,
			},
		)

		monkeypatch.setattr('MakroLyzer.structure_modules.structureAnalysisMain.readXYZ.readXYZ', lambda p: iter(['frame1']))
		monkeypatch.setattr('MakroLyzer.structure_modules.structureAnalysisMain.estimateFrames.EstimateFrames.estimateFramesXYZ', lambda p: 1)

		class FakeGraphManager:
			def __init__(self, data=None, boxSize=None, vib_factor=None):
				self.data = data

			def select_subgraph_nodes(self, selections):
				return {0}

			def subgraph(self, nodes):
				return {'selected_nodes': nodes}

		monkeypatch.setattr('MakroLyzer.structure_modules.structureAnalysisMain.graphs.GraphManager', FakeGraphManager)

		structureAnalysisMain.main(args)

		assert selected_analyzer.graphs[0].data == {'selected_nodes': {0}}
		assert full_graph_analyzer.graphs[0].data == 'frame1'
