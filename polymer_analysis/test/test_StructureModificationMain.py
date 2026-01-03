import tempfile
from pathlib import Path
from unittest.mock import Mock

import pytest

from MakroLyzer.modify_modules import modifier_registry
from MakroLyzer.modify_modules import structureModificationMain
from MakroLyzer.modify_modules.SubgraphPrint import SubgraphPrintModifier


class TestModifierRegistry:
    def test_create_subgraph_print_none_and_instance(self, tmp_path):
        # No flag -> factory returns None
        args = {}
        assert modifier_registry.create_subgraph_print(args) is None

        # With flag -> returns SubgraphPrintModifier instance
        temp_file = tmp_path / "Subgraph.xyz"
        args = {"subgraph_coords": True}
        instance = modifier_registry.create_subgraph_print(args)
        assert isinstance(instance, SubgraphPrintModifier)


class TestStructureModificationMain:
    def test_main_uses_registry_and_runs_modifiers(self, monkeypatch, tmp_path):
        # Prepare args expected by structureModificationMain
        xyz_path = str(tmp_path / "traj.xyz")

        args = {
            'xyzFile': xyz_path,
            'lmpFile': None,
            'nthStep': 1,
            'BoxSize': None,
            # enable a test modifier flag used in the patched registry
            'test_flag': True,
        }

        # Create a mock modifier with initialize_output and run
        modifier_mock = Mock()
        modifier_mock.initialize_output = Mock()
        modifier_mock.run = Mock()

        # Factory returns our mock modifier
        def factory(a, **ctx):
            return modifier_mock

        # Patch the registry used by structureModificationMain to only contain our factory
        monkeypatch.setattr(structureModificationMain, 'MODIFIERS_REGISTRATION', {'test_modifier': factory})

        # Patch readXYZ and estimateFrames to yield two frames
        def fake_read(path):
            yield 'frame1'
            yield 'frame2'

        monkeypatch.setattr(structureModificationMain, 'read', None, raising=False)
        monkeypatch.setattr('MakroLyzer.modify_modules.structureModificationMain.readXYZ.readXYZ', lambda p: fake_read(p))
        monkeypatch.setattr('MakroLyzer.modify_modules.structureModificationMain.estimateFrames.EstimateFrames.estimateFramesXYZ', lambda p: 2)

        # Patch GraphManager to accept the frame and return a dummy graph
        monkeypatch.setattr('MakroLyzer.modify_modules.structureModificationMain.graphs.GraphManager', lambda frame, boxSize=None: Mock())

        # Run main
        structureModificationMain.main(args)

        # validate factory produced modifier and lifecycle methods were called
        modifier_mock.initialize_output.assert_called()
        # run should be called once per yielded frame
        assert modifier_mock.run.call_count == 2
