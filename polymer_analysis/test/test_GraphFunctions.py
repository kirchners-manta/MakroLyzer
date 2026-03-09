import pytest
import numpy as np

from MakroLyzer.input_handling import readXYZ
from MakroLyzer.input_handling import subgraphSelection
from MakroLyzer import graphs

@pytest.fixture
def sample_data1():
    return 'test_structures/Graph/01.xyz'

@pytest.fixture
def sample_data2():
    return 'test_structures/Graph/02.xyz'

@pytest.fixture
def sample_data3():
    return 'test_structures/Graph/03.xyz'

@pytest.fixture
def sample_data4():
    return 'test_structures/Graph/04.xyz'

@pytest.fixture
def pattern_data04():
    return 'test_structures/Graph/04pattern.txt'

@pytest.fixture
def sample_data10():
    return 'test_structures/Hbonds/03.xyz'

@pytest.fixture
def sample_data11():
    return 'test_structures/Hbonds/04.xyz'

@pytest.fixture
def sample_data14():
    return 'test_structures/Graph/14.xyz'

@pytest.fixture
def sample_data18():
    return 'test_structures/Graph/18.xyz'

@pytest.fixture
def sample_data19():
    return 'test_structures/Graph/19.xyz'

@pytest.fixture
def sample_data20():
    return 'test_structures/Graph/20.xyz'

@pytest.fixture
def sample_data21():
    return 'test_structures/Graph/21.xyz'

@pytest.fixture
def sample_data23():
    return 'test_structures/Graph/23.xyz'

@pytest.fixture
def sample_data25():
    return 'test_structures/Graph/25.xyz'

@pytest.fixture
def sample_data26():
    return 'test_structures/Graph/26.xyz'
    
def test_remove_1order(sample_data1):
    xyz = next(readXYZ.readXYZ(sample_data1))
    testGraph = graphs.GraphManager(xyz)
    
    # number of nodes
    assert testGraph.number_of_nodes() == 26
    newGraph = testGraph.remove_1order()
    assert newGraph.number_of_nodes() == 8
    
    # type of nodes
    for node in newGraph.nodes():
        assert newGraph.nodes[node]['element'] == 'C'
        
def test_update_degree(sample_data1):
    xyz = next(readXYZ.readXYZ(sample_data1))
    testGraph = graphs.GraphManager(xyz)
    testGraph.surrounding()
    
    # degree of C nodes
    for node in testGraph.nodes():
        if testGraph.nodes[node]['element'] == 'C':
            assert testGraph.nodes[node]['degree'] == 4
            
    newGraph = testGraph.remove_1order()
    newGraph.update_degree()
    
    for node in newGraph.nodes():
        if newGraph.nodes[node]['element'] == 'C' and newGraph.nodes[node]['surroundingAtoms'] == 'C_CCHH':
            assert newGraph.nodes[node]['degree'] == 2
        elif newGraph.nodes[node]['element'] == 'C' and newGraph.nodes[node]['surroundingAtoms'] == 'C_CHHH':
            assert newGraph.nodes[node]['degree'] == 1

def test_update_coordinates(sample_data1):
    xyz = next(readXYZ.readXYZ(sample_data1))
    testGraph = graphs.GraphManager(xyz)

    original = testGraph.get_coordinates(0).copy()
    updated = xyz.copy()
    updated[['x', 'y', 'z']] = updated[['x', 'y', 'z']] + 2.5

    testGraph.update_coordinates(updated)
    assert testGraph.get_coordinates(0) == pytest.approx(original + 2.5)

    bad_data = updated.copy()
    bad_data.loc[0, 'atom'] = 'O' if bad_data.loc[0, 'atom'] != 'O' else 'C'
    with pytest.raises(ValueError):
        testGraph.update_coordinates(bad_data)
        
def test_longest_path(sample_data1):
    xyz = next(readXYZ.readXYZ(sample_data1))
    testGraph = graphs.GraphManager(xyz)
    testGraph.surrounding()
    
    # get longest path 
    longestPath = testGraph.find_longest_path()
    assert len(longestPath) == 10
    path = (['H_C', 'C_CHHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CHHH', 'H_C'])
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path[i]
    
    newGraph = testGraph.remove_1order()
    newGraph.update_degree()
    longestPath = newGraph.find_longest_path()
    assert len(longestPath) == 8
    path = (['C_CHHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CHHH'])
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path[i]
        
def test_longest_path_cycle(sample_data14):
    xyz = next(readXYZ.readXYZ(sample_data14))
    testGraph = graphs.GraphManager(xyz)
    testGraph.surrounding()
    
    # get longest path 
    longestPath = testGraph.find_longest_path()
    assert len(longestPath) == 4
    path = (['C_CC', 'C_CC', 'C_CC', 'C_CC'])
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path[i]
        
    it = iter(testGraph.nodes())    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 4
        
def test_longest_path_cycle2(sample_data18):
    xyz = next(readXYZ.readXYZ(sample_data18))
    testGraph = graphs.GraphManager(xyz)
    testGraph = testGraph.remove_1order()
    testGraph.surrounding()
    
    
    path = (['C_CC', 'C_CC', 'C_CCC', 'C_CC', 'C_CCC', 'C_CC', 'C_CC', 'C_CCC', 'C_CC', 'C_CC', 'C_CCC'])
    # get longest path 
    longestPath = testGraph.find_longest_path()
    assert len(longestPath) == 11
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path[i]
        
    it = iter(testGraph.nodes())    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 11
    
    it = iter(testGraph.nodes())    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 11
        
def test_longest_path_cycle3(sample_data19):
    xyz = next(readXYZ.readXYZ(sample_data19))
    testGraph = graphs.GraphManager(xyz)
    testGraph = testGraph.remove_1order()
    testGraph.surrounding()
    
    longestPath = testGraph.find_longest_path()
    assert len(longestPath) == 28
    
def test_longest_path_cycle4(sample_data20):
    xyz = next(readXYZ.readXYZ(sample_data20))
    testGraph = graphs.GraphManager(xyz)
    testGraph = testGraph.remove_1order()
    testGraph.surrounding()
    
    it = iter(testGraph.nodes())
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 42
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 42
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 42
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 42
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 42
    
def test_longest_path_cycle5(sample_data25):
    xyz = next(readXYZ.readXYZ(sample_data25))
    testGraph = graphs.GraphManager(xyz)
    testGraph = testGraph.remove_1order()
    testGraph.surrounding()
    
    it = iter(testGraph.nodes())
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 12
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 12
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 12
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 12
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 12    
    
def test_longest_path_cycle6(sample_data26):
    xyz = next(readXYZ.readXYZ(sample_data26))
    testGraph = graphs.GraphManager(xyz)
    testGraph = testGraph.remove_1order()
    testGraph.surrounding()
    
    path1 = (['C_CCC', 'C_CC', 'C_CC', 'C_CC', 'C_CC', 'C_CC'])
    path2 = (['C_CC', 'C_CCC', 'C_CC', 'C_CC', 'C_CC', 'C_CC'])
    path3 = (['C_CC', 'C_CC', 'C_CCC', 'C_CC', 'C_CC', 'C_CC'])
    path4 = (['C_CC', 'C_CC', 'C_CC', 'C_CCC', 'C_CC', 'C_CC'])
    path5 = (['C_CC', 'C_CC', 'C_CC', 'C_CC', 'C_CCC', 'C_CC'])
    path6 = (['C_CC', 'C_CC', 'C_CC', 'C_CC', 'C_CC', 'C_CCC'])
    
    it = iter(testGraph.nodes())
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 6
    
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path5[i]  
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 6
    
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path6[i]  
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 6
    
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path1[i]    

    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 6 
    
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path2[i]  
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 6     
    
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path3[i]  
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 6
    
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path2[i]  
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 6
    
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path6[i]  
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 6
    
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path6[i]  
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 6
    
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path6[i]  
    
    startAtom = next(it)
    longestPath = testGraph.find_longest_path(startAtom=startAtom)
    assert len(longestPath) == 6 
    
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path6[i]  
    

def test_surrounding(sample_data2):
    xyz = next(readXYZ.readXYZ(sample_data2))
    testGraph = graphs.GraphManager(xyz)
    testGraph.surrounding()
    
    # check surrounding atoms
    for node in testGraph.nodes():
        if testGraph.nodes[node]['element'] == 'C':
            assert testGraph.nodes[node]['surroundingAtoms'] == 'C_MgNOP'
        elif testGraph.nodes[node]['element'] == 'O':
            assert testGraph.nodes[node]['surroundingAtoms'] == 'O_CH'
        elif testGraph.nodes[node]['element'] == 'N':
            assert testGraph.nodes[node]['surroundingAtoms'] == 'N_CHH'
        elif testGraph.nodes[node]['element'] == 'P':
            assert testGraph.nodes[node]['surroundingAtoms'] == 'P_CHH'
        elif testGraph.nodes[node]['element'] == 'Mg':
            assert testGraph.nodes[node]['surroundingAtoms'] == 'Mg_CH'
        

def test_find_patterns(sample_data4):
    xyz = next(readXYZ.readXYZ(sample_data4))
    testGraph = graphs.GraphManager(xyz)
    testGraph.find_and_tag_patterns([['C_CCC', 'C_CCC', 'C_CC', 'C_CC', 'C_CC', 'C_CC', 'C_CC', 'C_CC'], ['C_CC', 'C_CCC', 'C_CC', 'C_CC', 'C_CC', 'C_CC', 'C_CC'], ['C_C']])
    
    # iterate over nodes and check fragment ID
    i = 0
    for node in testGraph.nodes():
        if i in range(0,4):
            assert testGraph.nodes[node]['fragmentID'] == 4
        elif i in range(4, 20):
            assert testGraph.nodes[node]['fragmentID'] == 0
        elif i in range(20, 36):
            assert testGraph.nodes[node]['fragmentID'] == 1
        elif i in range(36, 52):
            assert testGraph.nodes[node]['fragmentID'] == 2
        elif i in range(52, 66):
            assert testGraph.nodes[node]['fragmentID'] == 3
        i += 1
    
# PBC
def test_pbc(sample_data21):
    xyz = next(readXYZ.readXYZ(sample_data21))
    testGraph = graphs.GraphManager(xyz, boxSize=12)
    
    refformula = [('C8H18', 1)]
    formula = testGraph.get_chemicalFormulas()
    
    for i in range(len(refformula)):
        assert formula[i][0] == refformula[i][0]
        assert formula[i][1] == refformula[i][1]
        
# get_vectors_and_positions_along_path
def test_get_vectors_and_positions_along_path(sample_data23):
    xyz = next(readXYZ.readXYZ(sample_data23))
    testGraph = graphs.GraphManager(xyz)
    
    # remove 1-order nodes
    testGraph = testGraph.remove_1order()
    testGraph.update_degree()
    
    # get subgraphs
    subgraphs = testGraph.get_subgraphs()
    
    for subgraph in subgraphs:
        longestPath = subgraph.find_longest_path()
        if len(longestPath) < 2:
            continue
        
        # get vectors and positions along path
        pathDict = subgraph.get_vectors_and_positions_along_path(longestPath, unitSize=2)
        print(pathDict)
        
        assert len(pathDict) == 2
        
        position = list(pathDict.keys())
        vector = list(pathDict.values())
        
        assert len(vector) == 2
        assert vector[0][0][0] == pytest.approx(0.911446, abs=1e-5)
        assert vector[0][0][1] == pytest.approx(0.410612, abs=1e-5)
        assert vector[0][0][2] == pytest.approx(-0.025753, abs=1e-5)
        
        assert len(position) == 2
        assert position[0][0] == pytest.approx(-3.023285, abs=1e-5)
        assert position[0][1] == pytest.approx(-0.12005, abs=1e-5)
        assert position[0][2] == pytest.approx(-0.021475, abs=1e-5)
        
    
# test AminoAcid Backbone search
def test_AAbackbone(sample_data10):
    xyz = next(readXYZ.readXYZ(sample_data10))
    testGraph = graphs.GraphManager(xyz)

    backbone = testGraph.AminoAcidBackbone()
    assert backbone is not None
    assert len(backbone) == 37

    for i, atom in enumerate(backbone):
        if i % 3 == 0:
            assert testGraph.nodes[backbone[i]]['element'] == 'N'
        else:
            assert testGraph.nodes[backbone[i]]['element'] == 'C'

def test_AAbackbone_2(sample_data11):
    xyz = next(readXYZ.readXYZ(sample_data11))
    testGraph = graphs.GraphManager(xyz)

    backbone = testGraph.AminoAcidBackbone()
    assert backbone is not None
    assert len(backbone) == 37

    for i, atom in enumerate(backbone):
        if i % 3 == 0:
            assert testGraph.nodes[backbone[i]]['element'] == 'N'
        else:
            assert testGraph.nodes[backbone[i]]['element'] == 'C'


def test_functionalize_pe_co():
    graph = graphs.GraphManager()
    # Three-carbon chain, middle C has two H neighbors to functionalize
    graph.add_node(0, index=0, element='C', x=0.0, y=0.0, z=0.0)
    graph.add_node(1, index=1, element='C', x=1.54, y=0.0, z=0.0)
    graph.add_node(2, index=2, element='C', x=3.08, y=0.0, z=0.0)
    graph.add_node(3, index=3, element='H', x=1.54, y=1.0, z=0.0)
    graph.add_node(4, index=4, element='H', x=2.0, y=-0.5, z=0.2)

    graph.add_edge(0, 1)
    graph.add_edge(1, 2)
    graph.add_edge(1, 3)
    graph.add_edge(1, 4)

    before_nodes = graph.number_of_nodes()
    before_h = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'H')

    graph.functionalizePE("random", "CO", percentage=100, neighbor_exclusion=1)

    after_nodes = graph.number_of_nodes()
    after_h = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'H')
    o_nodes = [n for n in graph.nodes() if graph.nodes[n]['element'] == 'O']

    assert after_nodes == before_nodes - 1
    assert after_h == before_h - 2
    assert len(o_nodes) == 1
    assert o_nodes[0] in graph.neighbors(1)


def _build_linear_pe_chain(num_c):
    graph = graphs.GraphManager()
    h_index = num_c
    for i in range(num_c):
        x = 1.54 * i
        graph.add_node(i, index=i, element='C', x=x, y=0.0, z=0.0)
        if i > 0:
            graph.add_edge(i - 1, i)
        # Two hydrogens per carbon for functionalization geometry.
        graph.add_node(h_index, index=h_index, element='H', x=x, y=1.0, z=0.5)
        graph.add_edge(i, h_index)
        h_index += 1
        graph.add_node(h_index, index=h_index, element='H', x=x, y=-1.0, z=0.0)
        graph.add_edge(i, h_index)
        h_index += 1
    return graph


def test_functionalize_pe_ester_random():
    np.random.seed(0)
    graph = _build_linear_pe_chain(6)
    before_c = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'C')
    before_o = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'O')
    before_n = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'N')

    graph.functionalizePE("random", "ester", percentage=100, neighbor_exclusion=1)

    after_c = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'C')
    after_o = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'O')
    after_n = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'N')

    assert after_c == before_c - 1
    assert after_o == before_o + 2
    assert after_n == before_n


def test_functionalize_pe_amide_random():
    np.random.seed(0)
    graph = _build_linear_pe_chain(6)
    before_c = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'C')
    before_o = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'O')
    before_n = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'N')

    graph.functionalizePE("random", "amide", percentage=100, neighbor_exclusion=1)

    after_c = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'C')
    after_o = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'O')
    after_n = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'N')

    assert after_c == before_c - 1
    assert after_o == before_o + 1
    assert after_n == before_n + 1


def test_functionalize_pe_co_periodic():
    graph = _build_linear_pe_chain(6)
    before_o = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'O')
    before_h = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'H')

    graph.functionalizePE("periodic", "CO", distance=2)

    after_o = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'O')
    after_h = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'H')

    assert after_o == before_o + 2
    assert after_h == before_h - 4


def test_parse_selection_list():
    selections = subgraphSelection.parse_selection_list(['C2H6', 'CHO'])
    assert selections[0]['mode'] == 'exact'
    assert selections[0]['counts'] == {'C': 2, 'H': 6}
    assert selections[1]['mode'] == 'elements'
    assert selections[1]['elements'] == {'C', 'H', 'O'}

def test_get_mass():
    graph = graphs.GraphManager()
    graph.add_node(0, element='C', x=0.0, y=0.0, z=0.0)
    graph.add_node(1, element='H', x=1.0, y=0.0, z=0.0)
    graph.add_node(2, element='O', x=0.0, y=1.0, z=0.0)

    # C + H + O
    assert graph.get_mass() == pytest.approx(12.011 + 1.008 + 15.999)


def test_select_subgraph_nodes():
    graph = graphs.GraphManager()
    # Subgraph 1: C2H6
    graph.add_node(0, element='C')
    graph.add_node(1, element='C')
    graph.add_nodes_from([(2, {'element': 'H'}), (3, {'element': 'H'}), (4, {'element': 'H'}),
                          (5, {'element': 'H'}), (6, {'element': 'H'}), (7, {'element': 'H'})])
    graph.add_edge(0, 1)
    for h in [2, 3, 4]:
        graph.add_edge(0, h)
    for h in [5, 6, 7]:
        graph.add_edge(1, h)

    # Subgraph 2: C1H4
    graph.add_node(8, element='C')
    graph.add_nodes_from([(9, {'element': 'H'}), (10, {'element': 'H'}),
                          (11, {'element': 'H'}), (12, {'element': 'H'})])
    for h in [9, 10, 11, 12]:
        graph.add_edge(8, h)

    # Subgraph 3: CHO (element-only selection)
    graph.add_node(13, element='C')
    graph.add_node(14, element='O')
    graph.add_node(15, element='H')
    graph.add_edge(13, 14)
    graph.add_edge(14, 15)

    # Subgraph 4: N2 (should be excluded)
    graph.add_node(16, element='N')
    graph.add_node(17, element='N')
    graph.add_edge(16, 17)

    selections = subgraphSelection.parse_selection_list(['C2H6', 'C1H4', 'CHO'])
    selected_nodes = graph.select_subgraph_nodes(selections)
    assert selected_nodes == set(range(0, 16))
    
    # no N in set
    for n in selected_nodes:
        assert graph.nodes[n]['element'] != 'N'
