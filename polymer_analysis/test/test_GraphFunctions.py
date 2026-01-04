import pytest

from MakroLyzer.input_handling import readXYZ
from MakroLyzer.input_handling import readXYZ
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

    graph.functionalizePE(100, 'CO')

    after_nodes = graph.number_of_nodes()
    after_h = sum(1 for n in graph.nodes() if graph.nodes[n]['element'] == 'H')
    o_nodes = [n for n in graph.nodes() if graph.nodes[n]['element'] == 'O']

    assert after_nodes == before_nodes - 1
    assert after_h == before_h - 2
    assert len(o_nodes) == 1
    assert o_nodes[0] in graph.neighbors(1)