import pytest

from src.MakroLyzer.input_handling import readXYZ
from src.MakroLyzer.input_handling import readXYZ
from src.MakroLyzer.structure_modules import graphs
from src.MakroLyzer.structure_modules.endToEndDistance import end_to_end_dist
from src.MakroLyzer.structure_modules.dihedrals import get_all_dihedrals, get_CisTrans
from src.MakroLyzer.structure_modules.radiusOfGyration import get_radius_of_gyration
from src.MakroLyzer.structure_modules.anisotropy import get_anisotropy_factor
from src.MakroLyzer.structure_modules.asphericityParameter import get_asphericity_parameter
from src.MakroLyzer.structure_modules.hbonds import get_Hbonds
from src.MakroLyzer.structure_modules.countSubgraphs import count_subgraphs, count_rings
from src.MakroLyzer.structure_modules.orderParameter import get_S_from_Q

@pytest.fixture
def sample_data1():
    return 'test_structures/01.xyz'

@pytest.fixture
def sample_data2():
    return 'test_structures/02.xyz'

@pytest.fixture
def sample_data3():
    return 'test_structures/03.xyz'

@pytest.fixture
def sample_data4():
    return 'test_structures/04.xyz'

@pytest.fixture
def pattern_data04():
    return 'test_structures/04pattern.txt'

@pytest.fixture
def sample_data5():
    return 'test_structures/05.xyz'

@pytest.fixture
def sample_data5cis():
    return 'test_structures/05cis.xyz'

@pytest.fixture
def sample_data5cis2():
    return 'test_structures/05cis2.xyz'

@pytest.fixture
def sample_data5trans():
    return 'test_structures/05trans.xyz'

@pytest.fixture
def sample_data6():
    return 'test_structures/06.xyz'

@pytest.fixture
def sample_data7():
    return 'test_structures/07.xyz'

@pytest.fixture
def sample_data8():
    return 'test_structures/Hbonds/01.xyz'

@pytest.fixture
def sample_data9():
    return 'test_structures/Hbonds/02.xyz'

@pytest.fixture
def sample_data10():
    return 'test_structures/Hbonds/03.xyz'

@pytest.fixture
def sample_data11():
    return 'test_structures/Hbonds/04.xyz'

@pytest.fixture
def sample_data12():        
    return 'test_structures/Hbonds/05.xyz'

@pytest.fixture
def sample_data13():
    return 'test_structures/13.xyz'

@pytest.fixture
def sample_data14():
    return 'test_structures/14.xyz'

@pytest.fixture
def sample_data15():
    return 'test_structures/15.xyz'

@pytest.fixture
def sample_data16():
    return 'test_structures/Hbonds/06.xyz'

@pytest.fixture
def sample_data18():
    return 'test_structures/18.xyz'

@pytest.fixture
def sample_data19():
    return 'test_structures/19.xyz'

@pytest.fixture
def sample_data20():
    return 'test_structures/20.xyz'

@pytest.fixture
def sample_data21():
    return 'test_structures/21.xyz'

@pytest.fixture
def sample_data22():
    return 'test_structures/22.xyz'

@pytest.fixture
def sample_data23():
    return 'test_structures/23.xyz'

@pytest.fixture
def sample_data25():
    return 'test_structures/25.xyz'

@pytest.fixture
def sample_data26():
    return 'test_structures/26.xyz'

@pytest.fixture
def sample_data30():
    return 'test_structures/30.xyz'

@pytest.fixture
def sample_data31():
    return 'test_structures/31.xyz'

@pytest.fixture
def sample_data32():
    return 'test_structures/32.xyz'
    
def test_end_to_end(sample_data1):
    xyz = next(readXYZ.readXYZ(sample_data1))
    testGraph = graphs.GraphManager(xyz)
    distance = end_to_end_dist(testGraph)
    correctDistance = 5.595 
    assert distance[0] == pytest.approx(correctDistance, abs=1e-3)
    
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
        
    
# Radius of gyration tests
def test_radius_of_gyration(sample_data7):
    xyz = next(readXYZ.readXYZ(sample_data7))
    testGraph = graphs.GraphManager(xyz)
    
    print(testGraph.size())
    
    # test center of mass
    com = testGraph.get_com()
    assert com[0] == pytest.approx(-0.6814, abs=1e-3) 
    assert com[1] == pytest.approx(-0.8874, abs=1e-3)
    assert com[2] == pytest.approx(0.0000, abs=1e-3)
    
    # Calculate radius of gyration
    Rg_subgraphs, R_whole = get_radius_of_gyration(testGraph)
    
    # Check radius of gyration for the whole graph
    assert R_whole == pytest.approx(3.6197, abs=1e-3)  
    
    # Check radius of gyration for subgraphs
    assert len(Rg_subgraphs) == 4 
    for R in Rg_subgraphs:
        assert R == pytest.approx(0.0, abs=1e-4)
        
    
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
            