# MakroLyzer

Domain terminology for polymer structures and their construction.

## Language

**Monolignol**:
A lignin precursor: p-coumaryl alcohol (H), coniferyl alcohol (G), or sinapyl alcohol (S).

**Lignin unit**:
The monolignol-derived portion of a lignin structure, whose bonding and functional groups can differ from those of the free precursor.

**Interunit linkage**:
A covalent connection between specified positions on two lignin units, such as beta-O-4 or beta-1. A structural motif can contain several interunit linkages.

**Beta-1 linkage**:
A carbon-carbon interunit connection between the beta carbon of one lignin unit and ring carbon 1 of another. In native-lignin motif labels, beta-1 also identifies the characteristic coupling that forms spirodienone; see [Abu-Omar et al., Figure 1](https://pubs.rsc.org/en/content/articlehtml/2021/ee/d0ee02870c).

**Conventional acyclic beta-1 structure**:
A 1,2-diarylpropane-1,3-diol-type lignin structure containing a beta-1 connection without the alpha-O-alpha closure of spirodienone. This form is readily generated from spirodienone under acidic treatment; see [Imamura et al., 2024, introduction](https://doi.org/10.1515/hf-2023-0100).

**Structural motif**:
A named chemical substructure spanning lignin units and their interunit linkages, such as resinol, phenylcoumaran, or spirodienone. Different motifs can share lignin units.

**Reactive site**:
An atom or functional group that can participate in a motif-forming reaction in the current lignin structure. Its suitability depends on the other groups and existing connections required by that reaction.

**Spirodienone motif**:
The local core formed by two lignin units with beta-1/alpha-O-alpha connections and a dienone ring, commonly labeled spirodienone (beta-1). Neighboring structure can be omitted using wavy continuation bonds, as in [Imamura et al., Figure 1](https://www.degruyterbrill.com/document/doi/10.1515/hf-2023-0100/asset/graphic/j_hf-2023-0100_fig_001.jpg).

**Spirodienone trimer model**:
An extended representation of the spirodienone core with an additional beta-O-4-linked aryl partner, as in [Elder et al., Figure 1](https://research.fs.usda.gov/download/treesearch/54786.pdf). Elder's third aryl fragment is truncated relative to an intact phenylpropanoid unit.

**Spirodienone-beta-1**:
The MakroLyzer name for a trimer containing the spirodienone core and its beta-O-4-linked partner, with three complete H/G/S-derived units.

**Lignin ensemble**:
A collection of lignin molecules representing a sample's distribution of molecular structures and molar masses. Its composition and linkage statistics describe the collection as a whole.

**Lignin strand**:
A single connected lignin molecule, which may be branched. Its branch arms are parts of the same strand rather than additional strands.

**Total lignin-unit count**:
The number of monolignol-derived units across all strands in an ensemble, counting each unit once, including units belonging to composite motifs.

**Lignin project**:
An editable MakroLyzer lignin structure together with its lignin-unit identities, atom connectivity, coordinates, motif assignments, and construction metadata.

**Number-average molar mass, Mn**:
The arithmetic mean of molecular molar masses when each molecule has equal number weight, expressed in g/mol. For a collection containing one instance of each generated molecule, Mn = sum(Mi) / N.

**Weight-average molar mass, Mw**:
The mean molar mass weighted by each species' contribution to the sample mass, expressed in g/mol. For a collection containing one instance of each generated molecule, Mw = sum(Mi**2) / sum(Mi).

**Molar-mass dispersity, Đ**:
The dimensionless ratio Mw/Mn describing the spread of a sample's molar-mass distribution; see [IUPAC Recommendations 2009](https://publications.iupac.org/pac/pdf/2009/pdf/8102x0351.pdf).
_Avoid_: PDI as a separate independent material property.
