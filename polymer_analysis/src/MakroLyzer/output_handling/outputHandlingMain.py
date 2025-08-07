import numpy as np
import pandas as pd

def main(results):
    """
    Main function to handle the output of the PolyLyzer analysis.

    Parameters:
    results (dict): Dictionary containing the results of the analysis.
    """
    
    # Chemical Formulas #
    # Per frame we have the chemical formulas of each subgraph together with its count.
    if 'formulas' in results and results['formulas']:
        formulas_file = results['formulas_file']
        with open(formulas_file, 'w') as file:
            file.write("Frame,Polymer,Count\n")
            for frame, subgraph_formulas in enumerate(results['formulas']):
                for subgraph, count in subgraph_formulas:
                    file.write(f"{frame},{subgraph},{count}\n")
                    
    # Number of Subgraphs #
    # Per frame we have the number of subgraphs in the polymer.
    if 'noSub' in results and results['noSub']:
        noSub_file = results['noSub_file']
        with open(noSub_file, 'w') as file:
            file.write("Frame, Number of Subgraphs\n")
            for frame, no_sub in enumerate(results['noSub']):
                file.write(f"{frame},{no_sub}\n")

    
    # End-to-End Distances #
    # Per frame we have the end-to-end distance of the subgraphs.
    if 'distances' in results and results['distances']:
        distances_file = results['distances_file']
        with open(distances_file, 'w') as file:
            file.write("Frame, End-to-End Distance (subgraphs) / Å\n")
            for frame, subgraph_distances in enumerate(results['distances']):
                subgraph_str = ','.join(f"{val:.3f}" for val in subgraph_distances)
                file.write(f"{frame},{subgraph_str}\n")
                
                
    # Dihedrals #
    # Per frame we have the diheadral counts for each angle like this: dihedrals = sorted(dihedrals, key=lambda x: x[0])
    # The angles can either be in the range of 0-180 or -180-180 degrees.
    if 'dihedrals' in results and results['dihedrals']:
        dihedrals_file = results['dihedrals_file']

        # Build a DataFrame for each frame
        dfs = []
        for i, frame_data in enumerate(results['dihedrals']):
            df = pd.DataFrame(frame_data, columns=["Angle", f"Frame {i}"])
            dfs.append(df)

        # Merge all DataFrames on the "Angle" column
        df_merged = dfs[0]
        for df in dfs[1:]:
            df_merged = pd.merge(df_merged, df, on="Angle", how="outer")

        # Sort by angle
        df_merged = df_merged.sort_values("Angle").fillna(0).astype({col: int for col in df_merged.columns if col != "Angle"})

        # Save to CSV
        df_merged.to_csv(dihedrals_file, index=False)
            
                
    # Cis-Trans Counts #
    # Per frame we have the overall cis-trans counts like this: [('Cis', cis), ('Trans', trans)]
    if 'cisTrans' in results and results['cisTrans']:
        cisTrans_file = results['cisTrans_file']
        with open(cisTrans_file, 'w') as file:
            file.write("Frame,Cis count,Trans count\n")
            for frame, ct in enumerate(results['cisTrans']):
                # Convert to a dict for safe access
                ct_dict = dict(ct)
                cis_count = ct_dict.get('Cis', 0)
                trans_count = ct_dict.get('Trans', 0)
                file.write(f"{frame},{cis_count},{trans_count}\n")
                
    
    # Radius of Gyration #
    # Per frame we have the radius of gyration for each subgraph and the overall radius of gyration.
    if 'Rg' in results and results['Rg']:
        rg_file = results['Rg_file']
        with open(rg_file, 'w') as file:
            file.write("Frame, Rg(whole graph) / Å, Rg(subgraphs) / Å\n")
            for frame, (subgraph_rgs, whole_rg) in enumerate(results['Rg']):
                subgraph_str = ','.join(f"{val:.3f}" for val in subgraph_rgs)
                file.write(f"{frame},{whole_rg:.3f},{subgraph_str}\n")
                
                
    # Hydrogen Bonds #
    # Per frame we have the element type, distance, and number of hydrogen bonds.
    if 'hbonds' in results and results['hbonds']:
        hbonds_file = results['hbonds_file']
        with open(hbonds_file, 'w') as file:
            file.write("Frame, Element Type, H-Acceptor dist / Å, Donor-Acceptor dist / Å, Angle cutoff / °, Number of Hydrogen Bonds\n")
            for frame, hbonds in enumerate(results['hbonds']):
                for element_type, HAcceptor_dist, DonorAcceptor_dist, Angle_cut, count in hbonds:
                    file.write(f"{frame},{element_type},{HAcceptor_dist:.3f},{DonorAcceptor_dist:.3f},{Angle_cut:.3f},{count}\n")

    # Subgraph Coordinates #
    # If subgraph coordinates are requested, write them to a file.
    if results.get('subgraph_coords'):
        base = results['subgraph_coords_file']
        if base.endswith('.xyz'):
            base = base[:-4]
        for frame_idx, frame_subgraphs in enumerate(results['subgraph_coords']):
            out_path = f"{base}_frame_{frame_idx}.xyz"
            with open(out_path, 'w') as f:
                for sg in frame_subgraphs:
                    # number of atoms in this subgraph
                    f.write(f" {len(sg)}\n")
                    f.write("\n")
                    for element, x, y, z in sg:
                        # space‐separated element and coords
                        f.write(f"{element} {x:.3f} {y:.3f} {z:.3f}\n")
                        
    # Anisotropy Factor #
    # If anisotropy factor is requested, write it to a file.
    if results.get('anisotropy_factor'):
        anisotropy_file = results['anisotropy_file']
        with open(anisotropy_file, 'w') as file:
            file.write("Frame, Anisotropy Factor\n")
            for frame, anisotropy in enumerate(results['anisotropy_factor']):
                file.write(f"{frame},{anisotropy:.3f}\n")
                
    # Asphericity Parameter #
    # If asphericity parameter is requested, write it to a file.
    if results.get('asphericity_parameter'):
        asphericity_file = results['asphericity_file']
        with open(asphericity_file, 'w') as file:
            file.write("Frame, Asphericity Parameter\n")
            for frame, asphericity in enumerate(results['asphericity_parameter']):
                file.write(f"{frame},{asphericity:.3f}\n")
                
    # Order Parameter #
    # If order parameter is requested, write it to a file.
    if results.get('orderParameter'):
        order_file = results['orderParameter_file']
        with open(order_file, 'w') as file:
            file.write("Frame, Order Parameter S\n")
            for frame, order_param in enumerate(results['orderParameter']):
                file.write(f"{frame},{order_param:.3f}\n")

        