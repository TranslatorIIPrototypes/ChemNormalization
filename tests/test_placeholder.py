from src.ChemNormalization import ChemNormalization
import pandas as pd
import pytest

def test_SMILES_and_Normalization():
    # instantiate the class that does all the work
    cn = ChemNormalization()

    # alter the run parameters. turning off not output puts the app in test mode
    cn._do_KGX = 0
    cn._do_redis = 0

    # load the data frame with some simple SMILES
    df: pd.DataFrame = cn.get_simplified_smiles_for_chemicals()

    # convert the ambiguous series returned and get the name
    test = pd.Series(df.loc[df['chem_id'] == 'CHEBI:33609', 'name'] == 'elemental boron')

    # spot check the data
    assert(test.any())

    # convert the ambiguous series returned and get the simple smiles
    test = pd.Series(df.loc[df['chem_id'] == 'CHEBI:64002', 'simplified_SMILES'] == 'C=CCN1CCc2c(cc(O)c(O)c2Cl)C(c2cccc(C)c2)C1')

    # spot check the data
    assert(test.any())

    # load the data frame with some node normalized data
    df_n = cn.normalize_node_data(df)

    test = pd.Series(df_n.loc[df_n['chem_id'] == 'CHEMBL.COMPOUND:CHEMBL480049', 'name'] == 'CHEMBL480049')

    # spot check the data
    assert(test.any())
