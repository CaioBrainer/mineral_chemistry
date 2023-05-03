# Projeto de programa para cálculo e projeção de química mineral
# V 0.0.1

########################################################################################################################
#                                            Biotite Elements Normalization                                            #
########################################################################################################################
import pandas as pd


def normalization(epma_table, all_ferric=True):
    """Anions proportion calculation
        :param all_ferric: True for al Fe2+ or False for Fe3+ calculation
        :param epma_table: Biotite EPMA tabular data.
    """

    df = pd.DataFrame()
    if all_ferric:
        df['Si'] = epma_table['SiO2'] / (60.0843/2)
        df['Ti'] = epma_table['TiO2'] / (79.8788/2)
        df['Al'] = epma_table['Al2O3'] / 33.9867
        df['Fe2+'] = epma_table['FeO'] / 71.8464
        df['Mn'] = epma_table['MnO'] / 70.9374
        df['Mg'] = epma_table['MgO'] / 40.3044
        df['Ca'] = epma_table['CaO'] / 56.0794
        df['Na'] = epma_table['Na2O'] / 61.9774
        df['K'] = epma_table['K2O'] / 94.196
        df['total_anion'] = df.sum(axis='columns')
        