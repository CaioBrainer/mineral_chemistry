# Projeto de programa para cálculo e projeção de química mineral
# V 0.0.1

########################################################################################################################
#                                            Biotite Elements Normalization                                            #
########################################################################################################################
import pandas as pd


class Biotites:

    def __init__(self, dataframe):
        self.epma_table = dataframe
        pass


def normalization(epma_table, all_ferrous=True, n_oxygens=22):
    """Anions proportion calculation
        :param n_oxygens: Number of oxygens
        :param all_ferrous: True for al Fe2+ or False for Fe3+ calculation
        :param epma_table: Biotite EPMA tabular data.
    """

    df = pd.DataFrame()
    # if all_ferrous:
    df['Fe2+'] = (epma_table['FeO'] / 71.8464)
    df['Si'] = (epma_table['SiO2'] / 60.0843) * 2
    df['Ti'] = (epma_table['TiO2'] / 79.8788) * 2
    df['Al'] = (epma_table['Al2O3'] / 101.961) * 3
    df['Mn'] = epma_table['MnO'] / 70.9374
    df['Mg'] = epma_table['MgO'] / 40.3044
    df['Ca'] = epma_table['CaO'] / 56.0794
    df['Na'] = epma_table['Na2O'] / 61.9774
    df['K'] = epma_table['K2O'] / 94.196
    df['sum_oxygen_units'] = df.sum(axis='columns')
    print(df['sum_oxygen_units'])
    # else:
    #     df['Fe3+'] = (epma_table['FeO'] * 159.688 / (2 * 71.844)) * 3
    #     df['Si'] = (epma_table['SiO2'] / 60.0843)
    #     df['Ti'] = (epma_table['TiO2'] / 79.8788) * 2
    #     df['Al'] = (epma_table['Al2O3'] / 101.961) * 3
    #     df['Mn'] = epma_table['MnO'] / 70.9374
    #     df['Mg'] = epma_table['MgO'] / 40.3044
    #     df['Ca'] = epma_table['CaO'] / 56.0794
    #     df['Na'] = epma_table['Na2O'] / 61.9774
    #     df['K'] = epma_table['K2O'] / 94.196
    #     df['sum_oxygen_units'] = df.sum(axis='columns')

    major_elements = ['Si', 'Ti', 'Al', 'Fe2+', 'Mn',
                      'Mg', 'Ca', 'Na', 'K']

    for index, value in enumerate(df.index):
        for column in major_elements:
            df.loc[index, column] = df[column].loc[index] / df['sum_oxygen_units'].loc[index] * n_oxygens

    df['Si'] /= 2
    df['Ti'] /= 2
    df['Al'] *= (2 / 3)
    df['Na'] *= 2
    df['K'] *= 2

    return df[major_elements]


########################################################################################################################
#                                           Normalização de elementos testes                                           #
########################################################################################################################


pd.set_option('display.max_columns', None)
df = pd.read_csv("example_table_biotite.csv")
# df.drop(columns=['Unnamed: 0'], inplace=True)
# df.to_csv("example_table_biotite.csv", index=False)

print(df.head(50))
df2 = normalization(df)
print(df2.head(10))
