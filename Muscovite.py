# Projeto de programa para cálculo e projeção de química mineral
# V 0.0.1

########################################################################################################################
#                                            Muscovite Elements Normalization                                          #
########################################################################################################################
import pandas as pd


class Muscovites:
    def __init__(self, dataframe):
        self.epma_table = dataframe
        self.temp_table = pd.DataFrame()
        self.norm_df = pd.DataFrame()

    def normalization(self, all_ferrous=True, n_oxygens=22):
        """Anions proportion calculation
            :param n_oxygens: Number of oxygens. Default=22
            :param all_ferrous: True for all Fe2+ or False for aff Fe3+ calculation
        """

        # if all_ferrous:
        self.temp_table['Fe2+'] = (self.epma_table['FeO'] / 71.8464)
        self.temp_table['Si'] = (self.epma_table['SiO2'] / 60.0843) * 2
        self.temp_table['Ti'] = (self.epma_table['TiO2'] / 79.8788) * 2
        self.temp_table['Al'] = (self.epma_table['Al2O3'] / 101.961) * 3
        self.temp_table['Mn'] = self.epma_table['MnO'] / 70.9374
        self.temp_table['Mg'] = self.epma_table['MgO'] / 40.3044
        self.temp_table['Ca'] = self.epma_table['CaO'] / 56.0794
        self.temp_table['Na'] = self.epma_table['Na2O'] / 61.9774
        self.temp_table['K'] = self.epma_table['K2O'] / 94.196
        self.temp_table['sum_oxygen_units'] = self.temp_table.sum(axis='columns')

        major_elements = ['Si', 'Ti', 'Al', 'Fe2+', 'Mn', 'Mg', 'Ca', 'Na', 'K']

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

        # major_elements = ['Si', 'Ti', 'Al', 'Fe3+', 'Mn', 'Mg', 'Ca', 'Na', 'K']

        for index, value in enumerate(self.temp_table.index):
            for column in major_elements:
                self.temp_table.loc[index, column] = self.temp_table[column].loc[index] / \
                                                     self.temp_table['sum_oxygen_units'].loc[index] * n_oxygens

        self.temp_table['Si'] /= 2
        self.temp_table['Ti'] /= 2
        self.temp_table['Al'] *= (2 / 3)
        self.temp_table['Na'] *= 2
        self.temp_table['K'] *= 2

        self.norm_df = self.temp_table[major_elements]

        return self.norm_df

    def tables(self):
        tables = [self.epma_table, self.norm_df]
        concat_tables = pd.concat(tables, axis=1)
        return concat_tables
        pass


########################################################################################################################
#                                           Normalização de elementos testes                                           #
########################################################################################################################


pd.set_option('display.max_columns', None)
df = pd.read_csv("example_table_muscovite.csv")
# df.drop(columns=['Unnamed: 0'], inplace=True)
# df.to_csv("example_table_muscovite.csv", index=False)

# print(df.head(50))
muscovitas = Muscovites(df)
df2 = muscovitas.normalization()
df3 = muscovitas.tables()
# print(df2.head(10))
print(df3.set_index('Group').T)
