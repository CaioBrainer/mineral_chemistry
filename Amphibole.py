# Projeto de programa para cálculo e projeção de química mineral
# V 0.0.1
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


########################################################################################################################
#                                            Amphibole Elements Normalization                                          #
########################################################################################################################

class Amphibole:
    def __init__(self, dataframe):
        self.epma_dataframe = dataframe
        self.temp_table = pd.DataFrame()
        self.epma_normalized = pd.DataFrame()
        
    pass
    
    def normalization(self):   
        anion_proportions = pd.DataFrame()
        cation_proportions = pd.DataFrame()

        # Anion normalized calculation
        anion_proportions['Si_anion'] = (self.epma_dataframe['SiO2'] * 2) / 60.0843
        anion_proportions['Ti_anion'] = (self.epma_dataframe['TiO2'] * 2) / 79.8788
        anion_proportions['Al_anion'] = (self.epma_dataframe['Al2O3'] * 3) / 101.9602
        anion_proportions['Fe3_anion'] = (self.epma_dataframe['Fe2O3'] * 3) / 60.0843
        anion_proportions['Fe2_anion'] = self.epma_dataframe['FeO'] / 71.8464
        anion_proportions['Mn_anion'] = self.epma_dataframe['MnO'] / 70.9374
        anion_proportions['Mg_anion'] = self.epma_dataframe['MgO'] / 40.3044
        anion_proportions['Ca_anion'] = self.epma_dataframe['CaO'] / 56.0794
        anion_proportions['Na_anion'] = self.epma_dataframe['Na2O'] / 61.9774
        anion_proportions['K_anion'] = self.epma_dataframe['K2O'] / 94.196
        anion_proportions['total_anion'] = anion_proportions.sum(axis='columns')

        # Anion normalized calculation
        cation_proportions['Si_cation'] = self.epma_dataframe['SiO2'] / 60.0843
        cation_proportions['Ti_cation'] = self.epma_dataframe['TiO2'] / 79.8788
        cation_proportions['Al_cation'] = (self.epma_dataframe['Al2O3'] * 2) / 101.9602
        cation_proportions['Fe3+_cation'] = (self.epma_dataframe['Fe2O3'] * 2) / 60.0843
        cation_proportions['Fe2+_cation'] = self.epma_dataframe['FeO'] / 71.8464
        cation_proportions['Mn_cation'] = self.epma_dataframe['MnO'] / 70.9374
        cation_proportions['Mg_cation'] = self.epma_dataframe['MgO'] / 40.3044
        cation_proportions['Ca_cation'] = self.epma_dataframe['CaO'] / 56.0794
        cation_proportions['Na_cation'] = (self.epma_dataframe['Na2O'] * 2) / 61.9774
        cation_proportions['K_cation'] = (self.epma_dataframe['K2O'] * 2) / 94.196
        cation_proportions['total_cationFM+Ca'] = cation_proportions.sum(axis='columns') - cation_proportions['Na_cation'] - cation_proportions['K_cation']
        cation_proportions['total_cationFM'] = cation_proportions.sum(axis='columns') - cation_proportions['Ca_cation'] - cation_proportions['Na_cation'] \
                                      - cation_proportions['K_cation'] - cation_proportions['total_cationFM+Ca']
    
        """S13 normalized based on anions"""
        df_s13 = pd.DataFrame()
    
        df_s13['Si_S13'] = 13 * anion_proportions['Si_anion'] / cation_proportions['total_cationFM']
        df_s13['Ti_S13'] = 13 * anion_proportions['Ti_anion'] / cation_proportions['total_cationFM']
        df_s13['Al_S13'] = 13 * anion_proportions['Al_anion'] / cation_proportions['total_cationFM']
        df_s13['Fe2+_S13'] = 13 * anion_proportions['Fe2_anion'] / cation_proportions['total_cationFM']
        df_s13['Fe3+_S13'] = 0
        df_s13['Mn_S13'] = 13 * anion_proportions['Mn_anion'] / cation_proportions['total_cationFM']
        df_s13['Mg_S13'] = 13 * anion_proportions['Mg_anion'] / cation_proportions['total_cationFM']
        df_s13['Ca_S13'] = 13 * anion_proportions['Ca_anion'] / cation_proportions['total_cationFM']
        df_s13['Na_S13'] = 13 * anion_proportions['Na_anion'] / cation_proportions['total_cationFM']
        df_s13['K_S13'] = 13 * anion_proportions['K_anion'] / cation_proportions['total_cationFM']
        df_s13['total_S13'] = df_s13.sum(axis='columns')

        """S15 normalized based on anions"""
        df_s15 = pd.DataFrame()

        df_s15['Si_S15'] = 15 * anion_proportions['Si_anion'] / cation_proportions['total_cationFM+Ca']
        df_s15['Ti_S15'] = 15 * anion_proportions['Ti_anion'] / cation_proportions['total_cationFM+Ca']
        df_s15['Al_S15'] = 15 * anion_proportions['Al_anion'] / cation_proportions['total_cationFM+Ca']
        df_s15['Fe2+_S15'] = 15 * anion_proportions['Fe2_anion'] / cation_proportions['total_cationFM+Ca']
        df_s15['Fe3+_S15'] = 0
        df_s15['Mn_S15'] = 15 * anion_proportions['Mn_anion'] / cation_proportions['total_cationFM+Ca']
        df_s15['Mg_S15'] = 15 * anion_proportions['Mg_anion'] / cation_proportions['total_cationFM+Ca']
        df_s15['Ca_S15'] = 15 * anion_proportions['Ca_anion'] / cation_proportions['total_cationFM+Ca']
        df_s15['Na_S15'] = 15 * anion_proportions['Na_anion'] / cation_proportions['total_cationFM+Ca']
        df_s15['K_S15'] = 15 * anion_proportions['K_anion'] / cation_proportions['total_cationFM+Ca']
        df_s15['total_S15'] = df_s15.sum(axis='columns')

        # Fe3+ calculation

        # S13
        df_s13['Fe3+_S13'] = 2 * (23 - df_s13['total_S13'])

        for index, value in enumerate(df_s13['Fe3+_S13']):
            if value < 0:
                df_s13['Fe3+_S13'].loc[index] = 0
            else:
                df_s13['Fe3+_S13'].loc[index] = value
        # S15
        df_s15['Fe3+_S15'] = 2 * (23 - df_s15['total_S15'])

        for index, value in enumerate(df_s15['Fe3+_S15']):
            if value < 0:
                df_s15['Fe3+_S15'].loc[index] = 0
            else:
                df_s15['Fe3+_S15'].loc[index] = value

        # Fe3+ and Fe2+ calculation S13
        for index, value in enumerate(df_s13['Fe3+_S13']):
            if value > df_s13['Fe2+_S13'].loc[index]:
                df_s13.loc[index, 'Fe3+_S13'] = df_s13['Fe2+_S13'].loc[index]
            else:
                df_s13.loc[index, 'Fe3+_S13'] = df_s13['Fe3+_S13'].loc[index]

        df_s13['Fe2+_S13'] = df_s13['Fe2+_S13'] - df_s13['Fe3+_S13']

        # Fe3+ and Fe2+ calculation S15
        for index, value in enumerate(df_s15['Fe3+_S15']):
            if value > df_s15['Fe2+_S15'].loc[index]:
                df_s15['Fe3+_S15'].loc[index] = df_s15['Fe2+_S15'].loc[index]
            else:
                df_s15['Fe3+_S15'].loc[index] = df_s15['Fe3+_S15'].loc[index]

        df_s15['Fe2+_S15'] = df_s15['Fe2+_S15'] - df_s15['Fe3+_S15']

        # Cations based on 23 Oxygens
        cat_based13 = pd.DataFrame()
        cat_based15 = pd.DataFrame()

        # S13 Cations + Recheck
        cat_based13['Si'] = 23 * cation_proportions['Si_cation'] / anion_proportions['total_anion'] * 2
        cat_based13['Ti'] = 23 * cation_proportions['Ti_cation'] / anion_proportions['total_anion'] * 2
        cat_based13['Al'] = 23 * cation_proportions['Al_cation'] / anion_proportions['total_anion'] * 1.5
        cat_based13['Fe2+'] = 23 * df_s13['Fe2+_S13'] * (cation_proportions['total_cationFM'] / 13) / \
                              anion_proportions['total_anion']
        cat_based13['Fe3+'] = 23 * df_s13['Fe3+_S13'] * (cation_proportions['total_cationFM'] / 13) / \
                              anion_proportions['total_anion'] * 1.5
        cat_based13['Mn'] = 23 * cation_proportions['Mn_cation'] / anion_proportions['total_anion']
        cat_based13['Mg'] = 23 * cation_proportions['Mg_cation'] / anion_proportions['total_anion']
        cat_based13['Ca'] = 23 * cation_proportions['Ca_cation'] / anion_proportions['total_anion']
        cat_based13['Na'] = 23 * cation_proportions['Na_cation'] / anion_proportions['total_anion'] * 0.5
        cat_based13['K'] = 23 * cation_proportions['K_cation'] / anion_proportions['total_anion'] * 0.5
        cat_based13['Total'] = cat_based13.sum(axis='columns')
        cat_based13['adjustment'] = 23 / cat_based13['Total']

        # S15 Cations + Recheck
        cat_based15['Si'] = cat_based13['Si']
        cat_based15['Ti'] = cat_based13['Ti']
        cat_based15['Al'] = cat_based13['Al']
        cat_based15['Fe2+'] = 23 * df_s15['Fe2+_S15'] * (cation_proportions['total_cationFM+Ca'] / 15) / \
                              anion_proportions['total_anion']
        cat_based15['Fe3+'] = 23 * df_s15['Fe3+_S15'] * (cation_proportions['total_cationFM+Ca'] / 15) / \
                              anion_proportions['total_anion'] * 1.5
        cat_based15['Mn'] = cat_based13['Mn']
        cat_based15['Mg'] = cat_based13['Mg']
        cat_based15['Ca'] = cat_based13['Ca']
        cat_based15['Na'] = cat_based13['Na']
        cat_based15['K'] = cat_based13['K']
        cat_based15['Total'] = cat_based15.sum(axis='columns')
        cat_based15['adjustment'] = 23 / cat_based15['Total']

        # Final Formulas

        formulae_s13 = pd.DataFrame()
        formulae_s15 = pd.DataFrame()

        # S13 Cations
        formulae_s13['Si'] = (cat_based13['Si'] * cat_based13['adjustment']) / 2
        formulae_s13['Ti'] = (cat_based13['Ti'] * cat_based13['adjustment']) / 2
        formulae_s13['Al'] = (cat_based13['Al'] * cat_based13['adjustment']) / 1.5
        formulae_s13['Fe2+'] = (cat_based13['Fe2+'] * cat_based13['adjustment'])
        formulae_s13['Fe3+'] = (cat_based13['Fe3+'] * cat_based13['adjustment']) / 1.5
        formulae_s13['Mn'] = (cat_based13['Mn'] * cat_based13['adjustment'])
        formulae_s13['Mg'] = (cat_based13['Mg'] * cat_based13['adjustment'])
        formulae_s13['Ca'] = (cat_based13['Ca'] * cat_based13['adjustment'])
        formulae_s13['Na'] = (cat_based13['Na'] * cat_based13['adjustment']) / 0.5
        formulae_s13['K'] = (cat_based13['K'] * cat_based13['adjustment']) / 0.5

        # S15 Cations
        formulae_s15['Si'] = (cat_based15['Si'] * cat_based15['adjustment']) / 2
        formulae_s15['Ti'] = (cat_based15['Ti'] * cat_based15['adjustment']) / 2
        formulae_s15['Al'] = (cat_based15['Al'] * cat_based15['adjustment']) / 1.5
        formulae_s15['Fe2+'] = (cat_based15['Fe2+'] * cat_based15['adjustment'])
        formulae_s15['Fe3+'] = (cat_based15['Fe3+'] * cat_based15['adjustment']) / 1.5
        formulae_s15['Mn'] = (cat_based15['Mn'] * cat_based15['adjustment'])
        formulae_s15['Mg'] = (cat_based15['Mg'] * cat_based15['adjustment'])
        formulae_s15['Ca'] = (cat_based15['Ca'] * cat_based15['adjustment'])
        formulae_s15['Na'] = (cat_based15['Na'] * cat_based15['adjustment']) / 0.5
        formulae_s15['K'] = (cat_based15['K'] * cat_based15['adjustment']) / 0.5

        # return formulae_s13, formulae_s15
        # formulae_s13, formulae_s15 = final_formulaes(recheck13, recheck15)
    
        final_formulaes = formulae_s13.copy()
    
        for index, value in enumerate(formulae_s13):
            if (formulae_s13['Ca'].loc[index] + formulae_s13['Na'].loc[index]) < 1.34:
                final_formulaes.loc[index] = formulae_s15.loc[index]
            else:
                final_formulaes.loc[index] = formulae_s13.loc[index]
    
        # final_formulaes['Al_iv'] = 0
        for index, value in enumerate(final_formulaes['Si']):
            # Al iv
            if (final_formulaes['Si'].loc[index] + final_formulaes['Al'].loc[index]) >= 8:
                final_formulaes.loc[index, 'Al_iv'] = 8 - final_formulaes['Si'].loc[index]
                if final_formulaes['Al_iv'].loc[index] < 0:
                    final_formulaes.loc[index, 'Al_iv'] = 0
    
            else:
                final_formulaes.loc[index, 'Al_iv'] = final_formulaes['Al'].loc[index]
            # impedir Al negativo
            # Al vi
            final_formulaes['Al_vi'] = final_formulaes['Al'] - final_formulaes['Al_iv']
            # Si
            if final_formulaes['Si'].loc[index] > 8:
                final_formulaes['Si'].loc[index] = 8
            else:
                final_formulaes['Si'].loc[index] = final_formulaes['Si'].loc[index]
    
        final_formulaes['Soma_S2'] = final_formulaes.sum(axis='columns') - final_formulaes['Ca'] - \
                                     final_formulaes['Na'] - final_formulaes['K'] - final_formulaes['Al_vi'] - \
                                     final_formulaes['Al_iv']
    
        # Total in site (T):
        temp_df = pd.DataFrame()
        temp_df['Total_(T)'] = final_formulaes['Si'] + \
                               final_formulaes['Al'] + \
                               final_formulaes['Fe3+'] + \
                               final_formulaes['Ti']
        # Excess in site (T)
        temp_df['Excess_(T)'] = 0
        for index, value in enumerate(temp_df['Total_(T)']):
            if value > 8:
                temp_df.loc[index, 'Excess_(T)'] = value - 8
                if temp_df['Excess_(T)'].loc[index] > 5:
                    temp_df.loc[index, 'Excess_(T)'] = 5
            else:
                temp_df['Excess_(T)'].loc[index] = 0
    
        # Total in site C
        temp_df['Total_(C)'] = final_formulaes['Fe2+'] + \
                               final_formulaes['Mg'] + \
                               final_formulaes['Mn'] + \
                               temp_df['Excess_(T)']
    
        temp_df['Excess_(C)'] = 0
        for index, value in enumerate(temp_df['Total_(C)']):
            if value > 5:
                temp_df.loc[index, 'Excess_(C)'] = value - 5
                if temp_df.loc[index, 'Excess_(C)'] > 2:
                    temp_df.loc[index, 'Excess_(C)'] = 2
            else:
                temp_df.loc[index, 'Excess_(C)'] = 0
    
        # Total in site B
        temp_df['Total_(B)'] = final_formulaes['Ca'] + \
                               formulae_s13['Na'] + \
                               temp_df['Excess_(C)']
        # Excess in B
        temp_df['Excess_(B)'] = 0
        for index, value in enumerate(temp_df['Total_(B)']):
            if value > 2:
                temp_df.loc[index, 'Excess_(B)'] = value - 2
    
        # Excess in B restricted to Na
        temp_df['Excess_(B)_Na'] = 0
        for index, value in enumerate(temp_df['Excess_(B)']):
            if value > final_formulaes['Na'].loc[index]:
                temp_df.loc[index, 'Excess_(B)_Na'] = final_formulaes['Na'].loc[index]
            else:
                temp_df.loc[index, 'Excess_(B)_Na'] = temp_df['Excess_(B)'].loc[index]
    
        # Total in site A
        temp_df['Total_(A)'] = temp_df['Excess_(B)_Na'] + final_formulaes['K']
    
        # Excess in site A
        temp_df['Excess_(A)'] = 0
        for index, value in enumerate(temp_df['Total_(A)']):
            if value > 1:
                temp_df.loc[index, 'Excess_(A)'] = value - 1
            else:
                temp_df.loc[index, 'Excess_(A)'] = 0
    
        # Excess in site A if K > 1:
        temp_df['Excess_(A)K'] = 0
        for index, value in enumerate(final_formulaes['K']):
            if value > 1:
                temp_df.loc[index, 'Excess_(A)K'] = temp_df['Excess_(B)_Na'].loc[index]
            else:
                temp_df.loc[index, 'Excess_(A)K'] = temp_df['Excess_(A)'].loc[index]
    
        # Verifying if (Na + K) > 1:
        temp_df['(Na + K)'] = temp_df['Total_(A)']
        for index, value in enumerate(temp_df['(Na + K)']):
            if value > 1:
                temp_df.loc[index, '(Na + K)'] = 1
            else:
                temp_df.loc[index, '(Na + K)'] = value
            if final_formulaes['K'].loc[index] > 1:
                final_formulaes.loc[index, '(Na + K)'] = final_formulaes['K'].loc[index]
    
        # Na in site (B):
        temp_df['Na_(B)'] = final_formulaes['Na'] - temp_df['Excess_(B)_Na'] + temp_df['Excess_(A)K']
    
        # Ca + Na in site (B):
        temp_df['(Ca+Na)(B)'] = temp_df['Na_(B)'] + final_formulaes['Ca']
    
        # Mg/(Mg + Fe2+) ratio
        temp_df['Mg/(Mg + Fe2+)'] = final_formulaes['Mg'] / (final_formulaes['Mg'] + final_formulaes['Fe2+'])
    
        # Fe3+/(Fe3+ + Al vi) ratio
        temp_df['Fe3+/(Fe3+ + Al_vi)'] = final_formulaes['Fe3+'] / \
                                         (final_formulaes['Fe3+'] + final_formulaes['Al_vi'])
    
        final_formulaes['(Na+K)(A)'] = temp_df['(Na + K)']
        final_formulaes['Na_(B)'] = temp_df['Na_(B)']
        final_formulaes['(Ca+Na)(B)'] = temp_df['(Ca+Na)(B)']
        final_formulaes['Mg/(Mg + Fe2+)'] = temp_df['Mg/(Mg + Fe2+)']
        final_formulaes['Fe3+/(Fe3+ + Al_vi)'] = temp_df['Fe3+/(Fe3+ + Al_vi)']

        self.epma_normalized = final_formulaes
        return final_formulaes

    def leake1997(self):
        """
        :return: Leake (1997) Diagram for (Na+K)(A) > 0.5 apfu
        """
        plt.subplots(figsize=(8, 6))
        plt.xlim(7.5, 4.5)
        plt.ylim(0, 1)
        plt.plot([4.5, 7.5, 0, 1], [0.5, 0.5, 0.5, 0.5], color='black', zorder=0, linewidth=0.5)
        plt.plot([6.5, 6.5, 6.5, 6.5], [6.5, 6.5, 0, 1], color='black', zorder=0, linewidth=0.5)
        plt.plot([5.5, 5.5, 5.5, 5.5], [5.5, 5.5, 0, 1], color='black', zorder=0, linewidth=0.5)
        sns.scatterplot(x=self.epma_normalized["Si"], y=self.epma_normalized["Mg/(Mg + Fe2+)"],
                        hue=self.epma_dataframe["Sample"], s=100)
        plt.ylabel("$\mathrm{Mg/}\mathrm{(Mg}+\mathrm{Fe}^{2})$", fontsize=16)
        plt.xlabel('Si', fontsize=16)
        plt.text(7.15, 0.75, 'Edenite')
        plt.text(7.20, 0.25, 'Fe-edenite')
        plt.text(6.2, 0.75, 'Pargasite')
        plt.text(6.35, 0.25, 'Mg-hastingsite')
        plt.text(5.35, 0.75, 'Mg-sadanagaite')
        plt.text(5.35, 0.25, 'Fe-Tschermakite')
        plt.show()


########################################################################################################################
#                                           Normalização de elementos testes                                           #
########################################################################################################################


df = pd.read_csv("liliana_amp.csv")
pd.set_option('display.max_columns', None)
df['Fe2O3'] = 0

print(df.head())

meu_anfibolios = Amphibole(df)
# print(meu_anfibolios.epma_dataframe)
meu_anfibolios.normalization()
meu_anfibolios.leake1997()

