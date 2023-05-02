# Projeto de programa para cálculo e projeção de química mineral
# V 0.0.1
import pandas as pd

'''###############################################################
                   Normalização de elementos
###############################################################'''


# Anfibólio

def normalization(epma_table):
    """Anions proportion calculation

                :param epma_table: Amphibole EPMA tabular data.
                :returns: Anion and Cation normalized proportions dataframes
    """

    def norm():
        """Anions and Cations proportion calculation
            :returns: Anion and Cation normalized proportions epma_table
        """

        df_anion = pd.DataFrame()
        df_cation = pd.DataFrame()

        # Anion normalized calculation
        df_anion['Si_anion'] = (epma_table['SiO2'] * 2) / 60.0843
        df_anion['Ti_anion'] = (epma_table['TiO2'] * 2) / 79.8788
        df_anion['Al_anion'] = (epma_table['Al2O3'] * 3) / 101.9602
        df_anion['Fe3_anion'] = (epma_table['Fe2O3'] * 3) / 60.0843
        df_anion['Fe2_anion'] = epma_table['FeO'] / 71.8464
        df_anion['Mn_anion'] = epma_table['MnO'] / 70.9374
        df_anion['Mg_anion'] = epma_table['MgO'] / 40.3044
        df_anion['Ca_anion'] = epma_table['CaO'] / 56.0794
        df_anion['Na_anion'] = epma_table['Na2O'] / 61.9774
        df_anion['K_anion'] = epma_table['K2O'] / 94.196
        df_anion['total_anion'] = df_anion.sum(axis='columns')

        # Anion normalized calculation
        df_cation['Si_cation'] = epma_table['SiO2'] / 60.0843
        df_cation['Ti_cation'] = epma_table['TiO2'] / 79.8788
        df_cation['Al_cation'] = (epma_table['Al2O3'] * 2) / 101.9602
        df_cation['Fe3+_cation'] = (epma_table['Fe2O3'] * 2) / 60.0843
        df_cation['Fe2+_cation'] = epma_table['FeO'] / 71.8464
        df_cation['Mn_cation'] = epma_table['MnO'] / 70.9374
        df_cation['Mg_cation'] = epma_table['MgO'] / 40.3044
        df_cation['Ca_cation'] = epma_table['CaO'] / 56.0794
        df_cation['Na_cation'] = (epma_table['Na2O'] * 2) / 61.9774
        df_cation['K_cation'] = (epma_table['K2O'] * 2) / 94.196
        df_cation['total_cationFM+Ca'] = df_cation.sum(axis='columns') - df_cation['Na_cation'] - df_cation['K_cation']
        df_cation['total_cationFM'] = df_cation.sum(axis='columns') - df_cation['Ca_cation'] - df_cation['Na_cation'] \
                                      - df_cation['K_cation'] - df_cation['total_cationFM+Ca']

        return df_anion, df_cation

    an_prop, cat_prop = norm()

    def normalization_13(anion_proportions, cat_proportions):
        """S13 normalized based on anions"""

        anion = anion_proportions
        cation = cat_proportions
        df_s13 = pd.DataFrame()

        df_s13['Si_S13'] = 13 * anion['Si_anion'] / cation['total_cationFM']
        df_s13['Ti_S13'] = 13 * anion['Ti_anion'] / cation['total_cationFM']
        df_s13['Al_S13'] = 13 * anion['Al_anion'] / cation['total_cationFM']
        df_s13['Fe2+_S13'] = 13 * anion['Fe2_anion'] / cation['total_cationFM']
        df_s13['Fe3+_S13'] = 0
        df_s13['Mn_S13'] = 13 * anion['Mn_anion'] / cation['total_cationFM']
        df_s13['Mg_S13'] = 13 * anion['Mg_anion'] / cation['total_cationFM']
        df_s13['Ca_S13'] = 13 * anion['Ca_anion'] / cation['total_cationFM']
        df_s13['Na_S13'] = 13 * anion['Na_anion'] / cation['total_cationFM']
        df_s13['K_S13'] = 13 * anion['K_anion'] / cation['total_cationFM']
        df_s13['total_S13'] = df_s13.sum(axis='columns')

        return df_s13

    def normalization_15(anion_proportions, cat_proportions):
        """S15 normalized based on anions"""

        anion = anion_proportions
        cation = cat_proportions
        df_s15 = pd.DataFrame()

        df_s15['Si_S15'] = 15 * anion['Si_anion'] / cation['total_cationFM+Ca']
        df_s15['Ti_S15'] = 15 * anion['Ti_anion'] / cation['total_cationFM+Ca']
        df_s15['Al_S15'] = 15 * anion['Al_anion'] / cation['total_cationFM+Ca']
        df_s15['Fe2+_S15'] = 15 * anion['Fe2_anion'] / cation['total_cationFM+Ca']
        df_s15['Fe3+_S15'] = 0
        df_s15['Mn_S15'] = 15 * anion['Mn_anion'] / cation['total_cationFM+Ca']
        df_s15['Mg_S15'] = 15 * anion['Mg_anion'] / cation['total_cationFM+Ca']
        df_s15['Ca_S15'] = 15 * anion['Ca_anion'] / cation['total_cationFM+Ca']
        df_s15['Na_S15'] = 15 * anion['Na_anion'] / cation['total_cationFM+Ca']
        df_s15['K_S15'] = 15 * anion['K_anion'] / cation['total_cationFM+Ca']
        df_s15['total_S15'] = df_s15.sum(axis='columns')

        return df_s15

    s13, s15 = normalization_13(an_prop, cat_prop), normalization_15(an_prop, cat_prop)

    def oxid(s13_normalized, s15_normalized):
        # Fe3+ calculation

        # S13
        s13_normalized['Fe3+_S13'] = 2 * (23 - s13_normalized['total_S13'])

        for index, value in enumerate(s13_normalized['Fe3+_S13']):
            if value < 0:
                s13_normalized['Fe3+_S13'].loc[index] = 0
            else:
                s13_normalized['Fe3+_S13'].loc[index] = value
        # S15
        s15_normalized['Fe3+_S15'] = 2 * (23 - s15_normalized['total_S15'])

        for index, value in enumerate(s15_normalized['Fe3+_S15']):
            if value < 0:
                s15_normalized['Fe3+_S15'].loc[index] = 0
            else:
                s15_normalized['Fe3+_S15'].loc[index] = value

        # Fe3+ and Fe2+ calculation S13
        for index, value in enumerate(s13_normalized['Fe3+_S13']):
            if value > s13_normalized['Fe2+_S13'].loc[index]:
                s13_normalized.loc[index, 'Fe3+_S13'] = s13_normalized['Fe2+_S13'].loc[index]
            else:
                s13_normalized.loc[index, 'Fe3+_S13'] = s13_normalized['Fe3+_S13'].loc[index]

        s13_normalized['Fe2+_S13'] = s13_normalized['Fe2+_S13'] - s13_normalized['Fe3+_S13']

        # Fe3+ and Fe2+ calculation S15
        for index, value in enumerate(s15_normalized['Fe3+_S15']):
            if value > s15_normalized['Fe2+_S15'].loc[index]:
                s15_normalized['Fe3+_S15'].loc[index] = s15_normalized['Fe2+_S15'].loc[index]
            else:
                s15_normalized['Fe3+_S15'].loc[index] = s15_normalized['Fe3+_S15'].loc[index]

        s15_normalized['Fe2+_S15'] = s15_normalized['Fe2+_S15'] - s15_normalized['Fe3+_S15']

        return s13, s15

    s13_ox, s15_ox = oxid(s13, s15)

    # Cations based on 23 Oxygens
    def cat_based_on_23o(an, cat, oxid_state13, oxid_state15):
        df_13cat = pd.DataFrame()
        df_15cat = pd.DataFrame()

        # S13 Cations
        df_13cat['Si'] = 23 * cat['Si_cation'] / an['total_anion']
        df_13cat['Ti'] = 23 * cat['Ti_cation'] / an['total_anion']
        df_13cat['Al'] = 23 * cat['Al_cation'] / an['total_anion']
        df_13cat['Fe2+'] = 23 * oxid_state13['Fe2+_S13'] * (cat['total_cationFM'] / 13) / an['total_anion']
        df_13cat['Fe3+'] = 23 * oxid_state13['Fe3+_S13'] * (cat['total_cationFM'] / 13) / an['total_anion']
        df_13cat['Mn'] = 23 * cat['Mn_cation'] / an['total_anion']
        df_13cat['Mg'] = 23 * cat['Mg_cation'] / an['total_anion']
        df_13cat['Ca'] = 23 * cat['Ca_cation'] / an['total_anion']
        df_13cat['Na'] = 23 * cat['Na_cation'] / an['total_anion']
        df_13cat['K'] = 23 * cat['K_cation'] / an['total_anion']

        # S15 Cations
        df_15cat['Si'] = df_13cat['Si']
        df_15cat['Ti'] = df_13cat['Ti']
        df_15cat['Al'] = df_13cat['Al']
        df_15cat['Fe2+'] = 23 * oxid_state15['Fe2+_S15'] * (cat['total_cationFM+Ca'] / 15) / an['total_anion']
        df_15cat['Fe3+'] = 23 * oxid_state15['Fe3+_S15'] * (cat['total_cationFM+Ca'] / 15) / an['total_anion']
        df_15cat['Mn'] = df_13cat['Mn']
        df_15cat['Mg'] = df_13cat['Mg']
        df_15cat['Ca'] = df_13cat['Ca']
        df_15cat['Na'] = df_13cat['Na']
        df_15cat['K'] = df_13cat['K']

        return df_13cat, df_15cat

    cat_based13, cat_based15 = cat_based_on_23o(an_prop, cat_prop, s13_ox, s15_ox)

    def recheck_anion_sum(cb13, cb15):

        # S13 Cations
        cb13['Si'] = cb13['Si'] * 2
        cb13['Ti'] = cb13['Ti'] * 2
        cb13['Al'] = cb13['Al'] * 1.5
        cb13['Fe2+'] = cb13['Fe2+']
        cb13['Fe3+'] = cb13['Fe3+'] * 1.5
        cb13['Mn'] = cb13['Mn']
        cb13['Mg'] = cb13['Mg']
        cb13['Ca'] = cb13['Ca']
        cb13['Na'] = cb13['Na'] * 0.5
        cb13['K'] = cb13['K'] * 0.5
        cb13['Total'] = cb13.sum(axis='columns')
        cb13['adjustment'] = 23 / cb13['Total']

        # S15 Cations
        cb15['Si'] = cb15['Si'] * 2
        cb15['Ti'] = cb15['Ti'] * 2
        cb15['Al'] = cb15['Al'] * 1.5
        cb15['Fe2+'] = cb15['Fe2+']
        cb15['Fe3+'] = cb15['Fe3+'] * 1.5
        cb15['Mn'] = cb15['Mn']
        cb15['Mg'] = cb15['Mg']
        cb15['Ca'] = cb15['Ca']
        cb15['Na'] = cb15['Na'] * 0.5
        cb15['K'] = cb15['K'] * 0.5
        cb15['Total'] = cb15.sum(axis='columns')
        cb15['adjustment'] = 23 / cb15['Total']

        return cb13, cb15

    recheck13, recheck15 = recheck_anion_sum(cat_based13, cat_based15)

    def final_formulaes(_s13, _s15):
        """
        :param _s15: tabular data of Recheck sum anions 13 cations
        :param _s13: tabular data of Recheck sum anions 13 cations
        :return: formulae based on 23 Oxygens for 13 and 15 cations
        """

        formulae_s13 = pd.DataFrame()
        formulae_s15 = pd.DataFrame()

        # S13 Cations
        formulae_s13['Si'] = (_s13['Si'] * _s13['adjustment']) / 2
        formulae_s13['Ti'] = (_s13['Ti'] * _s13['adjustment']) / 2
        formulae_s13['Al'] = (_s13['Al'] * _s13['adjustment']) / 1.5
        formulae_s13['Fe2+'] = (_s13['Fe2+'] * _s13['adjustment'])
        formulae_s13['Fe3+'] = (_s13['Fe3+'] * _s13['adjustment']) / 1.5
        formulae_s13['Mn'] = (_s13['Mn'] * _s13['adjustment'])
        formulae_s13['Mg'] = (_s13['Mg'] * _s13['adjustment'])
        formulae_s13['Ca'] = (_s13['Ca'] * _s13['adjustment'])
        formulae_s13['Na'] = (_s13['Na'] * _s13['adjustment']) / 0.5
        formulae_s13['K'] = (_s13['K'] * _s13['adjustment']) / 0.5

        # S15 Cations
        formulae_s15['Si'] = (_s15['Si'] * _s15['adjustment']) / 2
        formulae_s15['Ti'] = (_s15['Ti'] * _s15['adjustment']) / 2
        formulae_s15['Al'] = (_s15['Al'] * _s15['adjustment']) / 1.5
        formulae_s15['Fe2+'] = (_s15['Fe2+'] * _s15['adjustment'])
        formulae_s15['Fe3+'] = (_s15['Fe3+'] * _s15['adjustment']) / 1.5
        formulae_s15['Mn'] = (_s15['Mn'] * _s15['adjustment'])
        formulae_s15['Mg'] = (_s15['Mg'] * _s15['adjustment'])
        formulae_s15['Ca'] = (_s15['Ca'] * _s15['adjustment'])
        formulae_s15['Na'] = (_s15['Na'] * _s15['adjustment']) / 0.5
        formulae_s15['K'] = (_s15['K'] * _s15['adjustment']) / 0.5

        return formulae_s13, formulae_s15

    formulae_s13, formulae_s15 = final_formulaes(recheck13, recheck15)

    def selection(formulaes13, formulaes15):
        """

        :param formulaes13: tabular data based on 23 Oxygens for 13 cations
        :param formulaes15: tabular data based on 23 Oxygens for 15 cations
        :return: select S13 or S15 cations data and returns a normalized dataset
        """

        formulaes13 = formulaes13
        formulaes15 = formulaes15
        final_formulaes = formulaes13.copy()

        for index, value in enumerate(formulaes13):
            if (formulaes13['Ca'].loc[index] + formulaes13['Na'].loc[index]) < 1.34:
                final_formulaes.loc[index] = formulaes15.loc[index]
            else:
                final_formulaes.loc[index] = formulaes13.loc[index]

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

        final_formulaes['Total_(T)'] = temp_df['Total_(T)']
        final_formulaes['Excess_(T)'] = temp_df['Excess_(T)']
        final_formulaes['Total_(C)'] = temp_df['Total_(C)']
        final_formulaes['Excess_(C)'] = temp_df['Excess_(C)']
        final_formulaes['Total_(B)'] = temp_df['Total_(B)']
        final_formulaes['Excess_(B)'] = temp_df['Excess_(B)']
        final_formulaes['Excess_(B)_Na'] = temp_df['Excess_(B)_Na']
        final_formulaes['Total_(A)'] = temp_df['Total_(A)']

        return final_formulaes

    formula = selection(formulae_s13, formulae_s15)

    return formula

########################################################################################################################
#                                           Normalização de elementos testes                                           #
########################################################################################################################


import pandas as pd
df = pd.read_csv("example_table.csv")
pd.set_option('display.max_columns', None)
df['Fe2O3'] = 0
normalized_dataset = normalization(df)

print(normalized_dataset.head(20))

