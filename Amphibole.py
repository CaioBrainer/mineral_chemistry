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

    def normalization_anion():
        """Anions proportion calculation
            :returns: Anion normalized proportions epma_table
        """
        df_anion = pd.DataFrame()

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

        return df_anion

    def normalization_cation():
        """Cation proportion calculation
            :returns: Cation normalized proportions epma_table
        """
        df_cation = pd.DataFrame()

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
        df_cation['total_cationFM'] = df_cation.sum(axis='columns') - df_cation['Ca_cation'] - df_cation['Na_cation'] - df_cation['K_cation'] - df_cation['total_cationFM+Ca']

        return df_cation

    an_prop, cat_prop = normalization_anion(), normalization_cation()

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
                s13_normalized['Fe3+_S13'][index] = 0
            else:
                s13_normalized['Fe3+_S13'][index] = value
        # S15
        s15_normalized['Fe3+_S15'] = 2 * (23 - s15_normalized['total_S15'])

        for index, value in enumerate(s15_normalized['Fe3+_S15']):
            if value < 0:
                s15_normalized['Fe3+_S15'][index] = 0
            else:
                s15_normalized['Fe3+_S15'][index] = value

        # Fe3+ and Fe2+ calculation S13
        for index, value in enumerate(s13_normalized['Fe3+_S13']):
            if value > s13_normalized['Fe2+_S13'][index]:
                s13_normalized['Fe3+_S13'][index] = s13_normalized['Fe2+_S13'][index]
            else:
                s13_normalized['Fe3+_S13'][index] = s13_normalized['Fe3+_S13'][index]

        s13_normalized['Fe2+_S13'] = s13_normalized['Fe2+_S13'] - s13_normalized['Fe3+_S13']

        # Fe3+ and Fe2+ calculation S15
        for index, value in enumerate(s15_normalized['Fe3+_S15']):
            if value > s15_normalized['Fe2+_S15'][index]:
                s15_normalized['Fe3+_S15'][index] = s15_normalized['Fe2+_S15'][index]
            else:
                s15_normalized['Fe3+_S15'][index] = s15_normalized['Fe3+_S15'][index]

        s15_normalized['Fe2+_S15'] = s15_normalized['Fe2+_S15'] - s15_normalized['Fe3+_S15']

        return s13, s15

    s13_ox, s15_ox = oxid(s13, s15)

    # Cations based on 23 Oxygens S13 Cations
    def cat_based_on_23o_s13(an, cat, oxid_state):
        df_cat = pd.DataFrame()
        df_cat['Si'] = 23 * cat['Si_cation'] / an['total_anion']
        df_cat['Ti'] = 23 * cat['Ti_cation'] / an['total_anion']
        df_cat['Al'] = 23 * cat['Al_cation'] / an['total_anion']
        df_cat['Fe2+'] = 23 * oxid_state['Fe3+_S13'] * (cat['total_cationFM'] / 13) / an['total_anion']
        df_cat['Fe3+'] = 23 * oxid_state['Fe2+_S13'] * (cat['total_cationFM'] / 13) / an['total_anion']

        return df_cat

    final = cat_based_on_23o_s13(an_prop, cat_prop, s13_ox)

    return final

