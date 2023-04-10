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
        df_anion['Cr_anion'] = (epma_table['Cr2O3'] * 3) / 151.9902
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
        df_cation['Cr_cation'] = (epma_table['Cr2O3'] * 2) / 151.9902
        df_cation['Fe3+_cation'] = (epma_table['Fe2O3'] * 2) / 60.0843
        df_cation['Fe2+_cation'] = epma_table['FeO'] / 71.8464
        df_cation['Mn_cation'] = epma_table['MnO'] / 70.9374
        df_cation['Mg_cation'] = epma_table['MgO'] / 40.3044
        df_cation['Ca_cation'] = epma_table['CaO'] / 56.0794
        df_cation['Na_cation'] = (epma_table['Na2O'] * 2) / 61.9774
        df_cation['K_cation'] = (epma_table['K2O'] * 2) / 94.196
        df_cation['total_cation'] = df_cation.sum(axis='columns') - df_cation['Ca_cation'] - df_cation['Sr_cation']
        df_cation['total_cation + Ca + Sr'] = df_cation.sum(axis='columns')  # falta adicionar o Sr

        return df_cation

    return normalization_anion(), normalization_cation()


def normalization_pt3a(dataframe):
    """S13 normalized based on anions"""
    df_s13 = pd.DataFrame()

    df_s13['Si_S13'] = 13 * dataframe['Si_anion'] / dataframe['total_cation']
    df_s13['Ti_S13'] = 13 * dataframe['Ti_anion'] / dataframe['total_cation']
    df_s13['Al_S13'] = 13 * dataframe['Al_anion'] / dataframe['total_cation']
    df_s13['Fe3+_S13'] = 13 * dataframe['Fe3_anion'] / dataframe['total_cation']
    df_s13['Fe2+_S13'] = 13 * dataframe['Fe2_anion'] / dataframe['total_cation']
    df_s13['Mn_S13'] = 13 * dataframe['Mn_anion'] / dataframe['total_cation']
    df_s13['Mg_S13'] = 13 * dataframe['Mg_anion'] / dataframe['total_cation']
    df_s13['Ca_S13'] = 13 * dataframe['Ca_anion'] / dataframe['total_cation']
    df_s13['Na_S13'] = 13 * dataframe['Na_anion'] / dataframe['total_cation']
    df_s13['K_S13'] = 13 * dataframe['K_anion'] / dataframe['total_cation']
    df_s13['total_S13'] = df_s13.sum(axis='columns')

    return df_s13


def normalization_pt3b(dataframe):
    """S15 normalized based on anions"""

    df_s15 = pd.DataFrame()

    df_s15['Si_S15'] = 15 * dataframe['Si_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Ti_S15'] = 15 * dataframe['Ti_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Al_S15'] = 15 * dataframe['Al_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Cr_S15'] = 15 * dataframe['Cr_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Fe3+_S15'] = 15 * dataframe['Fe3_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Fe2+_S15'] = 15 * dataframe['Fe2_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Mn_S15'] = 15 * dataframe['Mn_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Mg_S15'] = 15 * dataframe['Mg_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Ca_S15'] = 15 * dataframe['Ca_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Na_S15'] = 15 * dataframe['Na_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['K_S3'] = 15 * dataframe['K_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['total_S15'] = df_s15.sum(axis='columns')

    return df_s15


def iron_ox_state(dataframe):
    """
    Iron oxidation state calculation
    """

    # Cálculo da oxidação para S13
    fe3 = 2 * (23-dataframe['total_S13'])
    if fe3 < 0:
        fe3 = 0
    else:
        fe3 = fe3

    if fe3 > dataframe['Fe2+_S13']:
        fe3 = dataframe['Fe2+_S13']
    else:
        fe3 = fe3

    fe2 = dataframe['Fe2+_S13'] - fe3

    # Cálculo da oxidação para S15
    fe3 = 2 * (23 - dataframe['total_S15'])
    if fe3 < 0:
        fe3 = 0
    else:
        fe3 = fe3

    if fe3 > dataframe['Fe2+_S15']:
        fe3 = dataframe['Fe2+_S15']
    else:
        fe3 = fe3

    fe2 = dataframe['Fe2+_S15'] - fe3


def cation_o23_s13(dataframe):
    cat_df = pd.DataFrame()

    cat_df['Si'] = 23 * dataframe['Si_cation'] / dataframe['total_anion']
    cat_df['Ti'] = 23 * dataframe['Ti_cation'] / dataframe['total_anion']
    cat_df['Al'] = 23 * dataframe['Al_cation'] / dataframe['total_anion']
    cat_df['Fe3+'] = 23 * dataframe['Fe3+_S13'] * (dataframe['total_cation'] / 13) / dataframe['total_anion']
    cat_df['Fe2+'] = 23 * dataframe['Fe2+_S13'] * (dataframe['total_cation'] / 13) / dataframe['total_anion']
    cat_df['Mn'] = 23 * dataframe['Si_cation'] / dataframe['total_anion']
    cat_df['Mg'] = 23 * dataframe['Ti_cation'] / dataframe['total_anion']
    cat_df['Ca'] = 23 * dataframe['Al_cation'] / dataframe['total_anion']
    cat_df['Na'] = 23 * dataframe['Si_cation'] / dataframe['total_anion']
    cat_df['K'] = 23 * dataframe['Ti_cation'] / dataframe['total_anion']
    cat_df['Total'] = cat_df.sum(axis='columns')

    return cat_df
