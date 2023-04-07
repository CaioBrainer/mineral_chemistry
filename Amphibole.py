# Projeto de programa para cálculo e projeção de química mineral
# V 0.0.1
import pandas as pd

'''###############################################################
                   Normalização de elementos
###############################################################'''


# Anfibólio

def normalization_pt1(dataframe):
    """Anions proportion calculation"""

    df_anion = pd.DataFrame()

    df_anion['Si_anion'] = (dataframe['SiO2'] * 2) / 60.0843
    df_anion['Ti_anion'] = (dataframe['TiO2'] * 2) / 79.8788
    df_anion['Al_anion'] = (dataframe['Al2O3'] * 3) / 101.9602
    df_anion['Cr_anion'] = (dataframe['Cr2O3'] * 3) / 151.9902
    df_anion['Fe3_anion'] = (dataframe['Fe2O3'] * 3) / 60.0843
    df_anion['Fe2_anion'] = dataframe['FeO'] / 71.8464
    df_anion['Mn_anion'] = dataframe['MnO'] / 70.9374
    df_anion['Mg_anion'] = dataframe['MgO'] / 40.3044
    df_anion['Ca_anion'] = dataframe['CaO'] / 56.0794
    df_anion['Na_anion'] = dataframe['Na2O'] / 61.9774
    df_anion['K_anion'] = dataframe['K2O'] / 94.196
    df_anion['total_anion'] = df_anion.sum(axis='columns')

    return df_anion


def normalization_pt2(dataframe):
    """Cation proportion calculation"""
    df_cation = pd.DataFrame()

    df_cation['Si_cation'] = dataframe['SiO2'] / 60.0843
    df_cation['Ti_cation'] = dataframe['TiO2'] / 79.8788
    df_cation['Al_cation'] = (dataframe['Al2O3'] * 2) / 101.9602
    df_cation['Cr_cation'] = (dataframe['Cr2O3'] * 2) / 151.9902
    df_cation['Fe3_cation'] = (dataframe['Fe2O3'] * 2) / 60.0843
    df_cation['Fe2_cation'] = dataframe['FeO'] / 71.8464
    df_cation['Mn_cation'] = dataframe['MnO'] / 70.9374
    df_cation['Mg_cation'] = dataframe['MgO'] / 40.3044
    df_cation['Ca_cation'] = dataframe['CaO'] / 56.0794
    df_cation['Na_cation'] = (dataframe['Na2O'] * 2) / 61.9774
    df_cation['K_cation'] = (dataframe['K2O'] * 2) / 94.196
    df_cation['total_cation'] = df_cation.sum(axis='columns') - df_cation['Ca_cation'] - df_cation['Sr_cation']
    df_cation['total_cation + Ca + Sr'] = df_cation.sum(axis='columns')  # falta adicionar o Sr

    return df_cation


def normalization_pt3a(dataframe):
    """S13 normalized based on anions"""
    df_s13 = pd.DataFrame()

    df_s13['Si_S13'] = 13 * dataframe['Si_anion'] / dataframe['total_cation']
    df_s13['Ti_S13'] = 13 * dataframe['Ti_anion'] / dataframe['total_cation']
    df_s13['Al_S13'] = 13 * dataframe['Al_anion'] / dataframe['total_cation']
    df_s13['Cr_S13'] = 13 * dataframe['Cr_anion'] / dataframe['total_cation']
    df_s13['Fe3_S13'] = 13 * dataframe['Fe3_anion'] / dataframe['total_cation']
    df_s13['Fe2_S13'] = 13 * dataframe['Fe2_anion'] / dataframe['total_cation']
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
    df_s15['Fe3_S15'] = 15 * dataframe['Fe3_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Fe2_S15'] = 15 * dataframe['Fe2_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Mn_S15'] = 15 * dataframe['Mn_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Mg_S15'] = 15 * dataframe['Mg_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Ca_S15'] = 15 * dataframe['Ca_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['Na_S15'] = 15 * dataframe['Na_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['K_S3'] = 15 * dataframe['K_anion'] / dataframe['total_cation + Ca + Sr']
    df_s15['total_S15'] = df_s15.sum(axis='columns')

    return df_s15


'''def iron_state():
    s13_fe3 = 2*F-total anions
    se s13_fe3 < 0 = 0, se não s13_fe3 = s13_fe3 = 2*F-s13total anions'''
