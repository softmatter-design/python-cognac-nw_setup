import numpy as np
import pandas as pd
import datetime

def hello():
    print("hello_world")


def get_array():
    return np.array([1, 2, 3])


def get_df():
    return pd.Series([1, 2, 3])


def get_date():
    print(datetime.datetime(2019, 12, 1, 1, 1, 1))