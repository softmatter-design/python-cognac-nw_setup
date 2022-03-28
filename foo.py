import nw_setup

def main():
    nw_setup.hello()

    arr = nw_setup.get_array()
    print(arr)

    df = nw_setup.get_df()
    print(df)

    nw_setup.get_date()

if __name__ == '__main__':
    main()