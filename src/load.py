from ChemNormalization import ChemNormalization


if __name__ == '__main__':
    # instantiate the class that does all the work
    cn = ChemNormalization()

    # call to load redis instances with normalized chemical substance data
    success: bool = cn.load()

    # check the return
    if not success:
        print(f'Failed to load redis with chemical normalization data.')
    else:
        print(f'Success')
