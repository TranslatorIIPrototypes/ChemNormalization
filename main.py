from src.ChemNormalization import ChemNormalization


if __name__ == '__main__':
    # instantiate the class that does all the work
    cn = ChemNormalization()

    # call to load redis instances with normalized chemical substance data
    success: bool = cn.load()

    # check the return
    if not success:
        cn.print_debug_msg(f'Failed to load chemical normalization data.', True)
    else:
        cn.print_debug_msg(f'Success', True)
