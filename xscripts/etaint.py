import os
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy.io import fits
import xspec as xs


THREADS = "2"

os.environ["OPENBLAS_NUM_THREADS"] = THREADS
os.environ["VECLIB_MAXIMUM_THREADS"] = THREADS
os.environ["MKL_NUM_THREADS"] = THREADS
os.environ["NUMEXPR_NUM_THREADS"] = THREADS
os.environ["OMP_NUM_THREADS"] = THREADS

xs.Xset.allowPrompting = False



if __name__ == '__main__':
    '''pyetaint loads an MCMC Chain to obtain the quantiles of the intrinsic feedback'''

    # Organize the Parser to get Chain file, burn-in and samples to be used.
    parser = argparse.ArgumentParser(prog='pyetaint',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description='Obtain XSET variables based on XSPEC MCMC Chain FITS files.')

    parser.add_argument("chain", help="Path to XSPEC Chain FITS file", type=str, default='chain.fits')
    parser.add_argument("model", help="Name of the Model component to evaluate", type=str, default='vkompthdk')
    parser.add_argument("--pathtomodel", help="Full Path to Model Library", type=str, default='')
    parser.add_argument("--burn", help="Samples to Burn In", type=int, default=0, nargs='?')
    parser.add_argument("--samples",help="Number of random samples to include in the output table", type=int, default=-1, nargs='?')
    args = parser.parse_args()

    # Use the parsed arguments to get the selected data from the Chain FITS file
    chainName = args.chain
    model_variant = args.model
    PATHTOMODEL = args.pathtomodel
    BurnIn = int(args.burn)
    Samples = int(args.samples)

    try:
        chain = fits.open(chainName)
        nFields = int(chain[1].header['TFIELDS'])
        ChainLength = int(chain[1].header['NAXIS2'])

        if Samples<0:
            Samples = ChainLength
            idx = np.arange(min(BurnIn,0), ChainLength)
        else:
            idx = np.random.randint(low=min(BurnIn,0), high=ChainLength, size=Samples)

        print()
        print('===============================================')
        print(' Loading Chain: {}'.format(chainName))
        print(' Chain Length: {}'.format(ChainLength))
        print(' Number of Fields: {}'.format(nFields))
        print('===============================================')
        print(' Burning in {} samples'.format(BurnIn))
        print(' Using {} samples for output purposes'.format(Samples))
        print('===============================================')
        print()
        _ = input(f'\nChain FITS file successfully loaded. Press any key to continue... ')
        print()
    except:
        print(f'\n\n   ERROR: Could not load {chainName} chain FITS file. Exiting... \n')
        exit(1)

    # Create the DataFrame for the CornerPlot and fill it with the Data and Titles.
    df = pd.DataFrame()
    titles, ttypes, tnums, tnames = [], [], [], []

    # Read FITS chain file. Store tnums, tnames, titles and data.
    for i in range(nFields):
        ttype = chain[1].header['TTYPE{}'.format(i+1)]
        titles.append(ttype)
        df[ttype] = chain[1].data[ttype][idx]

   
    # Define the model variants
    model_variants = ['vkompthbb', 'vkompthdk', 'vkdualbb', 'vkdualdk']

    # Get parameters for current model variant
    if model_variant not in model_variants:
        print(f'\n\n   ERROR: Could not understand model name given: {model_variant}. Exiting... \n')
        exit(1)
    elif model_variant[:3] == 'vko':
        params = ['kTs', 'kTe', 'gam', 'size', 'eta', 'af', 'DHext', 'reflag', 'norm']
    else:
        params = ['kTs1', 'kTs2', 'kTe1', 'kTe2', 'gam1', 'gam2', 'size1', 'size2',
                  'eta1', 'eta2', 'af', 'DHext1', 'DHext2', 'phi', 'reflag', 'norm']

    # Create the first set of params
    parsNums = [-1]*len(params)
    pars = []
    for i, param in enumerate(params):
        pars.append('0.01 -1 -1e99 -1e99 1e99 1e99')

    # Load Model
    try:
        xs.AllModels.lmod(model_variant,PATHTOMODEL+'/'+model_variant+'/')
        m = xs.Model(model_variant)
        m.setPars(pars)
        print(f'\nModel {model_variant} succesfully loaded.')
    except:
        print(f'\n\n   ERROR: Could not load {model_variant} model. Exiting... \n')
        exit(1)

    _ = input('\nReady to print best-fitting model scheme. Press any key to continue... ')

    print()
    lines = chain[1].header['COMMENT']
    start_index = None
    end_index = None
    for i, line in enumerate(lines):
        if "Current model list:" in line.strip():
            if start_index is None:
                start_index = i
        if "_______________________________________" in line.strip() :
                end_index = i

    if start_index is not None and end_index is not None and start_index < end_index:
        for line in lines[start_index : end_index + 1]:
            print(line)
    print()

    print('Column numbers and parameters in the Chain FITS file:\n')
    num_columns = 3
    num_rows = (len(titles) + num_columns - 1) // num_columns
    for row in range(num_rows):
        for col in range(num_columns):
            index = row + col * num_rows
            if index < len(titles):
                print(f"{index + 1}. {titles[index]:<20}", end="")
        print()
    print()

    print('\nReady to assign columns and frozen/thaw values to model parameters.')
    print('   For "free" parameters, enter corresponding column number in Chain as an INTEGER.')
    print('   For "frozen" or "thaw" parameters, enter the parameter value as a FLOAT.')
    print()

    _ = input(f'Press any key to continue... ')
    print()

    for i, param in enumerate(params):
        flag = True
        while flag:
            p = input(f' Parameter: "{param}". Enter column (INT) or value (FLOAT): ')
            try:
                parNum = int(p)
                parsNums[i] = parNum
                flag = False
            except:
                try:
                    frozen = float(p)
                    pars[i] = frozen
                    parsNums[i] = -1
                    flag = False
                except:
                    flag = True

    m.setPars(pars)
    m.show()

    print('Calculating ETA_INT values...')
    if model_variant[:3] == 'vko':
        eta_ints = np.zeros(len(df))
        for i in tqdm(range(len(df))):
            for j, parsNum in enumerate(parsNums[:-1]):
                if parsNum == -1:
                    continue
                else:
                    pars[j] = df.iloc[i, parsNum-1]

            m.setPars(pars)
            m.show()
            xset_dict = dict(xs.Xset.modelStrings)
            eta_ints[i] = float(xset_dict['ETA_INT'])

    else:
        eta_ints1 = np.zeros(len(df))
        eta_ints2 = np.zeros(len(df))
        for i in tqdm(range(len(df))):
            for j, parsNum in enumerate(parsNums[:-1]):
                if parsNum == -1:
                    continue
                else:
                    pars[j] = df.iloc[i, parsNum-1]

            m.setPars(pars)
            m.show()
            xset_dict = dict(xs.Xset.modelStrings)
            eta_ints1[i] = float(xset_dict['ETA_INT1'])
            eta_ints2[i] = float(xset_dict['ETA_INT2'])

    quants = [0.05,0.16,0.50,0.84,0.95]
    if model_variant[:3] == 'vko':
        quantiles = np.quantile(eta_ints, quants)
        median = np.quantile(eta_ints, 0.5)
        print('\n\n ETA_INT QUANTILES:')
        for i, quant in enumerate(quants):
            print('  {:.3f}     {:.5f}     {:+.5f}'.format(quants[i], quantiles[i], quantiles[i]-median))
        print('\n')

    else:
        quantiles1 = np.quantile(eta_ints1, quants)
        quantiles2 = np.quantile(eta_ints2, quants)
        median1, median2 = np.quantile(eta_ints1, 0.5), np.quantile(eta_ints2, 0.5)
        print('\n\n ETA_INT1 QUANTILES:')
        for i, quant in enumerate(quants):
            print('  {:.3f}     {:.5f}     {:+.5f}'.format(quants[i], quantiles1[i], quantiles1[i]-median1))
        print('\n')
        print('\n ETA_INT2 QUANTILES:')
        for i, quant in enumerate(quants):
            print('  {:.3f}     {:.5f}     {:+.5f}'.format(quants[i], quantiles2[i], quantiles2[i]-median2))
        print('\n')

    print('\nProgram finished succesfully.\n')
    exit(0)
